#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

"""
Convert a SAM file into GFF3 format (https://uswest.ensembl.org/info/website/upload/gff3.html)
 that also includes additional information that GMAP GFF3 for the 'mRNA' type:
   --- coverage
   --- identity
   --- matches
   --- mismatches
   --- indels


ex: from GMAP

6   cow_hereford    mRNA    109018744   109018861   .   +   .   \
     ID=myid.mrna1;Name=myid;Parent=myid.path1;coverage=29.5;identity=98.3;matches=116;mismatches=2;indels=0;unknowns=0

"""

#!/usr/bin/env python
import os, sys
import subprocess
from collections import Counter
from math import floor
from BCBio import GFF as BCBio_GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from cupcake.io.BioReaders import  GMAPSAMReader




class GFF3Writer_fromSAM:
    """Write GFF3 files starting with standard Biopython objects.
    """
    def __init__(self):
        pass

    # not sure how useful for this case?
    def _clean_feature(self, feature):
        quals = {}
        for key, val in feature.qualifiers.items():
            if not isinstance(val, (list, tuple)):
                val = [val]
            val = [str(x) for x in val]
            quals[key] = val
        feature.qualifiers = quals
        # Support for Biopython 1.68 and above, which removed sub_features
        if not hasattr(feature, "sub_features"):
            feature.sub_features = []
        clean_sub = [self._clean_feature(f) for f in feature.sub_features]
        feature.sub_features = clean_sub
        return feature

    def _get_phase(self, feature):
        if "phase" in feature.qualifiers:
            phase = feature.qualifiers["phase"][0]
        elif feature.type == "CDS":
            phase = int(feature.qualifiers.get("codon_start", [1])[0]) - 1
        else:
            phase = "."
        return str(phase)


    def _write_feature(self, feature, rec_id, out_handle, parent_id=None):
        """Write a feature with location information.
        """
        if feature.strand == 1:
            strand = '+'
        elif feature.strand == -1:
            strand = '-'
        else:
            strand = '.'

        if parent_id:
            attributes = parent_id
        else:
            attributes = 'transcript_id "' + feature.qualifiers.get("ID")[0] + '";'
        parent_id = attributes

        if feature.type:
            if feature.type == "gene":
                ftype = "transcript"
                attributes += ' gene_id "' + feature.qualifiers.get("ID")[0] + '"'
            else:
                ftype = feature.type # "exon" for this application
        else:
            ftype = "sequence_feature"

        parts = [str(rec_id),
                 feature.qualifiers.get("source", ["feature"])[0],
                 ftype,
                 str(feature.location.nofuzzy_start + 1), # 1-based indexing
                 str(feature.location.nofuzzy_end),
                 feature.qualifiers.get("score", ["."])[0],
                 strand,
                 self._get_phase(feature),
                 attributes]

        out_handle.write("\t".join(parts) + "\n")
        for sub_feature in feature.sub_features:
            self._write_feature(sub_feature, rec_id, out_handle, parent_id)
        return


    def write(self, rec, out_handle):
        for sf in rec.features:
            sf = self._clean_feature(sf)
            self._write_feature(sf, rec.id, out_handle)



def convert_sam_rec_to_gff3_rec(r, source, qid_index_dict=None, allow_non_primary=True):
    """
    :param r: GMAPSAMRecord record
	:param qid_seen: list of qIDs processed so far -- if redundant, we have to put a unique suffix
    :return SeqRecord ready to be written as GFF3
    """
    if r.sID == '*':
        print("Skipping {0} because unmapped.".format(r.qID), file=sys.stderr)
        return None
    t_len = sum(e.end-e.start for e in r.segments)
 
    if not allow_non_primary and r.flag.is_secondary:
        return None

    seq = Seq('A'*t_len)  # DO NOT CARE since sequence is not written in GFF3
    rec = SeqRecord(seq, r.sID)
    strand = 1 if r.flag.strand == '+' else -1

    indels = r.num_ins+r.num_del
    mismatches = r.num_nonmatches
    matches = r.num_mat_or_sub - r.num_nonmatches

 
    if qid_index_dict is not None and allow_non_primary:
        qid_index_dict[r.qID] += 1
        if r.flag.is_supplementary:
            r.qID += '_sup' + str(qid_index_dict[r.qID])
        elif r.flag.is_secondary:
            r.qID += '_sec' + str(qid_index_dict[r.qID])
        else:
            r.qID += '_unk' + str(qid_index_dict[r.qID])


    gene_qualifiers = {"source": source, "ID": r.qID, "Name": r.qID} # for gene record
#    mRNA_qualifiers = {"source": source, "ID": r.qID+'.mRNA', "Name": r.qID+'.mRNA', "Parent": r.qID,
#                       "coverage": "{0:.2f}".format(r.qCoverage*10**2) if r.qCoverage is not None else "NA",
#                       "identity": "{0:.2f}".format(r.identity*10**2),
#                       "matches": matches, "mismatches": mismatches, "indels": indels}

    # gene line, one per record
    top_feature = SeqFeature(FeatureLocation(r.sStart, r.sEnd), type="gene", strand=strand, qualifiers=gene_qualifiers)
    # mRNA line, one per record
    top_feature.sub_features = [] #top_feature.sub_features = [SeqFeature(FeatureLocation(r.sStart, r.sEnd), type="mRNA", strand=strand, qualifiers=mRNA_qualifiers)]

    # exon lines, as many exons per record
    for i,e in enumerate(r.segments):
        _id = "{0}.exon{1}".format(r.qID,i+1)
        exon_qual = {"source": source, "ID": _id, "Name": _id}
        top_feature.sub_features.append(SeqFeature(FeatureLocation(e.start, e.end), type="exon", strand=strand, qualifiers=exon_qual))
    rec.features = [top_feature]
    return rec

def convert_sam_to_gff3(sam_filename, output_gff3, source, q_dict=None, allow_non_primary=True):
    qid_index_dict = Counter()
    writer = GFF3Writer_fromSAM()
    with open(output_gff3, 'w') as f:
        for r0 in GMAPSAMReader(sam_filename, True, query_len_dict=q_dict):
            rec = convert_sam_rec_to_gff3_rec(r0, source, qid_index_dict, allow_non_primary)
            if rec is not None:
                writer.write(rec, f)

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser("Convert SAM to GFF3 format using BCBio GFF")
    parser.add_argument("sam_filename")
    parser.add_argument("-i", "--input_fasta", default=None, help="(Optional) input fasta. If given, coverage will be calculated.")
    parser.add_argument("-s", "--source", required=True, help="source name (ex: hg38, mm10)")

    args = parser.parse_args()

    if not args.sam_filename.endswith('.sam'):
        print("Only accepts files ending in .sam. Abort!", file=sys.stderr)
        sys.exit(-1)

    prefix = args.sam_filename[:-4]
    output_gff3 = prefix + '.gff3'

    q_dict = None
    if args.input_fasta is not None:
        q_dict = dict((r.id, len(r.seq)) for r in SeqIO.parse(open(args.input_fasta), 'fasta'))

    with open(output_gff3, 'w') as f:
        recs = [convert_sam_rec_to_gff3_rec(r0, args.source) for r0 in GMAPSAMReader(args.sam_filename, True, query_len_dict=q_dict)]
        BCBio_GFF.write([x for x in recs if x is not None], f)


    print("Output written to {0}.".format(output_gff3), file=sys.stderr)

if __name__ == "__main__":
    main()
