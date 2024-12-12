#!/usr/bin/env python

import sys, re
from Bio import SeqIO
from collections import defaultdict

# def import_repeat_finder(infile):
#     imprt = False
#     repeat_d = {}
#     for i, line in enumerate(infile):
#         if i == 7:
#             imprt = True
#         if line.strip() == "" and imprt == True:
#             break
#         if imprt == True:
#             data = line.strip().split()
#             seq = data[4]
#             hits = data[1]
#             repeat_d[seq] = hits
#     infile.close()
#     return repeat_d

# def get_exact_matches(repeats, DUSs):
#     out_d = []
#     for seq, hits in repeats.items():
#         for k,v in DUSs.items():
#             if seq == v:
#                 out_d.append([k, hits, seq, v])
#     return out_d

def generate_re_term(seq, i, base):
    if base == 'A':
        re_term =  '{}[TGC]{}'.format(seq[:i], seq[i+1:])
    if base == 'G':
        re_term =  '{}[TAC]{}'.format(seq[:i], seq[i+1:])
    if base == 'T':
        re_term =  '{}[AGC]{}'.format(seq[:i], seq[i+1:])
    if base == 'C':
        re_term =  '{}[TGA]{}'.format(seq[:i], seq[i+1:])
    return re_term

def get_single_nuc_variants(recs, DUSs):
    hits_d = defaultdict(int)
    for rec in recs:
        for k,v in DUSs.items():
            hits_d[(k,v)] += len(re.findall(v, str(rec.seq)))
            hits_d[(k,v)] += len(re.findall(v, str(rec.seq.reverse_complement())))
            for i, base in enumerate(v):
                re_term =  generate_re_term(v, i, base)
                re_hits = re.findall(re_term, str(rec.seq)) + re.findall(re_term, str(rec.seq.reverse_complement()))
                for hit in re_hits:
                    if not hit in list(DUSs.values()):
                        hits_d[(f"{k}-like",hit)] += 1
    return hits_d

def check_if_single_nuc_variant_already_defined(DUSs, seq, v):
    for v2 in DUSs.values():
        if v2 == v:
            continue
        if v2 == seq:
            return True
    return False

def fix_double_hits(hits_d):
    fixed_hits_d = {}
    multiple_hits_d = defaultdict(list)
    for name, seq in hits_d.keys():
        multiple_hits_d[seq].append(name)
    for k, v in multiple_hits_d.items():
        if len(v) == 1:
            fixed_hits_d[(v[0], k)] = hits_d[(v[0],k)]
        else: 
            new_name_list = [ x.replace("-like", "") for x in v ]
            new_name = "-".join(new_name_list) + "-like"
            fixed_hits_d[(new_name, k)] = hits_d[(v[0],k)]
    return fixed_hits_d

def main():
    DUSs = {
        'AT-DUS':'ATGCCGTCTGAA',
        'TG-wadDUS':'TGCCTGTCTGAA',
        'AG-DUS':'AGGCCGTCTGAA',
        'AG-mucDUS':'AGGTCGTCTGAA',
        'AG-simDUS':'AGGCTGCCTGAA',
        'AG-kingDUS':'AGGCAGCCTGAA',
        'AG-king3DUS':'AAGCAGCCTGCA',
        'AG-eikDUS':'AGGCTACCTGAA'
        }
    recs = SeqIO.parse(sys.argv[1], 'fasta')
    #repeats = import_repeat_finder(infile)
    pre_hits_d = get_single_nuc_variants(recs, DUSs)
    hits_d = fix_double_hits(pre_hits_d)
    outfile = open(sys.argv[2], 'w')
    sorted_hits_d = sorted(hits_d.items(), key=lambda x:x[1], reverse=True)
    for k, v in sorted_hits_d:
        outfile.write(f"{k[0]},{k[1]},{v}\n")
    outfile.close()

    
if __name__ == '__main__':
    main()
