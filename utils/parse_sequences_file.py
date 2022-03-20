import glob
from typing import NamedTuple

import lz4.frame

class Sequences(NamedTuple):
    k: int
    l: int
    node_minims: dict
    kmer_to_seq: dict
    kmer_abundance: dict
    kmer_origins: dict
    minim_shift: dict


def parse(prefix: str):
    kmer_to_seq = dict()
    kmer_abundance = dict()
    node_minims = dict()
    k, l = 0, 0
    kmer_origins = dict()
    minim_shift = dict()
    for filename in glob.glob(prefix+"*.sequences"):
        for line in lz4.frame.open(filename):
            line = line.decode()
            # ignore #'s lines except for getting the k value
            if line.startswith('#'):
                if line.startswith('# k = '):
                    k = int(line.split()[-1])
                if line.startswith('# l = '):
                    l = int(line.split()[-1])
                continue
            spl = line.split()
            seq_id = spl[0]
            minims = tuple(map(lambda x: int(x.strip('[').strip(']').replace(',','')),spl[1:-5]))
            if spl[-4] != "PLACEHOLDER": 
                abundance = spl[-4]
                if abundance != "*": abundance = int(abundance)
                kmer_abundance[seq_id] = abundance
            origin = spl[-3]
            seq = spl[-5]
            shifts = tuple(map(lambda x: int(x.strip('(').strip(')').replace(',','')),spl[-2:]))
            node_minims[seq_id] = minims
            kmer_to_seq[minims] = seq
            kmer_origins[minims] = origin
            minim_shift[seq_id] = shifts

    return k, l, node_minims, kmer_to_seq, kmer_abundance, kmer_origins, minim_shift

# File fromat
# <line> ::= <hashtag_line> | <sequence_line> <line_end>
# <hashtag_line> ::= # <parameter_line> | # <comment_line>
# <comment_line> ::= <string>
# <parameter_line> ::= "k = <int>" | "l = <int>"
# <sequence_line> ::= <seq_id> <minim_list> <abundance> <origin> <seq> <shifts>
# <seq_id> ::= <string>
# <minim_list> ::= <item> | <item> <minim_list>
# <item> ::= <int> | <int> <comma>
# <abundance> ::= <int> | PLACEHOLDER | *
# <origin> ::= <int> | *
# <seq> ::= <string>
# <shifts> ::= <int> <comma> <int>
