import sys
from typing import Tuple

from parse_sequences_file import parse


def find_overlap(source, sink):
    shift = source[3][0] if source[1] == "+" else source[3][1]
    assert(shift < len(source[2]))
    return shift


def gfa_segment(node, kmer_abundance) -> Tuple[str, str]:
    segment_id = node[0]
    segment_seq = node[2]
    segment_len = f"LN:i:{len(segment_seq)}"
    # The last element must always have newline char
    segment_kc = f"KC:i:{kmer_abundance[segment_id]}\n"

    return segment_id, "\t".join(["S", segment_id, segment_seq, segment_len, segment_kc])


def gfa_link(source, sink, len) -> Tuple[str, str]:
    link_from = source[0]
    link_fromorient = source[1]
    link_to = sink[0]
    link_toorient = sink[1]
    # The last element must always have newline char
    link_overlap = f"{len}M\n"

    link_id = "".join([link_from, link_fromorient, link_to, link_toorient])
    return link_id, "\t".join(["L", link_from, link_fromorient, link_to, link_toorient, link_overlap])


segments = {}
links = {}

k, l, node_minims, kmer_seq, kmer_abundance, origins1, minim_shift = parse(sys.argv[1])

def update_kmer_abundance(line, kmer_abundance):
    #S       217     *       LN:i:1  KC:i:137
    index = line.split()[1]
    for field in line.split():
        if field.startswith('KC'):
            abundance = int(field.split(':')[-1])
            kmer_abundance[index] = abundance


with open(sys.argv[2], "r") as gfa_file:
    for line in gfa_file:
        if line.startswith('S'):
            update_kmer_abundance(line, kmer_abundance)
            continue
        elif line.startswith('H'): continue
        elif not line.startswith('L'):
            print("weird GFA line: " + line)
            sys.exit(1)
    
        #L	2377	+	5976	+	0M
        spl = line.split()
        source_minims = node_minims[spl[1]]
        sink_minims = node_minims[spl[3]]
        source = (spl[1], spl[2], kmer_seq[source_minims], minim_shift[spl[1]])
        sink = (spl[3], spl[4], kmer_seq[sink_minims], minim_shift[spl[3]])
        shift = find_overlap(source, sink)
        overlap_length = len(source[2]) - shift
        overlap_length = min(overlap_length, len(sink[2])-1)

        source_id, segments[source_id] = gfa_segment(source, kmer_abundance)
        sink_id, segments[sink_id] = gfa_segment(sink, kmer_abundance)
        link_id, links[link_id] = gfa_link(source, sink, overlap_length)
    

output_filename = '.'.join(sys.argv[1].split('.')[:-1])+".complete.gfa"
print(f"Writing to {output_filename}")

with open(output_filename, 'w') as output:
    output.write("H\tVN:Z:1.0\n")
    output.writelines(segments.values())
    output.writelines(links.values())