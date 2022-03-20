def parse(filename):
    # read [.gfa] file, just get abundances for now
    #S       840     *       LN:i:1  KC:i:1
    kmer_abundance = dict()
    with open(filename) as f:
        for line in f:
            if line.startswith('S'):
                spl = line.split()
                seq_id = spl[1]
                for field in spl[2:]:
                    if field.startswith("KC:i"):
                        abundance = int(field.split(':')[-1])
                        kmer_abundance[seq_id] = abundance

    return kmer_abundance
