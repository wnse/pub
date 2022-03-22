from Bio import SeqIO

def filter_fa_by_len(fa, outfile, cutoff=500, name='scaffolds'):
    filter_fa = []
    max_len = 0
    seq_count = 0
    name_no = 1
    try:
        for seq in SeqIO.parse(fa, 'fasta'):
            tmp_len = len(seq.seq)
            if tmp_len > max_len:
                max_len = tmp_len
            if tmp_len >cutoff:
                seq.id = name + '_' + str(name_no)
                seq.description = seq.id + ':' + str(len(seq.seq))
                filter_fa.append(seq)
                name_no += 1
        if filter_fa:
            seq_count = SeqIO.write(filter_fa, outfile, 'fasta')
    except Exception as e:
        logging.error(f'filter_fa_by_len {e}')
    return seq_count, max_len