def read_fasta(fp): # this is a generator function
    "See https://github.com/lh3/readfq"
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng = ''.join(seqs), 0
            for l in fp: # read the quality
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq # yield a fasta record instead
                break
