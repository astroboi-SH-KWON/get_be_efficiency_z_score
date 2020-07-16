from Bio import SeqIO

import Logic
class LogicPreps:
    def __init__(self, init):
        self.init = init

    def get_longer_seq_by_winsize(self, abe_list):
        chr_dir = self.init[0]
        win_size = self.init[1]
        chr_arr = self.init[2]
        result_dict = {}
        for chr_id in chr_arr:
            for seq_record in SeqIO.parse(chr_dir + "chr" + chr_id + ".fa", "fasta"):
                chr_seq = seq_record.seq.upper()
                rev_chr_seq = seq_record.seq.reverse_complement().upper()
                for val_arr in abe_list:
                    be_seq = val_arr[0].upper()
                    if be_seq in chr_seq:
                        str_chr_seq = str(chr_seq)
                        split_chr_seq = str_chr_seq[
                                        str_chr_seq.index(be_seq) - win_size[0]: str_chr_seq.index(be_seq) + len(
                                            be_seq) + win_size[1]]
                        if be_seq in result_dict:
                            result_dict[be_seq].append([split_chr_seq, "chr" + chr_id + ".fa"])
                        else:
                            result_dict.update({be_seq: [[split_chr_seq, "chr" + chr_id + ".fa"]]})

                    if be_seq in rev_chr_seq:
                        str_chr_seq = str(rev_chr_seq)
                        split_chr_seq = str_chr_seq[
                                        str_chr_seq.index(be_seq) - win_size[0]: str_chr_seq.index(be_seq) + len(
                                            be_seq) + win_size[1]]
                        if be_seq in result_dict:
                            result_dict[be_seq].append([split_chr_seq, "chr" + chr_id + ".fa"])
                        else:
                            result_dict.update({be_seq: [[split_chr_seq, "chr" + chr_id + ".fa"]]})

        return result_dict