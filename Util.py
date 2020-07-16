import glob

import Logic

class Utils:
    def __init__(self):
        self.ext_txt = ".txt"
        self.ext_dat = ".dat"
        self.ext_xlsx = ".xlsx"

    def read_tb_txt(self, path):
        result_list = []
        with open(path, "r") as f:
            while True:
                tmp_line = f.readline().replace("\n", "")
                if tmp_line == "":
                    break

                tmp_arr = tmp_line.split("\t")
                result_list.append(tmp_arr)

        return result_list

    def make_tab_txt(self, path, data_dict):
        with open(path, "a") as f:
            idx = 1
            f.write("IDX\t30 bps = 4 bp + 20 bp protospacer + PAM(NGG) + 3 bp\t50 bps = 16 bp + (30 bps) + 4 bp\tchr\n")
            for be_seq, split_chr_seq_list in data_dict.items():
                for val_arr in split_chr_seq_list:
                    f.write(str(idx) + "\t" + be_seq + "\t" + val_arr[0] + "\t" + val_arr[1] + "\n")
                    idx += 1

    def make_tab_txt_exception_data(self, path, output_dict, input_list):
        with open(path, "a") as f:
            idx = 1
            f.write("IDX\tno match seq\n")
            for be_seq_arr in input_list:
                be_seq = be_seq_arr[0]
                if be_seq not in output_dict:
                    f.write(str(idx) + "\t" + be_seq + "\n")
                    idx += 1

    def make_tb_txt_w_efficiency_z_score(self, path, efficiency_z_score_list):
        with open(path, "a") as f:
            idx = 1
            f.write("IDX\t30 bps = 4 bp + 20 bp protospacer + PAM(NGG) + 3 bp\t50 bps = 16 bp + (30 bps) + 4 bp\tz_score\tchr\n")
            for be_seq_arr in efficiency_z_score_list:
                f.write(
                    str(idx) + "\t" + be_seq_arr[1] + "\t" + be_seq_arr[2] + "\t" + str(be_seq_arr[4]) + "\t" + be_seq_arr[
                        3] + "\n")
                idx += 1

