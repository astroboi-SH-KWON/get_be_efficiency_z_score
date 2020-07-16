import time
import os
import numpy as np
from Bio import SeqIO

import predict as be_efficiency_model
import Util
import LogicPrep
import Logic
############### start to set env ###############
# WORK_DIR = os.getcwd() + "/"
WORK_DIR = "D:/000_WORK/KimNahye/20200622/WORK_DIR/"
CHR_DIR = "D:/000_WORK/000_reference_path/human/hg38/Splited/"

INPUT_DIR = "input/"
OUTPUT_DIR = "output/"
CHR_ARR = ["X", "Y"]

BE_ARR = ["ABE", "BE4"]
# BE_ARR = ["ABE_test"]

WIN_SIZE = [16, 4]

PATH_SPLIT = "/input"

INIT = [CHR_DIR, WIN_SIZE]
############### end setting env ################

def get_be_efficiency_z_score():
    util = Util.Utils()

    for be_nm in BE_ARR:
        be_list = util.read_tb_txt(WORK_DIR + INPUT_DIR + be_nm + "_long_seq_list.txt")

        be_efficiency_model.init_model(base_editor=be_nm, celltype='HEK293T')  # BE4

        for be_arr in be_list:
            trgt_seq = be_arr[2]

            pred_d = be_efficiency_model.predict(trgt_seq)
            be_arr.append(pred_d['Predicted logit score'])

        util.make_tb_txt_w_efficiency_z_score(WORK_DIR + OUTPUT_DIR + be_nm + "_efficiency_z_score.txt", be_list)

def get_human_seq_by_windows():
    for idx in range(1, 23):
        CHR_ARR.append(str(idx))

    INIT.append(CHR_ARR)
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps(INIT)

    for be_nm in BE_ARR:
        be_list = util.read_tb_txt(WORK_DIR + INPUT_DIR + be_nm + ".txt")

        long_seq_dict = logic_prep.get_longer_seq_by_winsize(be_list)

        util.make_tab_txt_exception_data(WORK_DIR + OUTPUT_DIR + be_nm + "_exception_list.txt", long_seq_dict, be_list)

        util.make_tab_txt(WORK_DIR + OUTPUT_DIR + be_nm + "_long_seq_list.txt", long_seq_dict)

def test2():
    for seq_record in SeqIO.parse("D:/000_WORK/000_reference_path/human/hg19/chromosomes/" + "chr" + "15" + ".fa", "fasta"):
        chr_seq = seq_record.seq.upper()
        rev_chr_seq = seq_record.seq.reverse_complement().upper()
        print(chr_seq)
        print(rev_chr_seq)
        if "GAACCGCCGTAACTGCCGCGCCCGGGGGAG" in chr_seq:
            print("True")
        if "GAACCGCCGTAACTGCCGCGCCCGGGGGAG" in rev_chr_seq:
            print("-")


start_time = time.perf_counter()
print("start >>>>>>>>>>>>>>>>>>")
get_be_efficiency_z_score()
# get_human_seq_by_windows()
# test2()
# print(1 / (1 + np.exp(-(0.29512018470961776 * 2 + 0.5))))
end_time = time.perf_counter()
print("::::::::::: %.2f seconds ::::::::::::::" % (end_time - start_time))