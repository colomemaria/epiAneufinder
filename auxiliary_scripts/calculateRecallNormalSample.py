import sys
import seaborn as sns
import pandas as pd
import numpy as np
import scipy.stats
from collections import defaultdict
from matplotlib import pyplot as plt
import fargv

#Script for calculating the recall rate of a "normal" sample

def countNplicitiesFromTable(table):
    #Function that takes as input a csv file/output from epiAneufinder and calculates the pseudobulk dataset
    atac_dict=table.set_index(['seq', 'start', 'end']).T.to_dict('list')
    loss_atac=[]
    base_atac=[]
    gain_atac=[]
    for k in atac_dict:
        (chr, start, end) = k
        loss_atac.append(atac_dict[k].count(0))
        base_atac.append(atac_dict[k].count(1))
        gain_atac.append(atac_dict[k].count(2))
    return(loss_atac,base_atac,gain_atac)

def calculateRecall(loss_atac,base_atac,gain_atac):
    total_loss=np.sum(loss_atac)
    total_gain = np.sum(gain_atac)
    total_base = np.sum(base_atac)
    recall=total_base/(total_gain+total_loss+total_base)
    return(recall)




if __name__ =="__main__":
    "  "
    p = {"atac_input":"/home/katia/Helmholz/epiAneufinder/revisions/Satpathy_BoneMarrow_msCNV0/epiAneufinder_results/results_table.tsv"}
    p,_ = fargv.fargv(p)
    atac_loss, atac_base, atac_gain = countNplicitiesFromTable(pd.read_csv(p.atac_input, sep=" "))
    recall=calculateRecall(atac_loss, atac_base, atac_gain)
    print(recall)
