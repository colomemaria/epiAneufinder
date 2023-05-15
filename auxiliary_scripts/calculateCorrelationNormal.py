import sys
import seaborn as sns
import pandas as pd
import numpy as np
import copy
import scipy.stats
from matplotlib import pyplot as plt
from sklearn.metrics import r2_score, mean_absolute_error
#plt.style.use('seaborn-whitegrid')
#sns.set_theme()


#Script for calculating the correlation between a normal sample and an ideal diploid sample.
#Input is the table with the results from epiAneufinder.
#Additional outputs is the MSE and Variation score.

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

#Function for creating a dictionary from the epiAneufinder data
def createDictionaryFromTable(table):
    snu_dict=table.set_index(['seq', 'start', 'end']).T.to_dict('list')
    #print(snu_dict)
    return(snu_dict)

#Function for creating the ideal normal sample, by copying the dictionary created from the epiAneufinder data
def createNormalBase(dictionary):
    base_dict=copy.deepcopy(dictionary)
    for i in base_dict.values():
        for j in range(len(i)):
            if i[j] != "2":
                i[j]="2"
    #print(base_dict)
    return(base_dict)

def calculatePopulationSomies(atac_dict,base_dict):
    base_wgs = []
    gain_atac = []
    loss_atac = []
    base_atac = []
    common_keys = set(base_dict).intersection(atac_dict) #filtering for the common CNV locations between the two datasets
    sort_common_keys=sorted(common_keys)
    #print(sort_common_keys)
    counts=0
    for k in sort_common_keys:
        if k[0]!=0: #selecting for all chromosomes
        #if k[0]!=0:  # selecting for all chromosomes
            counts=counts+1
            #Calculating pseudobulk representation for the scWGS. 1 is loss, 2 is disomic and 3 is gain
            base_wgs.append(base_dict[k].count("2") / len(base_dict[k]))
            #Calculating pseudobulk representation for the scATAC. 0 is loss, 1 is disomic and 2 is gain
            #If the user changes notation it should be changed here as well
            loss_atac.append(atac_dict[k].count(0) / len(atac_dict[k]))
            base_atac.append(atac_dict[k].count(1) / len(atac_dict[k]))
            gain_atac.append(atac_dict[k].count(2) / len(atac_dict[k]))
    print("Count Bins:",counts)
    return(base_wgs, loss_atac, base_atac, gain_atac)

def createLinePlot(base_wgs, loss_atac, base_atac, gain_atac):
    new_base_wgs = [x * 2 for x in base_wgs]
    #print(new_base_wgs)
    new_base_atac = [x * 2 for x in base_atac]
    new_gain_atac = [x * 3 for x in gain_atac]
    wgs_plot=[sum(x) for x in zip(new_base_wgs)]
    atac_plot = [sum(x) for x in zip(new_gain_atac, new_base_atac, loss_atac)]
    atac_array=np.array(atac_plot)
    wgs_array=np.array(wgs_plot)

    difference_array = np.subtract(atac_array, wgs_array)
    squared_array = np.square(difference_array)
    mse = squared_array.mean()
    print("Variation: ", scipy.stats.variation(atac_array))
    print("Mean Square Error: ",squared_array.mean())

    x = list(range(len(wgs_plot)))
    #Borders between chromosomes for the plotting. One border list per sample
    borders_bone_marrow=[0,1801, 3655, 5078, 6288, 7528, 8746, 9937, 10999, 11851, 12863, 13813, 14842, 15327, 15909, 16485, 17014, 17634, 18146, 18627, 19088, 19336, 19609]
    borders_PBMC1=[0, 1325, 2530, 3521, 4239, 5056, 5895, 6652, 7333, 7938, 8647, 9324, 10031, 10331, 10776, 11217,11643, 12194, 12482, 12926, 13293, 13469, 13692]
    brain_rep1=[0, 2107, 4361, 6223, 7981, 9631, 11225, 12674, 14024, 15086, 16296, 17465, 18690, 19590, 20419, 21163, 21875, 22610, 23329, 23831, 24400, 240707, 25024]
    brain_rep2=[0, 2076, 4311, 6157, 7875, 9487, 11064, 12484, 13820, 14857, 16049, 17205, 18416, 19298, 20109, 20851, 21554, 22278, 22980, 23481, 24046, 24345, 24658]
    brain_rep3=[0, 2076, 4315, 6168, 7894, 9523, 11094, 12527, 13868, 14903, 16100, 17268, 18488, 19379, 20198, 20938, 21642, 22368, 23076, 23576, 24145, 24451, 24763]
    brain_rep4=[0, 2078, 4305, 6151, 7871, 9485, 11051, 12464, 13804, 14846, 16036, 17193, 18410, 19297, 20115, 20856, 21562, 22292, 23001, 23502, 24069, 24379, 24690]
    brain_multi_rep1=[0, 1940, 4008,5670, 7095, 8538, 9943, 11240,12460, 13410, 14520, 15574, 16685, 17443, 18186, 18888, 19567, 20273, 20928, 21420, 21970, 22236, 22544]
    brain_multi_rep2 = [0, 1974, 4078, 5766, 7268, 8745, 10175, 11490, 12732, 13693, 14817, 15889, 17025, 17813, 18578, 19290, 19978, 20696, 21362, 21861, 22422, 22699, 23013]
    brain_multi_rep3 = [0, 1935, 3969, 5613, 7053, 8469, 9878, 11160, 12361, 13296, 14401, 15453, 16553, 17314, 18055, 18752, 19431, 20135, 20784, 21278, 21820, 22079, 22389]
    plt.plot(x, wgs_plot, color='blue', label="Base", linewidth='2')
    plt.plot(x, atac_plot, color='orange', label="ATAC")
    for border in brain_multi_rep3:
        plt.axvline(border, color='gray')
    plt.title("Brain rep4 normal dataset")
    plt.xlim((0, len(atac_plot)))
    plt.ylim(1,3)
    plt.legend()
    plt.show()

if __name__ =="__main__":
    snu_full=pd.read_csv("/home/katia/Helmholz/epiAneufinder/revisions/multiome_brain_Greenleaf/hft_ctx_w21_dc1r3_r1_ms0/epiAneufinder_results/multiome_brain_dc1r3_r1_msCNV0_results_table_noChr.tsv", sep=" ")
    snu_dict = createDictionaryFromTable(snu_full)
    base_dict = createNormalBase(snu_dict)
    #print(base_dict)
    base_wgs, loss_atac, base_atac, gain_atac = calculatePopulationSomies(snu_dict,base_dict)
    createLinePlot(base_wgs, loss_atac, base_atac, gain_atac)
