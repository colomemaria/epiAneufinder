import sys
import seaborn as sns
import pandas as pd
import csv
import numpy as np
import scipy.stats
import fargv
from collections import defaultdict
from matplotlib import pyplot as plt
from sklearn.metrics import r2_score, mean_absolute_error, mutual_info_score, normalized_mutual_info_score
#plt.style.use('seaborn-whitegrid')
#sns.set_theme()


#Script for comparing the scATAC with the WGS dataset
#Inputs are the bed file containing the normalized count data from WGS data and the table with the results from epiAneufinder.
#The bed entries are intersect with the epiAneufinder windows and only common windows are kept.

standarize = lambda x: ((x - x.mean()) / (x.std() + .0000000000000001)) #The standarization formula for bringing the two datasets in the same scale
normalize = lambda x: ((2*(x-x.min())/(x.max()-x.min()))+1)


def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w


#create dictionary from bed file
def createDictionaryFromBed(bedfile):
    bed_dict=defaultdict(list)
    lines=bedfile.readlines()
    l=[]
    for line in lines:
        #print(line)
        if line.strip()[0] != "t":
            #chrom = "chr"+str(int(line.strip().split("\t")[0]))
            chrom=int(line.strip().split("\t")[0])
            start = int(line.strip().split("\t")[1])
            end = int(line.strip().split("\t")[2])
            somy = float(line.strip().split("\t")[3])
            if somy >=20:
                new_somy=20
            else:
                new_somy=somy

            l.append([chrom, start, end, new_somy])
            #d['chr1', 100000, 200000][0:10].count(1)

    for chrom, start, end,somy in l:
        bed_dict[chrom,start,end].append(somy)
        #print(bed_dict)
    #print((bed_dict))
    return(bed_dict)



#Function for creating a dictionary from the epiAneufinder data
def createDictionaryFromTable(table):
    snu_dict=table.set_index(['seq', 'start', 'end']).T.to_dict('list')
    return(snu_dict)


def calculatePopulationSomiesWGS(atac_dict,wgs_dict):
    #new_atac_dict = {k: [sum(atac_dict[k])/len(atac_dict[k])] for k in atac_dict.keys() if k[0]==4}
    new_atac_dict = {k: [sum(atac_dict[k]) / len(atac_dict[k])] for k in atac_dict.keys()}
    common_keys = set(wgs_dict).intersection(new_atac_dict) #filtering for the common CNV locations between the two datasets
    sort_common_keys=sorted(common_keys)
    #print(sort_common_keys)
    common_wgs_dict = {key: wgs_dict[key] for key in sort_common_keys}
    common_atac_dict = {key: new_atac_dict[key] for key in sort_common_keys}
    #with open('pGBM_keys.csv', 'wb') as f:
    #    w = csv.writer(f)
    #    w.writerows(" ".join(sort_common_keys))
    return(common_wgs_dict, common_atac_dict)

#Function for calculating different metrics between the two datasets and creating a line plot of the pseudobulk data
def createLinePlotAneufinder(common_wgs_dict, common_atac_dict,gaussian_sigma=0,filter_edges="constant"):
    wgs_list=([i for i in common_wgs_dict.values()])
    atac_list= [i for i in common_atac_dict.values()]
    wgs_array = np.hstack(wgs_list)
    smoothed_wgs_array = scipy.ndimage.gaussian_filter1d(wgs_array, gaussian_sigma, axis=- 1, order=0, output=None,
                                                   mode=filter_edges)
    atac_array = np.hstack(atac_list)
    print("Pearson Correlation : ", scipy.stats.pearsonr(standarize(atac_array),standarize(smoothed_wgs_array)))
    print("Spearman Correlation : ", scipy.stats.spearmanr(atac_array, smoothed_wgs_array)[0])
    print("Kendall Correlation : ", scipy.stats.kendalltau(atac_array, smoothed_wgs_array)[0])

    difference_array = np.subtract(atac_array, smoothed_wgs_array)
    squared_array = np.square(difference_array)
    mse = squared_array.mean()
    print("Mean Square Error: ",squared_array.mean())

    print("R square: ",r2_score(atac_array, smoothed_wgs_array))
    print("Mean Absolute Error: ",mean_absolute_error(atac_array,smoothed_wgs_array))
    x = list(range(len(smoothed_wgs_array)))

    plt.plot(x, (standarize(smoothed_wgs_array)), color='#98d1d1', label="WGS")
    #Borders between chromosomes for the plotting. One border list per sample
    borders_colorep1 = [0, 2109, 4392, 6241, 8068,9661, 11085, 12484, 13785, 14854, 16070, 17287, 18527, 19236, 20058, 20758, 21504, 22212, 22836, 23312, 23867, 24172, 24482]
    borders_pGBM2937=[0,2218,4507, 6361, 8204, 9908, 11561,13069, 14459, 15540, 16829, 18115, 19394, 20343, 21188,21940, 22691, 23448, 24178, 24713, 25308, 25637, 25963]
    borders_pGBM2932=[0, 2225,4519, 6373, 8216, 9920, 11572, 13086, 14476, 15558, 16850, 18136, 19417, 20367, 21212, 21967, 22722,23481, 24212, 24747, 25342, 25673, 26001]
    borders_pGBM4021=[0, 2215,4502, 6356, 8198, 9900, 11551, 13050, 14440, 15518, 16803, 18087, 19367, 20316, 21160, 21912, 22658, 23417, 24142, 24676, 25271, 25598, 25922]
    borders_pGBM3749=[0, 1961, 3911, 5560, 7034, 8468, 9881, 11199, 12527, 13437, 14521, 15597, 16676, 17431, 18123, 18771, 19386, 20117, 20744, 21238, 21757, 22064, 22365]
    borders_pGBM3402=[0,2198,4484, 6338, 8179, 9871, 11523,13019,14406, 15485, 16756, 18019, 19296, 20243, 21079, 21829, 22577, 23301, 24027, 24527, 25122, 25450, 25802]
    plt.plot(x,((atac_array)), color='#df979e', label="ATAC")
    plt.plot(x,((wgs_array)), color='#98d1d1', label="WGS")
    plt.title("pGBM sample 3402")
    for border in borders_pGBM3402:
        plt.axvline(border, color='gray')
    plt.xlim((0, len(atac_array)))
    plt.ylim(0.3,3)
    plt.legend()
    plt.show()
    return(scipy.stats.pearsonr(standarize(atac_array),standarize(smoothed_wgs_array)), smoothed_wgs_array, standarize(atac_array), standarize(smoothed_wgs_array))

if __name__ =="__main__":
    p = {"filter_edges": "nearest",
         "atac_input": "/home/katia/Helmholz/epiAneufinder/revisions/GSM4861367_COLO320HSR_rep1_atac_v3blacklist/epiAneufinder_results/results_table_noChr.tsv",
         "wgs_reads": "/home/katia/Helmholz/epiAneufinder/Colo320HSP/aneufinder/COLO320_HSR_WGS.medianNorm.bed",
         "smooth_sigma": .1,
         "title": "WGS vs scATAC"}
    p, _ = fargv.fargv(p)
    harmonics = {}
    fin = open(p.wgs_reads)
    snu_full=pd.read_csv(p.atac_input, sep=" ")
    snu_dict=createDictionaryFromTable(snu_full)
    bed_dict=createDictionaryFromBed(fin)
    common_wgs_dict, common_atac_dict = calculatePopulationSomiesWGS(snu_dict,bed_dict)
    
    #for sigma in range(20, 80, 5):
    #    p.smooth_sigma = sigma
    for _ in (0,):
        correlation, smoothed_wgs, standardATAC, standardWGS=createLinePlotAneufinder(common_wgs_dict, common_atac_dict,gaussian_sigma=p.smooth_sigma, filter_edges=p.filter_edges)
        harmonics[p.smooth_sigma] = correlation
        print(harmonics)

