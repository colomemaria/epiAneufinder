import pandas as pd
import sys
import fargv
import seaborn as sns
import matplotlib.pyplot as plt

#Script for calculating the full similarity percentage between the SNU601 and the different sub-samplings.
#Input are the table from epiAneufinder for the whole SNU601 dataset and the corresponding tables of the sub-sampligs.
#The user can decide if the disomic regions should be included in the calculation or not.
#The output is one file with the similarity percenatge per sub-sampling, that is used for the box-plot of figure 4a.

def remove_disomic(cnv_table):
    for col in cnv_table:
        cnv_table=cnv_table.loc[cnv_table[col]==2, col] = -10
        return(cnv_table)

def calculate_difference(full, percent):
    cell=[]
    for col in percent.columns:
        cell.append(col)
    full_cell_resize=full[full.columns.intersection(cell)]
    full_percent_resize = full_cell_resize[full_cell_resize.set_index(['seqnames','start','end']).index.isin(percent.set_index(['seqnames','start','end']).index)]
    full_percent_resize=full_percent_resize[percent.columns] #making sure the ordering of the columns is the saame
    full_percent_diff = (full_percent_resize.set_index(['seqnames','start','end'])- percent.set_index(['seqnames','start','end'])).reset_index()
    return(full_percent_diff,cell)

if __name__ =="__main__":
    params = {"input": "/home/katia/Helmholz/epiAneufinder/results/SNU601/final/SNU601_FINAL_full_finalTable.csv", "filepath":"/home/katia/Helmholz/epiAneufinder/results/subsampling/",
              "disomic":"False"}
    p, _ = fargv.fargv(params)
    if p.disomic=="True":
        full = pd.read_csv(p.input)
        print(full)
        a = 0.1
        b=1
        while a== 0.1:
            print("In the while disomic=True")
            percentage = pd.read_csv(p.filepath + str(a) + "/Final/SNU601_0" + str(b) + "_FINAL_finalTable.csv")
            print(percentage)
            for col in percentage:
                percentage.loc[percentage[col]==2, col] = "NA"
            print(percentage)
            for col in full:
                full.loc[full[col]==2, col] = "NA"
            print(full)
            #full=remove_disomic(full)
            full_percent_diff, cell = calculate_difference(full, percentage)
            percent_similarity = []
            for c in cell:
                if c != ("seqnames") and c != ("start") and c != ("end"):
                    percent_similarity.append(100 * ((full_percent_diff[c] == 0).sum() / len(full_percent_diff)))
                percent_similarity_df = pd.DataFrame(percent_similarity, columns=[str(a)])
                #print(percent_similarity_df)
            out = "table_FINAL_" + str(a) + "_similarity_nodisomic.csv"
            percent_similarity_df.to_csv(out, index=False, header=True, sep=",")
            a = round(a + 0.1, 1)
            b=b+1
            print("end iteration")
    else:
        full = pd.read_csv(p.input)
        a=0.1
        b=1
        while a<1:
            print("In the while disomic=flase")
            percentage = pd.read_csv(p.filepath + str(a) + "/Final/SNU601_0" + str(b) + "_FINAL_finalTable.csv")
            print(p.filepath + str(a) + "/Final/SNU601_0" + str(b) + "_FINAL_finalTable.csv")
            full_percent_diff, cell=calculate_difference(full,percentage)
            percent_similarity = []
            for c in cell:
                if c != ("seqnames") and c != ("start") and c != ("end"):
                    percent_similarity.append(100 * ((full_percent_diff[c] == 0).sum() / len(full_percent_diff)))
                percent_similarity_df = pd.DataFrame(percent_similarity,columns=[str(a)])
                print(percent_similarity_df)
            out = "table_FINAL_"+str(a)+"_similarity.csv"
            percent_similarity_df.to_csv(out, index=False, header=True, sep=",")
            a=round(a+0.1,1)
            b=b+1
            print("end iteration")


