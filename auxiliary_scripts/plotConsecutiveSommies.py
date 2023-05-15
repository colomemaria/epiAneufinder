import matplotlib.pyplot as plt
import pandas as pd
from itertools import groupby
import numpy as np

#Script for plotting the consecutive somies per sample.
#One figure for gain/loss/disomic regions
#Used to calculate the maximum and minimum CNV size

wd="/home/katia/Helmholz/epiAneufinder/revisions/SNU601_br15_minsizeCNV0/epiAneufinder_results/"
temp=pd.read_csv(wd+"results_table.tsv", sep=" ")
df = temp.iloc[: , 3:]
count0=[]
count1=[]
count2=[]
for i in range(len(df.columns)):
    groups = groupby(df.iloc[:, i])
    result = ([(label, sum(1 for _ in group)) for label, group in groups])
    for i in result:
        if i[0] == 1:
            count1.append(i[1])
        if i[0] == 2:
            count2.append(i[1])
        if i[0] == 0:
            count0.append(i[1])
plt.hist(count0, bins=np.arange(min(count0), max(count0) + 10, 10))
plt.ylabel('Number')
plt.xlabel('Size')
plt.title("Consecutive bins with loss (min="+str(min(count0))+" ,max="+str(max(count0))+")")
plt.savefig(wd+'consecutiveLossHist.png')
plt.close()
plt.hist(count1, bins=np.arange(min(count1), max(count1) + 10, 10))
plt.ylabel('Number')
plt.xlabel('Size')
plt.title("Consecutive disomic bins (min="+str(min(count1))+" ,max="+str(max(count1))+")")
plt.savefig(wd+'/consecutiveDisomicHist.png')
plt.close()
plt.hist(count2, bins=np.arange(min(count2), max(count2) + 10, 10))
plt.ylabel('Number')
plt.xlabel('Size')
plt.title("Consecutive bins with gain (min="+str(min(count2))+" ,max="+str(max(count2))+")")
plt.savefig(wd+'consecutiveGainHist.png')
plt.close()
