import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import fargv
import sys
#sns.set_style(style="white")
sns.set(font_scale = 2)
sns.set_style("whitegrid", {'axes.grid' : False})

#Script for creating the boxplots of figure 4, using as inputs the results from the cnvcall_agreement.py script, for all subsamplings.
#The user will have to modify the paths of the files in order to run the script.

tmp_loss = []
tmp_gain = []
tmp_disomy = []
a = 0.1
b=1
while a<1:
    f_gain=open("/home/katia/Helmholz/epiAneufinder/results/subsampling/"+str(a)+"/Final/SNU601_0"+str(b)+"_FINAL_finalTablegains_lost.txt")
    f_loss=open("/home/katia/Helmholz/epiAneufinder/results/subsampling/"+str(a)+"/Final/SNU601_0"+str(b)+"_FINAL_finalTablelosses_lost.txt")
    f_disomy=open("/home/katia/Helmholz/epiAneufinder/results/subsampling/"+str(a)+"/Final/SNU601_0"+str(b)+"_FINAL_finalTabledisomy_lost.txt")
    l_lines=f_loss.readlines()
    g_lines = f_gain.readlines()
    d_lines = f_disomy.readlines()
    for line in l_lines:
        if line.strip() != "-10000.0":
            tmp_loss.append(line.strip() + '\t'+str(100*a)+"%")
    for line in g_lines:
        if line.strip() != "-10000.0":
            tmp_gain.append(line.strip() + '\t'+str(100*a)+"%")
    for line in d_lines:
        if line.strip() != "-10000.0":
            tmp_disomy.append(line.strip() + '\t'+str(100*a)+"%")
    a = round(a + 0.1, 1)
    b=b+1
fout_losses=open("/home/katia/Helmholz/epiAneufinder/results/subsampling/table_lost_losses_Final.txt","w")
fout_gains=open("/home/katia/Helmholz/epiAneufinder/results/subsampling/table_lost_gains_Final.txt","w")
fout_disomy=open("/home/katia/Helmholz/epiAneufinder/results/subsampling/table_lost_disomy_Final.txt","w")
fout_losses.write("\n".join(tmp_loss))
fout_gains.write("\n".join(tmp_gain))
fout_disomy.write("\n".join(tmp_disomy))

Gain = pd.read_csv("/home/katia/Helmholz/epiAneufinder/results/subsampling/table_lost_gains_Final.txt", sep='\t',
                  names=["Lost gains","Percent"])
Loss = pd.read_csv("/home/katia/Helmholz/epiAneufinder/results/subsampling/table_lost_losses_Final.txt", sep='\t',
                  names=["Lost losses","Percent"])
Disomy = pd.read_csv("/home/katia/Helmholz/epiAneufinder/results/subsampling/table_lost_disomy_Final.txt", sep='\t',
                  names=["Lost disomy","Percent"])
Gain = Gain[Gain.Percent != "10.0%"]

Loss = Loss[Loss.Percent != "10.0%"]
Disomy = Disomy[Disomy.Percent != "10.0%"]
ax1=sns.boxplot(y="Percent", x="Lost gains", data=Gain, showmeans=True,meanprops={"markeredgecolor": "yellow",
                       "markersize": "10"})
ax1.set_box_aspect(1)
print(Gain.groupby(["Percent"])["Lost gains"].describe())
plt.show()
ax2=sns.boxplot(y="Percent", x="Lost losses", data=Loss,showmeans=True,meanprops={"markeredgecolor": "yellow",
                       "markersize": "10"})
ax2.set_box_aspect(1)
print(Loss.groupby(["Percent"])["Lost losses"].describe())
plt.show()
ax3=sns.boxplot(x="Lost disomy", y="Percent", data=Disomy,showmeans=True,meanprops={"markeredgecolor": "yellow",
                       "markersize": "10"})
ax3.set_box_aspect(1)
print(Disomy.groupby(["Percent"])["Lost disomy"].describe())
plt.show()