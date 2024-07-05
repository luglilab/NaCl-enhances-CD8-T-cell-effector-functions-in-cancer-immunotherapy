import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


# load data file
d = pd.read_csv("/Users/simonepuccio/Documents/HumanitasProjects/SP025/input_volcanoplot.txt",
                sep="\t",
                header=0)

d.loc[(d['log2FoldChange'] >= 0.5) & (d['padj'] < 0.05), 'color'] = "green"  # upregulated
d.loc[(d['log2FoldChange'] < -0.5) & (d['padj'] < 0.05), 'color'] = "black"   # downregulated
d['color'].fillna('lightgrey', inplace=True) # intermediate
d['logpv']=-(np.log10(d['padj']))



plt.text(-3.43566, 21, "Gzmk")
plt.text(-1.712629, 20.970616, "Klf2")
plt.text(-2.030216, 9.970616, "Eomes")
plt.text(-1.551007, 8.101824, "Cxcr3")
plt.text(-1.269935, 7.510042, "Ccl4")
plt.text(-1.207603, 7.498941, "Ccl3")
plt.text(-0.985961, 6.573489, "Tox")
plt.text(-1.932466, 5.446117, "Prdm1")
plt.text(-0.80915, 5.055024, "Ccr5")
plt.text(-0.768472, 3.133766, "Tnf")
plt.text(-1.271413, 1.484424, "Maf")


plt.text(2.118559, 27.575118, "Klrd1")
plt.text(1.380999, 23.261219, "Irf8")
plt.text(1.033581, 19.527244, "Stat3")
plt.text(1.53703, 18.543634, "Klrc1")
plt.text(1.794269, 17.318759, "Prf1")
plt.text(1.961144, 15.844664, "Gzmc")
plt.text(0.681256, 15.806875, "Il2rb")
plt.text(0.992855, 15.108463, "Id2")
plt.text(2.332166, 14.815309, "Il7r")
plt.text(4.656493, 13.531653, "Gzmf")
plt.text(1.124743, 9.485452, "Satb1")
plt.text(2.422941, 7.777284, "Csf1")
plt.text(3.679162, 5.478862, "Gzme")
plt.text(1.764106, 2.420362, "Cx3cr1")
plt.text(0.283454, 1.164489, "Ifngr1")
plt.text(0.730974, 4.856985, "Runx2")




plt.scatter(d['log2FoldChange'], d['logpv'], c=d['color'],s=10,
            alpha=1)
plt.xlim(-7, 7)
plt.ylim(0, 30)
plt.xlabel('log2 Fold Change',fontsize=15, fontname="sans-serif", fontweight="bold")
plt.ylabel('-log10(Q-value)', fontsize=15, fontname="sans-serif", fontweight="bold")
plt.xticks(fontsize=12, fontname="sans-serif")
plt.yticks(fontsize=12, fontname="sans-serif")
plt.axhline(1.3010, color='black', lw=1,linestyle="dashed")
plt.axvline(0.5, color='black', lw=1,linestyle="dashed")
plt.axvline(-0.5, color='black', lw=1,linestyle="dashed")
# plt.show()
plt.savefig('/Users/simonepuccio/Documents/HumanitasProjects/SP025/Mouse_volcanoplotCD8.pdf', format='pdf', bbox_inches='tight', dpi=300)

