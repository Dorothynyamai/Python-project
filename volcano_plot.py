#A script to plot volcano plots

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
import random

#read the csv file
deseq_results = pd.read_csv('ctnnb_TPMdeseq_analysis.csv').dropna()
deseq_results['nlog10'] = -np.log10(deseq_results.padj)



#Determine the upregulated and downregulated genes based on a specified threshold for log2 fold change and p-adj

upregulated_genes = deseq_results[(deseq_results['log2FoldChange'] > 0.5) & (deseq_results['padj'] < 0.05)].symbol.tolist()

upcount = len(upregulated_genes)

#print("Number of upregulated genes", upcount)
#not_diff = deseq_results[deseq_results['padj'] > 0.05]

downregulated_genes = deseq_results[(deseq_results['log2FoldChange'] < -0.5) & (deseq_results['padj'] < 0.05)].symbol.tolist()

downcount = len(downregulated_genes)

#print("Number of downregulated genes", downcount)

notdiff_genes = deseq_results[((deseq_results['log2FoldChange'] < 0.5) & (deseq_results['log2FoldChange'] > -0.5) & (deseq_results['padj'] < 0.05))|(deseq_results['padj'] > 0.05)].symbol.tolist()

notdiff = len(notdiff_genes)
#print(notdiff)


#print("Number of not differentially expressed genes", notdiff)
#print(upregulated_genes)



def map_color(a):
    log2FoldChange, symbol, padj, nlog10 = a

    if symbol in notdiff_genes:
        return 'Not_diff_expr'+' '+str(notdiff)
    elif symbol in upregulated_genes:
        return 'upregulated'+' '+str(upcount)
    elif symbol in downregulated_genes:
        return 'downregulated'+' '+str(downcount)




deseq_results['key'] = deseq_results[['log2FoldChange', 'symbol', 'padj', 'nlog10']].apply(map_color, axis = 1)


#print(deseq_results)
#make the plot

plt.figure(figsize = (6,6))

ax = sns.scatterplot(data = deseq_results, x = 'log2FoldChange', y = 'nlog10',
                    hue = 'key', hue_order = [ 'upregulated'+' '+str(upcount),'downregulated'+' '+str(downcount), 'Not_diff_expr'+' '+str(notdiff)],
                    palette = ['red','green', 'grey'])


ax.axhline(1, zorder = 0, c = 'k', lw = 1, ls = '--')
ax.axvline(-1, zorder = 0, c = 'k', lw = 1, ls = '--')
ax.axvline(1, zorder = 0, c = 'k', lw = 1, ls = '--')

#mark selected genes

mark_list=['Cdkn2A','Bst2','Grem1','Bst2','Foxj1','Gabrb3','Cers1']
texts = []
texts = [plt.text(row['log2FoldChange'], row['nlog10'], row['symbol'], fontsize = 12, weight = 'bold') for index, row in deseq_results.iterrows() if row['symbol'] in mark_list]

# Adjust the text
adjust_text(texts, arrowprops=dict(arrowstyle='-', color='k'))


#for index,row in deseq_results.iterrows():
#    if row['symbol'] in mark_list:
#        plt.annotate(row['symbol'],(row['log2FoldChange'],row['nlog10']), textcoords="offset points", xytext=(5, 5),
 #                    ha='left', bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.5))



plt.legend(loc = 1, bbox_to_anchor = (1.5,1), frameon = True, prop = {'weight':'bold'})

for axis in ['bottom', 'left']:
    ax.spines[axis].set_linewidth(2)

for axis in ['top', 'right']:
    ax.spines[axis].set_linewidth(2)

#ax.spines['top'].set_visible(True)
#ax.spines['right'].set_visible(False)

ax.tick_params(width = 2)

plt.xticks(size = 12, weight = 'bold')
plt.yticks(size = 12, weight = 'bold')

plt.xlabel("$log_{2}$ fold change", size = 15)
plt.ylabel("-$log_{10}$ FDR", size = 15)

plt.savefig('volcano_ctnnb_labeled_TPM.png', dpi = 300, bbox_inches = 'tight', facecolor = 'white')


plt.show()


