# Comparisons calculations that we could make

import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from statannot import add_stat_annotation

df = pd.read_csv('data.csv')

def RegionStats(DataFrame):

    metrics = df['variable'].unique()
    metrics = np.delete(metrics,0)
    i = 0
    Met = []
    Reg = []
    Pv = []
    for metric in metrics:
        # visualizing the intensity distributions
        fig, ax = plt.subplots()
        ax = sns.boxplot(data=df.loc[(df['TIMEPT'] == 'BL') & (df['variable'] == metric)], x = 'Response', y = 'value')
        
        if 'VOL' in metric:
            ax.set_title('Average Tumour Volume')
            #ax.set_xlabel('Response')
            ax.set_xticklabels(labels = ['No Response','Response'])
            ax.set_ylabel('Volume ($mm^3$)')
            plt.savefig('Figure{0}.png'.format(i),bbox_inches = 'tight')
        elif 'mu' in metric:
            x = metric.split(' ')[0]
            ax.set_title('Average Hounsfield Units in {0}'.format(x))
            #ax.set_xlabel('Response')
            ax.set_xticklabels(labels = ['No Response','Response'])
            ax.set_ylabel('Hounsfield Units')
            plt.savefig('Figure{0}.png'.format(i),bbox_inches = 'tight')
        elif 'S' in metric:
            x = metric.split(' ')[0]
            ax.set_title('Average Entropy in {0}'.format(x))
            #ax.set_xlabel('Response')
            ax.set_xticklabels(labels = ['No Response','Response'])
            ax.set_ylabel('Entropy Units?')
            plt.savefig('Figure{0}.png'.format(i),bbox_inches = 'tight')
        # t-test between the rim and perivascular intensity calculations
        responders = df.loc[(df['TIMEPT'] == 'BL') & (df['variable'] == metric) & (df['Response'] == True)]['value'].values
        non_responders = df.loc[(df['TIMEPT'] == 'BL') & (df['variable'] == metric) & (df['Response'] == False)]['value'].values
        responders = responders[~np.isnan(responders)]
        non_responders = non_responders[~np.isnan(non_responders)]
        ttest = stats.ttest_ind(responders, non_responders)
        # keep here to inspect result
        Met.append(metric)
        Reg.append(metric.split(' ')[0])
        Pv.append(ttest.pvalue)
        i+=1
    results = pd.DataFrame(list(zip(Met,Reg,Pv)),
               columns =['Metrics', 'Region','P-value'])
    return results
results = RegionStats(df)
results.to_csv('Results.csv')