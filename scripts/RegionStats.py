# Comparisons calculations that we could make

import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from statannot import add_stat_annotation
import matplotlib.pyplot as plt

from statannotations.Annotator import Annotator




df = pd.read_csv('data.csv')

def RegionStats(DataFrame):

    metrics = df['variable'].unique()
    metrics = np.delete(metrics,0)
    i = 0
    Met = []
    Reg = []
    PosResp_count = []
    NonResp_count = []
    NegResp_count = []
    PosResp_Norm_pv = []
    NonResp_Norm_pv = []
    NegResp_Norm_pv = []
    Mw1_pv = []
    Mw2_pv = []
    Mw3_pv = []

    for metric in metrics:
        # visualizing the intensity distributions
        plt.style.use('dark_background')
        fig, ax = plt.subplots()
        ax = sns.boxplot(data=df.loc[(df['TIMEPT'] == 'BL') & (df['variable'] == metric) & (df['Response'] != 9)], x = 'Response', y = 'value', order = ['r','n','g'])
        x = "Response"
        y = "value"
        order = ['r','n','g']
        pairs=[('r','n'),('r','g'),('n','g')]
        test = df.loc[(df['TIMEPT'] == 'BL') & (df['variable'] == metric) & (df['Response'] != 9)]
        annotator = Annotator(ax, pairs, data=test, x=x, y=y, order=order)
        annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
        annotator.apply_and_annotate()

        if 'VOL' in metric:
            ax.set_title('Baseline tumour volume', fontsize = fs+2)
            ax.set_xlabel('', fontsize = fs)
            ax.set_xticklabels(labels = ['Volume Reduced', 'No Response','Volume Increased'], fontsize = fs)
            ax.set_ylabel('Volume ($cm^3$)', fontsize = fs)
            plt.savefig('Figure{0}.png'.format(i),bbox_inches = 'tight')
        elif 'mu' in metric:
            x = metric.split(' ')[0]
            ax.set_title('CT intensity in {0}'.format(x), fontsize = fs+2)
            ax.set_xlabel('', fontsize = fs)
            ax.set_xticklabels(labels = ['Volume Reduced', 'No Response','Volume Increased'], fontsize = fs)
            ax.set_ylabel('CT intensity (HU)', fontsize = fs)
            ax.set_ylim(-900,400)
            plt.savefig('Figure{0}.png'.format(i),bbox_inches = 'tight')
        elif 'S' in metric:
            x = metric.split(' ')[0]
            ax.set_title('Entropy in {0}'.format(x), fontsize = fs+2)
            ax.set_xlabel('', fontsize = fs)
            ax.set_xticklabels(labels = ['Volume Reduced', 'No Response','Volume Increased'], fontsize = fs)
            ax.set_ylabel('Entropy (1)', fontsize = fs)
            ax.set_ylim(1,3.0)
            plt.savefig('Figure{0}.png'.format(i),bbox_inches = 'tight')



        # t-test between the rim and perivascular intensity calculations
        positive_responders = df.loc[(df['TIMEPT'] == 'BL') & (df['variable'] == metric) & (df['Response'] == 'r')]['value'].values
        non_responders = df.loc[(df['TIMEPT'] == 'BL') & (df['variable'] == metric) & (df['Response'] == 'n')]['value'].values
        negative_responders = df.loc[(df['TIMEPT'] == 'BL') & (df['variable'] == metric) & (df['Response'] == 'g')]['value'].values

        positive_responders = positive_responders[~np.isnan(positive_responders)]
        non_responders = non_responders[~np.isnan(non_responders)]
        negative_responders = negative_responders[~np.isnan(negative_responders)]


        #ttest = stats.ttest_ind(positive_responders, negative_responders, equal_var=False)
        manwit1 = stats.mannwhitneyu(positive_responders,non_responders)
        manwit2 = stats.mannwhitneyu(positive_responders,negative_responders)
        manwit3 = stats.mannwhitneyu(non_responders,negative_responders)
        #kruk = stats.kruskal(positive_responders, non_responders, negative_responders)

        # keep here to inspect result
        Met.append(metric)
        Reg.append(metric.split(' ')[0])

        PosResp_count.append(len(positive_responders))
        NonResp_count.append(len(non_responders))
        NegResp_count.append(len(negative_responders))

        PosResp_Norm_pv.append(stats.shapiro(positive_responders).pvalue)
        NonResp_Norm_pv.append(stats.shapiro(non_responders).pvalue)
        NegResp_Norm_pv.append(stats.shapiro(negative_responders).pvalue)


        Mw1_pv.append(manwit1.pvalue)
        Mw2_pv.append(manwit2.pvalue)
        Mw3_pv.append(manwit3.pvalue)
        i+=1
    results = pd.DataFrame(list(zip(Met,Reg,PosResp_count,NonResp_count,NegResp_count,PosResp_Norm_pv,NonResp_Norm_pv,NegResp_Norm_pv,Mw1_pv,Mw2_pv,Mw3_pv)),
               columns =['Metrics', 'Region','Positive Responder Sample Size','Non-Responder Sample Size','Negative Responder Sample Size',
               'Responder Normality P-value', 'Non-Responder Normality P-value','Negative Responder Normality P-value', 'Reduction/No Response P-value', 'Reduction/Increase P-value', 'No Response/Increase P-value'])
    return results


fs = 12

results = RegionStats(df)
results.to_csv('Results.csv')

        
