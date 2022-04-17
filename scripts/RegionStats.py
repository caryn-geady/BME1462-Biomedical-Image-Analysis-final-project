# Comparisons calculations that we could make

import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from statannotations.Annotator import Annotator

def CreateBP(DataFrame):
    ''' (Dataframe) --> None

    Creates box plots of volume, intensity, and entropy from the 3 response groups defined by the ResponseAssignment.py script
    '''

    # filter DF to only contain baseline image data
    DataFrame = DataFrame.loc[(DataFrame['TIMEPT'] == 'BL')]
    
    # Divide DF into volume (d1), intensity (d2), and entropy (d3)
    d1= DataFrame[DataFrame['variable'].str.contains('VOL')]
    d2 = DataFrame[DataFrame['variable'].str.contains('mu')]
    d3 = DataFrame[DataFrame['variable'].str.contains('S')]

    # Assign regions for intensity and entropy data
    # Region 1 = whole lesion, 2 = core, 3 = internal rim, 4 = external rim
    d2a = d2[d2['variable'].str.contains('lesion')].assign(Region = 1)
    d2b = d2[d2['variable'].str.contains('core')].assign(Region = 2)
    d2c = d2[d2['variable'].str.contains('interior')].assign(Region = 3)
    d2d = d2[d2['variable'].str.contains('exterior')].assign(Region = 4)

    d3a = d3[d3['variable'].str.contains('lesion')].assign(Region = 1)
    d3b = d3[d3['variable'].str.contains('core')].assign(Region = 2)
    d3c = d3[d3['variable'].str.contains('interior')].assign(Region = 3)
    d3d = d3[d3['variable'].str.contains('exterior')].assign(Region = 4)

    # create intensity and entropy dataframes
    mu = pd.concat([d2a,d2b,d2c,d2d])
    S = pd.concat([d3a,d3b,d3c,d3d])

    # Figure 1 - Baseline lesion volume comparison
    fig, ax = plt.subplots()

    ax = sns.boxplot(data=d1, 
                    x = 'Response', y = 'value', 
                    order = ['Growth','Reduction','No Change'], 
                    ax = ax)
    x = "Response"
    y = "value"
    pairs=[('Reduction','No Change'),('Reduction','Growth'),('No Change','Growth')]
    annotator = Annotator(ax, pairs, data=d1, x=x, y=y)
    annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
    annotator.apply_and_annotate()
    ax.set_title('Baseline Tumor Volume', fontsize = fs+2)
    ax.set_xlabel('Region: Whole Lesion', fontsize = fs)
    ax.set_xticklabels(labels = '', fontsize = fs)
    ax.set_ylabel('Volume ($cm^3$)', fontsize = fs)
    plt.savefig('Figure1.png',bbox_inches = 'tight')

    # Figure 2 - Mean CT intensity comparison
    fig, ax = plt.subplots()
    pairs=[
        [(1,'Reduction'),(1,'No Change')],
        [(1,'Reduction'),(1,'Growth')],
        [(1,'No Change'),(1,'Growth')],

        [(2,'Reduction'),(2,'No Change')],
        [(2,'Reduction'),(2,'Growth')],
        [(2,'No Change'),(2,'Growth')],

        [(3,'Reduction'),(3,'No Change')],
        [(3,'Reduction'),(3,'Growth')],
        [(3,'No Change'),(3,'Growth')],

        [(4,'Reduction'),(4,'No Change')],
        [(4,'Reduction'),(4,'Growth')],
        [(4,'No Change'),(4,'Growth')],
        ]
    hue_plot_params = {
        'data': mu,
        'x': 'Region',
        'y': 'value',
        'hue': 'Response'
    }
    ax = sns.boxplot(ax = ax,**hue_plot_params)
    annotator = Annotator(ax,pairs,**hue_plot_params)
    annotator.configure(test = "Mann-Whitney", text_format='star', loc='inside').apply_and_annotate()
    ax.set_title('Mean CT Intensity', fontsize = fs+2)
    ax.set_xlabel('Region', fontsize = fs)
    ax.set_xticklabels(labels = ['Whole Lesion', 'Core','Interior Rim', 'Exterior Rim'], fontsize = fs)
    ax.set_ylabel('CT intensity (HU)', fontsize = fs)
    ax.set_ylim(-900,500)
    ax.get_legend().remove()
    plt.savefig('Figure2.png',bbox_inches = 'tight')

    # Figure 3 - Entropy comparison
    fig, ax = plt.subplots()
    pairs=[
        [(1,'Reduction'),(1,'No Change')],
        [(1,'Reduction'),(1,'Growth')],
        [(1,'No Change'),(1,'Growth')],

        [(2,'Reduction'),(2,'No Change')],
        [(2,'Reduction'),(2,'Growth')],
        [(2,'No Change'),(2,'Growth')],

        [(3,'Reduction'),(3,'No Change')],
        [(3,'Reduction'),(3,'Growth')],
        [(3,'No Change'),(3,'Growth')],

        [(4,'Reduction'),(4,'No Change')],
        [(4,'Reduction'),(4,'Growth')],
        [(4,'No Change'),(4,'Growth')],
        ]
    hue_plot_params = {
        'data': S,
        'x': 'Region',
        'y': 'value',
        'hue': 'Response'
    }
    ax = sns.boxplot(ax = ax,**hue_plot_params)
    annotator = Annotator(ax,pairs,**hue_plot_params)
    annotator.configure(test = "Mann-Whitney", text_format='star', loc='inside').apply_and_annotate()
    ax.set_title('Entropy', fontsize = fs+2)
    ax.set_xlabel('Region', fontsize = fs)
    ax.set_xticklabels(labels = ['Whole Lesion', 'Core','Interior Rim', 'Exterior Rim'], fontsize = fs)
    ax.set_ylabel('Entropy (1)', fontsize = fs)
    ax.set_ylim(1,3)
    ax.get_legend().remove()
    plt.savefig('Figure3.png',bbox_inches = 'tight')


def StatsTable(DataFrame):
    ''' (Dataframe) -> Dataframe

    Returns a dataframe containing the statistical test results between the 3 response groups defined by the ResponseAssignment.py script
    '''
    # filter DF to only contain baseline image data
    DataFrame = DataFrame.loc[(DataFrame['TIMEPT'] == 'BL')]

    # Divide DF into volume (d1), intensity (d2), and entropy (d3)
    d1= DataFrame[DataFrame['variable'].str.contains('VOL')].assign(Region = 1)
    d2 = DataFrame[DataFrame['variable'].str.contains('mu')]
    d3 = DataFrame[DataFrame['variable'].str.contains('S')]

    # Assign regions for intensity and entropy data
    # Region 1 = whole lesion, 2 = core, 3 = internal rim, 4 = external rim
    d2a = d2[d2['variable'].str.contains('lesion')].assign(Region = 1)
    d2b = d2[d2['variable'].str.contains('core')].assign(Region = 2)
    d2c = d2[d2['variable'].str.contains('interior')].assign(Region = 3)
    d2d = d2[d2['variable'].str.contains('exterior')].assign(Region = 4)

    d3a = d3[d3['variable'].str.contains('lesion')].assign(Region = 1)
    d3b = d3[d3['variable'].str.contains('core')].assign(Region = 2)
    d3c = d3[d3['variable'].str.contains('interior')].assign(Region = 3)
    d3d = d3[d3['variable'].str.contains('exterior')].assign(Region = 4)

    # create intensity and entropy dataframes
    mu = pd.concat([d2a,d2b,d2c,d2d])
    S = pd.concat([d3a,d3b,d3c,d3d])

    # define lists to store parameter values
    Metric = []
    Region = []
    Sample_Count_Growth = []
    Sample_Count_Reduction = []
    Sample_Count_NoChange = []
    Normality_Growth = []
    Normality_Reduction =[]
    Normality_NoChange = []
    MannWhitney_Growth_Reduction = []
    MannWhitney_Growth_NoChange = []
    MannWhitney_Reduction_NoChange = []


    # Loop through different regions for each metric

    # volume dataframe
    for metric in d1['variable'].unique():
        
        # keep track of the metric and region
        Metric.append('Lesion Volume')
        Region.append('lesion')

        # divide into different response groups
        Response_Growth = d1.loc[(df['variable'] == metric) & (df['Response'] == 'Growth')]['value'].values
        Response_Reduction = d1.loc[(df['variable'] == metric) & (df['Response'] == 'Reduction')]['value'].values
        Response_NoChange = d1.loc[(df['variable'] == metric) & (df['Response'] == 'No Change')]['value'].values

        # keep track of sample size for each group
        Sample_Count_Growth.append(len(Response_Growth))
        Sample_Count_Reduction.append(len(Response_Reduction))
        Sample_Count_NoChange.append(len(Response_NoChange))

        # Normality test for the different groups
        Normality_Growth.append(stats.shapiro(Response_Growth).pvalue)
        Normality_Reduction.append(stats.shapiro(Response_Reduction).pvalue)
        Normality_NoChange.append(stats.shapiro(Response_NoChange).pvalue)

        # mann-whitney u test for the different groups
        MannWhitney_Growth_Reduction.append(stats.mannwhitneyu(Response_Growth,Response_Reduction).pvalue)
        MannWhitney_Growth_NoChange.append(stats.mannwhitneyu(Response_Growth,Response_NoChange).pvalue)
        MannWhitney_Reduction_NoChange.append(stats.mannwhitneyu(Response_Reduction,Response_NoChange).pvalue)
    
    # intensity dataframe
    for metric in mu['variable'].unique():
        
        # keep track of the metric and region
        Metric.append('Intensity')
        Region.append(metric.split(' ')[0])

        # divide into different response groups
        Response_Growth = mu.loc[(df['variable'] == metric) & (df['Response'] == 'Growth')]['value'].values
        Response_Reduction = mu.loc[(df['variable'] == metric) & (df['Response'] == 'Reduction')]['value'].values
        Response_NoChange = mu.loc[(df['variable'] == metric) & (df['Response'] == 'No Change')]['value'].values

        # keep track of sample size for each group
        Sample_Count_Growth.append(len(Response_Growth))
        Sample_Count_Reduction.append(len(Response_Reduction))
        Sample_Count_NoChange.append(len(Response_NoChange))

        # Normality test for the different groups
        Normality_Growth.append(stats.shapiro(Response_Growth).pvalue)
        Normality_Reduction.append(stats.shapiro(Response_Reduction).pvalue)
        Normality_NoChange.append(stats.shapiro(Response_NoChange).pvalue)

        # mann-whitney u test for the different groups
        MannWhitney_Growth_Reduction.append(stats.mannwhitneyu(Response_Growth,Response_Reduction).pvalue)
        MannWhitney_Growth_NoChange.append(stats.mannwhitneyu(Response_Growth,Response_NoChange).pvalue)
        MannWhitney_Reduction_NoChange.append(stats.mannwhitneyu(Response_Reduction,Response_NoChange).pvalue)
    
    # entropy dataframe
    for metric in S['variable'].unique():
        
        # keep track of the metric and region
        Metric.append('Entropy')
        Region.append(metric.split(' ')[0])

        # divide into different response groups
        Response_Growth = S.loc[(df['variable'] == metric) & (df['Response'] == 'Growth')]['value'].values
        Response_Reduction = S.loc[(df['variable'] == metric) & (df['Response'] == 'Reduction')]['value'].values
        Response_NoChange = S.loc[(df['variable'] == metric) & (df['Response'] == 'No Change')]['value'].values

        # keep track of sample size for each group
        Sample_Count_Growth.append(len(Response_Growth))
        Sample_Count_Reduction.append(len(Response_Reduction))
        Sample_Count_NoChange.append(len(Response_NoChange))

        # Normality test for the different groups
        Normality_Growth.append(stats.shapiro(Response_Growth).pvalue)
        Normality_Reduction.append(stats.shapiro(Response_Reduction).pvalue)
        Normality_NoChange.append(stats.shapiro(Response_NoChange).pvalue)

        # mann-whitney u test for the different groups
        MannWhitney_Growth_Reduction.append(stats.mannwhitneyu(Response_Growth,Response_Reduction).pvalue)
        MannWhitney_Growth_NoChange.append(stats.mannwhitneyu(Response_Growth,Response_NoChange).pvalue)
        MannWhitney_Reduction_NoChange.append(stats.mannwhitneyu(Response_Reduction,Response_NoChange).pvalue)


    results = pd.DataFrame(list(zip(Metric,
                                    Region,
                                    Sample_Count_Growth, Sample_Count_Reduction, Sample_Count_NoChange, 
                                    Normality_Growth, Normality_Reduction, Normality_NoChange,
                                    MannWhitney_Growth_Reduction, MannWhitney_Growth_NoChange, MannWhitney_Reduction_NoChange)),
               columns =['Metrics',
                         'Region',
                         'Sample Size: Growth', 'Sample Size: Reduction','Sample Size: No Change',
                         'Shapiro Wilk P-Value: Growth','Shapiro Wilk P-Value: Reduction','Shapiro Wilk P-Value: No Change',
                         'Mann-Whitney P-Value: G-R', 'Mann-Whitney P-Value: G-NC', 'Mann-Whitney P-Value: R-NC'
               ])
    return results

# Read in output from Response Assignment script
df = pd.read_csv('data.csv')

# basic font size for the plots
fs = 12

# Call function to create box plots
CreateBP(df)

table = StatsTable(df)
table.to_csv('ResultsTable.csv')

        
