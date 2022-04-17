from cmath import isnan
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def Remove_Outliers(data, threshold=3.5):
    """ (DataFrame, float) --> DataFrame

    Median Absolute Deviation (MAD) based outlier detection based on the work by Iglewicz and Hoaglin
    B. Iglewicz and D. C. Hoaglin, Volume 16: How to Detect and Handle Outliers", The ASQC Basic References in Quality Control: Statistical Techniques, vol. 36, no. 3. 1993.
    """
    # Median of the data
    median = np.median(data)

    # Median absolute deviation
    med_abs_deviation = np.median(abs(data - median))

    # Modified z-score calculation (as recommended by Iglewicz and Hoaglin)
    modified_z_score = 0.6745 * (data - median)/ med_abs_deviation
    data[modified_z_score > threshold] = np.nan

    return data 

def RemoveMissingPairs(DataFrame):
    ''' (Dataframe) --> None

    Removes data pairs missing baseline or post treatment data points
    '''

    # loop through each unique patient
    for pt in DataFrame['USUBJID'].unique():
        criterion_1 = (DataFrame['USUBJID'] == pt)

        # loop through each lesion for each patient
        for lesion_index in DataFrame.loc[criterion_1]['LSNIND'].unique():
            criterion_2 = (DataFrame['LSNIND'] == lesion_index)
            criterion_3a = (DataFrame['TIMEPT'] == 'BL')
            criterion_3b = (DataFrame['TIMEPT'] == 'C2')

            # find the pre/post volumes and indices
            vol_pre = DataFrame.loc[criterion_1 & criterion_2 & criterion_3a]['LSNVOL'].values
            vol_post = DataFrame.loc[criterion_1 & criterion_2 & criterion_3b]['LSNVOL'].values
            vol_pre_index = DataFrame.index[criterion_1 & criterion_2 & criterion_3a]
            vol_post_index = DataFrame.index[criterion_1 & criterion_2 & criterion_3b]
            
            # flag missing pairs
            if np.size(vol_pre) == 0 or np.size(vol_post) == 0:
                DataFrame.loc[vol_pre_index, 'Response'] = np.nan
                DataFrame.loc[vol_post_index, 'Response'] = np.nan
            else:
                DataFrame.loc[vol_pre_index, 'Response'] = 9
                DataFrame.loc[vol_post_index, 'Response'] = 9

    # remove missing pairs            
    DataFrame = DataFrame.dropna()
    DataFrame = DataFrame.reset_index(drop=True)



def ResponseAssignment(DataFrame):
    ''' (dataframe) --> dataframe

    Add response flag based on ratio of baseline and post treatment tumour volumes
    '''

    # loop through each unique patient
    for pt in DataFrame['USUBJID'].unique():

        criterion_1 = (DataFrame['USUBJID'] == pt)

        # loop through each lesion for each patient
        for lesion_index in DataFrame.loc[criterion_1]['LSNIND'].unique():
            criterion_2 = (DataFrame['LSNIND'] == lesion_index)
            criterion_3a = (DataFrame['TIMEPT'] == 'BL')
            criterion_3b = (DataFrame['TIMEPT'] == 'C2')

            # find the pre/post volumes and indices
            vol_pre = DataFrame.loc[criterion_1 & criterion_2 & criterion_3a]['LSNVOL'].values
            vol_post = DataFrame.loc[criterion_1 & criterion_2 & criterion_3b]['LSNVOL'].values
            vol_pre_index = DataFrame.index[criterion_1 & criterion_2 & criterion_3a]
            vol_post_index = DataFrame.index[criterion_1 & criterion_2 & criterion_3b]

            # calculate volume ratio
            ratio = vol_post/vol_pre
                
            # create response flag in the dataframe
            if ratio <= 0.85:
                DataFrame.loc[vol_pre_index, 'Response'] = 'Reduction'
                DataFrame.loc[vol_post_index, 'Response'] = 'Reduction'
            elif 0.85 < ratio < 1.15:
                DataFrame.loc[vol_pre_index, 'Response'] = 'No Change'
                DataFrame.loc[vol_post_index, 'Response'] = 'No Change'
            else:
                DataFrame.loc[vol_pre_index, 'Response'] = 'Growth'
                DataFrame.loc[vol_post_index, 'Response'] = 'Growth'


# read in dataframe with segmentation output
df = pd.read_csv('results_5Apr2022.csv')

# filter to treatment arm only
df = df.loc[(df['ARM'] ==True)]

# filter any missing data
df = df.dropna()
df = df.reset_index(drop=True)

# change volume to mm^3 from cm^3
df['LSNVOL'] = df['LSNVOL'].div(1000)

# flag any lesions with 0 volume as NaN and remove
df['LSNVOL'] = df['LSNVOL'].replace(0,np.nan)
df = df.dropna()
df = df.reset_index(drop=True)

# separate baseline and post treatment groups to use median filter
df2 = df.loc[df['TIMEPT'] == 'BL']
df3 = df.loc[df['TIMEPT'] == 'C2']

# use median filter based on volume
df2['LSNVOL'] = Remove_Outliers(df2['LSNVOL'])
df3['LSNVOL'] = Remove_Outliers(df3['LSNVOL'])
df2 = df2.dropna()
df3 = df3.dropna()
df2 = df2.reset_index(drop=True)
df3 = df3.reset_index(drop=True)

# recombine into 1 dataframe
df = pd.concat([df2,df3])

# remove any pairs with missing data
RemoveMissingPairs(df)

# assign response based on volume ratio
ResponseAssignment(df)

# reshape dataframe to be more friendly for stats functions
headers = ['USUBJID','TIMEPT','ARM','LSNIND','Response']
df = pd.melt(df,id_vars=headers)
df.to_csv('data.csv')
            
