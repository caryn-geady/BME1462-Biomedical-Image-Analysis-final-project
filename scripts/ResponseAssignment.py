from cmath import isnan
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def remove_outliers2(data, threshold=3.5):
        """
        Median Absolute Deviation (MAD) based outlier detection
        Removes outliers and if selected fills with polynomial  interpolation
        fill: Boolean
        https://www.programcreek.com/python/?CodeExample=remove+outliers
        https://www.itl.nist.gov/div898/handbook/eda/section3/eda356.htm#MAD
        https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
        """
        median = np.median(data)

        med_abs_deviation = np.median(abs(data - median))
        # scale constant 0.6745
        modified_z_score = 0.6745 * (data - median)/ med_abs_deviation
        data[modified_z_score > threshold] = np.nan

        return data 
def RemoveMissingPairs(DataFrame):
    ''' (dataframe) --> dataframe
    Add response flag based on pre and post tumour volumes
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
            
            if np.size(vol_pre) == 0 or np.size(vol_post) == 0:
                DataFrame.loc[vol_pre_index, 'Response'] = np.nan
                DataFrame.loc[vol_post_index, 'Response'] = np.nan
            else:
                DataFrame.loc[vol_pre_index, 'Response'] = 9
                DataFrame.loc[vol_post_index, 'Response'] = 9

def ResponseAssignment(DataFrame):
    ''' (dataframe) --> dataframe
    Add response flag based on pre and post tumour volumes
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
            
            if np.size(vol_pre) == 0 or np.size(vol_post) == 0:
                DataFrame.loc[vol_pre_index, 'Response'] = np.nan
                DataFrame.loc[vol_post_index, 'Response'] = np.nan
                continue
            # ratio of post vs pre volume
            ratio = vol_post/vol_pre
                
            # create response flag in the dataframe
            if ratio <=0.1:
                DataFrame.loc[vol_pre_index, 'Response'] = 'r'
                DataFrame.loc[vol_post_index, 'Response'] = 'r'
            elif ratio <= 0.85:
                DataFrame.loc[vol_pre_index, 'Response'] = 'r'
                DataFrame.loc[vol_post_index, 'Response'] = 'r'
            elif 0.85 < ratio < 1.15:
                DataFrame.loc[vol_pre_index, 'Response'] = 'n'
                DataFrame.loc[vol_post_index, 'Response'] = 'n'
            elif 1.15 <= ratio <= 1.9:
                DataFrame.loc[vol_pre_index, 'Response'] = 'g'
                DataFrame.loc[vol_post_index, 'Response'] = 'g'
            else:
                DataFrame.loc[vol_pre_index, 'Response'] = 'g'
                DataFrame.loc[vol_post_index, 'Response'] = 'g'

# define dataframe with patient/lesion data
df = pd.read_csv('results_5Apr2022.csv')

df = df.loc[(df['ARM'] ==True)]
print(df.shape[0])
df['LSNVOL'] = df['LSNVOL'].div(1000)
df['LSNVOL'] = df['LSNVOL'].replace(0,np.nan)


# drop nan values
df = df.dropna()
df = df.reset_index(drop=True)

RemoveMissingPairs(df)

df = df.dropna()
df = df.reset_index(drop=True)

df2 = df.loc[df['TIMEPT'] == 'BL']
df3 = df.loc[df['TIMEPT'] == 'C2']
df2['LSNVOL'] = remove_outliers2(df2['LSNVOL'])
df2 = df2.dropna()
df2 = df2.reset_index(drop=True)

df3['LSNVOL'] = remove_outliers2(df3['LSNVOL'])
df3 = df3.dropna()
df3 = df3.reset_index(drop=True)

df = pd.concat([df2,df3])
ResponseAssignment(df)

df = df.dropna()
df = df.reset_index(drop=True)


#df = df.loc[df['TIMEPT'] == 'BL']

fs = 14
plt.style.use('dark_background')
#df['LSNVOL'].hist(bins=range(0,500,50), color = 'blue')
df['LSNVOL'].hist(bins=10, color = 'blue')
plt.xlabel('Volume ($cm^3$)', fontsize = fs)
plt.ylabel('Frequency', fontsize = fs)
#plt.xlim(0,400)
plt.xlim(0,10)
plt.xticks(fontsize = fs)
#plt.ylim(0,140)
plt.ylim(0,20)
plt.yticks(fontsize = fs)

#plt.title('Lesion volume histogram: prior to filtering', fontsize = fs+2)
plt.title('Lesion volume histogram: post filtering', fontsize = fs+2)
plt.show()

print(df.shape[0])
# reshape dataframe to be more friendly for stats functions
headers = ['USUBJID','TIMEPT','ARM','LSNIND','Response']
df = pd.melt(df,id_vars=headers)
df.to_csv('data.csv')
            
