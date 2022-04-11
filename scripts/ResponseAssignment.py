import pandas as pd


# define dataframe with patient/lesion data
df = pd.read_csv('results_5Apr2022.csv')

def ResponseAssignment(DataFrame):
    ''' (dataframe) --> dataframe
    Add response flag based on pre and post tumour volumes
    '''
    # loop through each unique patient
    for pt in df['USUBJID'].unique():
        criterion_1 = (df['USUBJID'] == pt)
        # loop through each lesion for each patient
        for lesion_index in df.loc[criterion_1]['LSNIND'].unique():
            criterion_2 = (df['LSNIND'] == lesion_index)
            criterion_3a = (df['TIMEPT'] == 'BL')
            criterion_3b = (df['TIMEPT'] == 'C2')

            # find the pre/post volumes and indices
            vol_pre = df.loc[criterion_1 & criterion_2 & criterion_3a]['LSNVOL'].values
            vol_post = df.loc[criterion_1 & criterion_2 & criterion_3b]['LSNVOL'].values
            vol_pre_index = df.index[criterion_1 & criterion_2 & criterion_3a]
            vol_post_index = df.index[criterion_1 & criterion_2 & criterion_3b]
            
            # ratio of post vs pre volume
            ratio = vol_post/vol_pre
            
            # create response flag in the dataframe
            if ratio <= 0.8:
                df.loc[vol_pre_index, 'Response'] = True
                df.loc[vol_post_index, 'Response'] = True
            else:
                df.loc[vol_pre_index, 'Response'] = False
                df.loc[vol_post_index, 'Response'] = False

ResponseAssignment(df)
headers = ['USUBJID','TIMEPT','ARM','LSNIND','Response']
df = pd.melt(df,id_vars=headers)
df.to_csv('data.csv')
            
