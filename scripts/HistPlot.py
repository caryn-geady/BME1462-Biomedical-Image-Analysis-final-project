import matplotlib.pyplot as plt
import numpy as np
from cmath import isnan
import pandas as pd


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

        data = data.dropna()
        return data 

# define dataframe with patient/lesion data
df = pd.read_csv('results_5Apr2022.csv')

# removing rows with missing values
# needed this and the size check in the function to get rid of all nan/zeros
df = df.dropna()
criterion_4 = (df['ARM'] ==True)
df['LSNVOL'] = remove_outliers2(df.loc[criterion_4]['LSNVOL'])


df = df.loc[criterion_4]['LSNVOL']
#df.hist(bins=range(0,500000,50000), color = 'green')
df.hist(bins=10, color = 'blue')
plt.xlabel('Volume ($mm^3$)')
plt.ylabel('Frequency')
plt.xlim(0,400000)
plt.ylim(0,120)
plt.show()


            
