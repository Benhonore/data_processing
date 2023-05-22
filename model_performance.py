import pandas as pd
import numpy as np
from sklearn.metrics import mean_absolute_error as mae
from sklearn.metrics import mean_squared_error as mse
from math import sqrt

### Passing this function a dataframe that contans true DFT shifts and IMPRESSION predicted shifts will print put a model performance summary for carbon and proton.

ref={'0':'10','1':'9','2':'8','3':'7','4':'6', '5':'5', '6':'4', '7': '3', '8': '2', '9': '1', '10': '0'}

def performance_summary(df):
    carbon_df=df[df['typestr']=='C']
    print('CARBON')
    print(f"MAE:{mae(carbon_df['shift'], carbon_df['predicted_shift'])}")
    print(f"RMSE:{sqrt(mse(carbon_df['shift'], carbon_df['predicted_shift']))}")
    print(f"MaxE:{np.max(abs(carbon_df['shift']-carbon_df['predicted_shift']))}")
    print('')
    proton_df=df[df['typestr']=='H']
    print('PROTON')
    print(f"MAE:{mae(proton_df['shift'], proton_df['predicted_shift'])}")
    print(f"RMSE:{sqrt(mse(proton_df['shift'], proton_df['predicted_shift']))}")
    print(f"MaxE:{np.max(abs(proton_df['shift']-proton_df['predicted_shift']))}")


### This function will rank the closeness in shift prediction for each conformer to each other conformer in a dataframe of IMPRESSION predicted shifts for multiple
### conformers for the atom type and molecule parsed.
    
def rank_conf_prediction(name, nucleus, shift_df):

    df=pd.DataFrame()
    for i in range(50):
        try:
            df=df.append(shift_df[shift_df['molecule_name']==f'{name}_conf_{i}'], ignore_index=True)
        except:
            print(f'{i}th conformer missing')
            continue

    num_confs=len(df["molecule_name"].unique())
    #print(f'{num_confs} conformers')
    if len(df['molecule_name'].unique()) < 10:
        count=len(df['molecule_name'].unique())
    else:
        count=10

    dic={}

    for i in df['molecule_name'].unique():
        pred_dt6_shifts = df[df['molecule_name']==i][df['typestr']==f'{nucleus}']['predicted_shift']


        errors=[]
        for j in df['molecule_name'].unique():
            true_dt6_shifts = df[df['molecule_name']==j][df['typestr']==f'{nucleus}']['shift']

            errors.append(mae(pred_dt6_shifts, true_dt6_shifts))
        dic[i]=errors

    c_df=pd.DataFrame(dic)
    c_df.insert(loc=0, column ='Name', value=c_df.columns)

    comp={}
    for i in c_df.columns[1:]:
        k=[]
        lis = list(c_df[i])
        for num in range(count):
            val = np.min(lis)
            lis.remove(val)
            index=c_df.index[c_df[i]==val][0]
            k.append(str(c_df.iloc[index]['Name']))

        comp[i] = k

    c=0
    for i in comp:
        if i == comp[i][0]:
            c+=1


    dd=pd.DataFrame(comp)

    f={}
    for i in dd.columns:
        rank=[]
        for j in dd.columns:
            if j in list(dd[i]):
                rank.append(ref[str(dd.index[dd[i]==j][0])])
            else:
                rank.append(0)
        f[i]=rank

    c_rank_df = pd.DataFrame(f)
    c_rank_df.insert(loc=0, column ='Name', value=c_rank_df.columns)

    return c, num_confs, comp, c_rank_df, c_df
