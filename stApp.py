import streamlit as st

import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import Draw

@st.cache
def load_df():
    return pd.read_csv('mydf.csv')


st.set_page_config(
    layout='wide',
    )


def main():
    #first things first, load the dataframe of molecules and their properties:
    #df = pd.read_csv('sample.smifi')
    df = load_df()
    #handy values
    mwmin = float(df['mw'].min()-1)
    mwmax = float(df['mw'].max()+1)
    clogpmin = float(df['clogp'].min()-1)
    clogpmax = float(df['clogp'].max()-1)


    #print out some explanation stuff in the sidebar:
    st.sidebar.title("Start here:")


    #and some intro text in the main frame:
    st.title('Welcome')
    st.write('Instructions on how to use :)')


    st.image('density.svg')
    

    


    ###now the app:
    
    #property sliders:    
    mw_min = st.slider('MW_min', 
                       min_value = mwmin,
                       max_value = mwmax,
                       #value = float(np.percentile(df['mw'], 5)),
                       value = (mwmax-mwmin)*0.05 + mwmin,
                       step=0.05
                       )
    mw_max = st.slider('MW_max', 
                       min_value = mwmin,
                       max_value = mwmax,
                       value = (mwmax-mwmin)*0.95 + mwmin,
                       )

    
    clogp_min = st.slider('cLogP min',
                          min_value = clogpmin,
                          max_value = clogpmax,
                          value = (clogpmax-clogpmin)*0.05 + clogpmin,
                         )
    clogp_max = st.slider('cLogP max',
                          min_value = clogpmin,
                          max_value = clogpmax,
                          value = (clogpmax-clogpmin)*0.95 + clogpmin
                         )

    
    
    mask = (df['mw'] <= mw_max) & (df['mw'] >= mw_min) \
            & (df['clogp'] <= clogp_max) & (df['clogp'] >= clogp_min)
    st.write('Number of molecules: ', mask.sum())


    #this is the main event. Based on the filters/sliders above:
    #1. select a random sample of N ligands that meet the selected filter.
    #2. turn them into molecules,
    #3. and draw!
    N = 24
    if st.button('Show sample'):

        ##1:
        mask = (df['mw'] <= mw_max) & (df['mw'] >= mw_min) \
            & (df['clogp'] <= clogp_max) & (df['clogp'] >= clogp_min)

        ##1.5: quick error check:
        flag = mask.sum()>0
        if not flag:
            st.write('Set the property filters again - there are no molecules that fit those parameters')

        ##good to go.
        else:
            sample = df[mask].sample(N)
            ##2:
            mols = [Chem.MolFromSmiles(i) for i in sample['smiles']]
            
            ##3:
            st.image(Draw.MolsToGridImage(mols, molsPerRow=6, legends=list(sample['zinc_id'])))
        

if __name__=="__main__":
    main()
