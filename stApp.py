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
    st.sidebar.title("WDMPLL?")
    st.sidebar.write("If you want to see your favourite molecular property included, drop a line at [@lewischewis](https://twitter.com/lewischewis) or ljmartin at hey dot com, or open a github issue")
    st.sidebar.write("""If you ask 'but why?' or 'but how?', see the readme at the [github page](https://github.com/ljmartin/what_do_mol_prop_look_like)""")
    st.sidebar.write('Click the ✖️ to close this bar and widen the view')


    #and some intro text in the main frame:
    st.title('What do molecular properties look like?')
    st.write("""The [Lipinski Ro5](https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five) helps people focus their drug discovery efforts on the molecules most likely to make good therapeutic drugs.""")
    st.write("""But, [increasingly](https://doi.org/10.1021/acs.jmedchem.8b00686), drug-like molecules break the Ro5, so it's helpful to push the boundaries of molecular properties when considering a molecule library. One way to get a feel for how far they can be pushed is to just stare at molecules in a certain property-space and decide if they look reasonable or not.""")
    st.write("""### Instructions""")
    st.write('There are sliders below that set the minimum or maximum Molecular Weight (MW) or calculated logP (cLogP). First, set a desired range. Then, click the "**Show Sample**" button. A small sample of 24 molecules satisfying the filters will be chosen and visualized. Just click it again to get a new batch.')

    st.write("""### Histograms """)
    st.write("If you set an unrealistic range, there won't be any molecules left. There are 500k molecules in the set, but the distribution isn't uniform. Here's a guide to help:")
    st.image('density.svg')
    

    

    st.write("""### Filters:""")
    ###now the app:
    
    #property sliders:    
    mw_min = st.slider('Molecular weight (MW) min:', 
                       min_value = mwmin,
                       max_value = mwmax,
                       #value = float(np.percentile(df['mw'], 5)),
                       value = (mwmax-mwmin)*0.05 + mwmin,
                       step=0.05
                       )
    mw_max = st.slider('Molecular weight (MW) max:', 
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

    st.write("""### Molecules:""")
    st.write('Number of molecules left: ', mask.sum())


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
            sample = df[mask].sample(min([N, mask.sum()]))
            ##2:
            mols = [Chem.MolFromSmiles(i) for i in sample['smiles']]
            
            ##3:
            st.image(Draw.MolsToGridImage(mols, molsPerRow=6, legends=list(sample['zinc_id'])))
        

if __name__=="__main__":
    main()
