# What do molecular properties look like?

view the app:
[link](https://share.streamlit.io/ljmartin/what_do_mol_prop_look_like/main/stApp.py)

### what
this a streamlit app centred around 500,000 molecules from the ZINC database. You can drag the sliders to filter the molecules by their physical chemistry properties (so far, only molecular weight and calculated logP). You then press a button to select a small random sample and visualize them for inspection.


### why
[ultra-large molecular docking](https://doi.org/10.1038/s41586-019-0917-9) is in vogue now. A typical virtual library can easily reach hundreds of millions of molecules, requiring high-performance computing clusters. One way to cut down on computation time is to pre-filter the library by properties you think are desirable. The [ZINC tranche browser](http://zinc15.docking.org/tranches/home/) is one way to do this. But how do you decide where to set the cut-off? The five rules from [Lipinsk, Lombardo, Dominy, and Feeney](https://doi.org/10.1016/S0169-409X(96)00423-1) were the original way, but more and more people are seeing the [rules broken by useful drugs](https://doi.org/10.1021/acs.jmedchem.8b00686).

So the rules are flexible - but where do you stop? This project attempts to teach an intuitive feel for molecular properties using the old-fashioned trick of looking - if you stare long enough at molecules in a certain property range, you get a feel for what is being included and what is being excluded by a given filtering rule.  

### how
molecules were downloaded from ZINC using the tranche browser - 2d, standard reactivity, and in-stock. ZINC does all the work here - The `.curl` file to do this is in the `./smilesfiles` dir but the files are a bit big for github so you'd have to download them yourself. I used a notebook (`setup.ipynb`) to downsample this to 500,000 molecules (`sample.smifi`). The notebook also uses `rdkit` to calculate molecular weight and logP (using Crippen algorithm). It's worth noting that some say the Crippen algorithm can exaggerate logP (citation needed).

The app itself is a [streamlit](https://www.streamlit.io/) python app. Best to look on their website for instructions on how to make streamlit apps, or inspect `stApp.py`. Streamlit currently (jan 2021) have a service that hosts apps. If that goes down, this app can be run locally using `streamlit run stApp.py`. One thing to note is that `rdkit` is installed by conda, which requires the `conda.txt` instead of a `requirements.txt`, which uses pip, as well as a `conda_channels.txt`. 


### add more properties?
open an issue or drop a line with a property and I'll add it!
