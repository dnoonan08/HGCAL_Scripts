import numpy as np
from mapSuperTC import superTCMap2x2

def superTCMerging(fulldf,mergedBranches = ['tc_pt','tc_eta','tc_phi','tc_energy','tc_simenergy','tc_x','tc_y','tc_cell'], useMaxPtLocation=True, superTCMap = None):

    fulldf["tc_superTC"] = np.where(fulldf.tc_subdet==5,fulldf.tc_cell,fulldf['tc_cell'].map(superTCMap))

    columns = ['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC']+mergedBranches
    superTCGroup = fulldf.reset_index().loc[:,columns].groupby(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'])


    if useMaxPtLocation:
        fullDFmax = fulldf.reset_index().sort_values('tc_pt', ascending=False).drop_duplicates(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'])

        fullDFmax.set_index(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'],inplace=True)

        superTC = superTCGroup.sum()
        if 'tc_cell' in mergedBranches: 
            superTC['tc_cell'] = fullDFmax['tc_cell']
        if 'tc_eta' in mergedBranches:
            superTC['tc_eta'] = fullDFmax['tc_eta']
        if 'tc_phi' in mergedBranches:
            superTC['tc_phi'] = fullDFmax['tc_phi']
        if 'tc_x' in mergedBranches:
            superTC['tc_x'] = fullDFmax['tc_x']
        if 'tc_y' in mergedBranches:
            superTC['tc_y'] = fullDFmax['tc_y']
        if 'tc_z' in mergedBranches:
            superTC['tc_z'] = fullDFmax['tc_z']

    else:
        groups = {x:"sum" for x in mergedBranches}
        if 'tc_eta' in mergedBranches:
            groups['tc_eta']="mean"
        if 'tc_phi' in mergedBranches:
            groups['tc_phi']="mean"
        if 'tc_x' in mergedBranches:
            groups['tc_x']="mean"
        if 'tc_y' in mergedBranches:
            groups['tc_y']="mean"
        if 'tc_z' in mergedBranches:
            groups['tc_z']="mean"
        superTC = superTCGroup.agg(groups)

    superTC = superTC[mergedBranches]
    superTC.reset_index(inplace=True)
    superTC.set_index(['entry'],inplace=True)
    return superTC

def superTCMerging_Fractional(fulldf,mergedBranches = ['tc_pt','tc_eta','tc_phi','tc_energy','tc_simenergy','tc_x','tc_y','tc_cell'], useMaxPtLocation=True, superTCMap = None):

    fulldfMerged = superTCMerging(fulldf,mergedBranches = mergedBranches, useMaxPtLocation=useMaxPtLocation, superTCMap = superTCMap)
    fulldf.reset_index(inplace=True)
    fulldf.set_index(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'],inplace=True)
    fulldfMerged.reset_index(inplace=True)
    fulldfMerged.set_index(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'],inplace=True)
    fulldf['supertc_pt'] = fulldfMerged['tc_pt']
    fulldf['tc_frac'] = fulldf['tc_pt']/fulldf['supertc_pt']

    
    columns = ['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC']+mergedBranches+['supertc_pt','tc_frac']
    fulldf.reset_index(inplace=True)
    fulldf = fulldf[columns]
    fulldf.set_index(['entry'],inplace=True)

    return fulldf
