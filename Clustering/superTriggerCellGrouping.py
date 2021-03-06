import numpy as np
import pandas as pd
from mapSuperTC import superTCMap2x2
from getGeometry import getSuperTCGeometry, getSuperTCEqualShareGeometry

def superTCMerging(fulldf,mergedBranches = ['tc_pt','tc_eta','tc_phi','tc_energy','tc_simenergy','tc_x','tc_y','tc_cell'], useMaxPtLocation=True,useMaxELocation=False, equalShare=False, superTCMap = None, geomVersion="V9"):

    fulldf["tc_superTC"] = np.where(fulldf.tc_subdet==5,fulldf.tc_cell,fulldf['tc_cell'].map(superTCMap))

    columns = ['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC']+mergedBranches
    superTCGroup = fulldf.reset_index().loc[:,columns].groupby(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'])


    if useMaxPtLocation or useMaxELocation:
        val = 'tc_pt'
        if useMaxELocation: val = 'tc_energy'

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

    elif equalShare:
        geomDF = getSuperTCEqualShareGeometry(geomVersion=geomVersion, superTCMap=superTCMap)

        superTC = superTCGroup.sum()
        superTC /= 4

        superTCList = []
        eventList = superTC.reset_index().entry.unique()
        eventList.sort()
        for e in eventList:
            geomDF['tc_pt'] = superTC.loc[e].tc_pt
            geomDF['tc_energy'] = superTC.loc[e].tc_energy
            geomDF['tc_simenergy'] = superTC.loc[e].tc_simenergy
            geomDF['entry'] = e
            superTCList.append(geomDF.copy())

        superTC = pd.concat(superTCList)            
        superTC.reset_index(inplace=True)
        superTC.set_index(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'],inplace=True)
        superTC.dropna(inplace=True)

    else:
        geomDF = getSuperTCGeometry(geomVersion=geomVersion, superTCMap=superTCMap)

        superTC = superTCGroup.sum()

        superTC.reset_index(inplace=True)
        superTC.set_index(['tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'],inplace=True)
        if 'tc_eta' in mergedBranches:
            superTC.tc_eta=geomDF.tc_eta
        if 'tc_phi' in mergedBranches:
            superTC.tc_phi=geomDF.tc_phi
        if 'tc_x' in mergedBranches:
            superTC.tc_x=geomDF.tc_x
        if 'tc_y' in mergedBranches:
            superTC.tc_y=geomDF.tc_y
        if 'tc_z' in mergedBranches:
            superTC.tc_z=geomDF.tc_z
        superTC.set_index(['entry'],append=True,inplace=True)
        
        # groups = {x:"sum" for x in mergedBranches}
        # if 'tc_eta' in mergedBranches:
        #     groups['tc_eta']="mean"
        # if 'tc_phi' in mergedBranches:
        #     groups['tc_phi']="mean"
        # if 'tc_x' in mergedBranches:
        #     groups['tc_x']="mean"
        # if 'tc_y' in mergedBranches:
        #     groups['tc_y']="mean"
        # if 'tc_z' in mergedBranches:
        #     groups['tc_z']="mean"
        # superTC = superTCGroup.agg(groups)

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
