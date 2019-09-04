import numpy as np
import pandas as pd
from mapSuperTC import superTCMap2x2, superTCMapFrom2x2To2x4
from scintillatorMapping import scintillatorMapToHGROC
from getGeometry import getTCGeometry, getSuperTCGeometry, getSuperTCEqualShareGeometry


def superTCMerging_CTC8_LDM(fulldf,mergedBranches = ['tc_pt','tc_eta','tc_phi','tc_energy','tc_simenergy','tc_x','tc_y','tc_cell'],geomVersion='V9'):

    geomDF = getTCGeometry(geomVersion=geomVersion).reset_index().drop_duplicates(['tc_subdet','tc_zside','tc_layer','tc_wafer']).set_index(['tc_subdet','tc_zside','tc_layer','tc_wafer'])[['HDM']]

    fulldf.reset_index(inplace=True)
    fulldf.set_index(['tc_subdet','tc_zside','tc_layer','tc_wafer'],inplace=True)


    fulldf['HDM'] = geomDF.HDM

    fulldf.loc[5,'HDM']=True

    fulldf.reset_index(inplace=True)
    fulldf.set_index(['entry','subentry'])

    fulldf["tc_superTC"] = np.where(fulldf.tc_subdet==5,fulldf['tc_cell'].map(scintillatorMapToHGROC),
                                    np.where(fulldf.HDM, fulldf['tc_cell'].map(superTCMap2x2), fulldf['tc_cell'].map(superTCMap2x2)))

    columns = ['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC','HDM']+mergedBranches
    superTCGroup = fulldf.reset_index().loc[:,columns].groupby(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'])

#    geomDF = getSuperTCEqualShareGeometry(geomVersion=geomVersion, superTCMap=superTCMap2x2)

    val = 'tc_energy'
    
    fullDFmax = fulldf.reset_index().sort_values(val, ascending=False).drop_duplicates(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'])

    fullDFmax.set_index(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'],inplace=True)

    superTC = superTCGroup.sum()
    superTC['HDM'] = superTCGroup.HDM.any()


    #use max energy TC for high density modules
    if 'tc_cell' in mergedBranches: 
        superTC.loc[fullDFmax.HDM,'tc_cell'] = fullDFmax.loc[fullDFmax.HDM,'tc_cell']
    if 'tc_eta' in mergedBranches:
        superTC.loc[fullDFmax.HDM,'tc_eta'] = fullDFmax.loc[fullDFmax.HDM,'tc_eta']
    if 'tc_phi' in mergedBranches:
        superTC.loc[fullDFmax.HDM,'tc_phi'] = fullDFmax.loc[fullDFmax.HDM,'tc_phi']
    if 'tc_x' in mergedBranches:
        superTC.loc[fullDFmax.HDM,'tc_x'] = fullDFmax.loc[fullDFmax.HDM,'tc_x']
    if 'tc_y' in mergedBranches:
        superTC.loc[fullDFmax.HDM,'tc_y'] = fullDFmax.loc[fullDFmax.HDM,'tc_y']
    if 'tc_z' in mergedBranches:
        superTC.loc[fullDFmax.HDM,'tc_z'] = fullDFmax.loc[fullDFmax.HDM,'tc_z']


    geomDF_STC = getSuperTCGeometry(geomVersion=geomVersion, superTCMapHDM=superTCMap2x2, superTCMapLDM=superTCMap2x2,superTCMapScint = scintillatorMapToHGROC)

    superTC.reset_index(inplace=True)
    superTC.set_index(['tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'],inplace=True)
    if 'tc_eta' in mergedBranches:
        superTC.loc[np.invert(superTC.HDM),"tc_eta"] = geomDF_STC.loc[np.invert(geomDF_STC.HDM),"tc_eta"]
    if 'tc_phi' in mergedBranches:
        superTC.loc[np.invert(superTC.HDM),"tc_phi"] = geomDF_STC.loc[np.invert(geomDF_STC.HDM),"tc_phi"]
    if 'tc_x' in mergedBranches:
        superTC.loc[np.invert(superTC.HDM),"tc_x"] = geomDF_STC.loc[np.invert(geomDF_STC.HDM),"tc_x"]
    if 'tc_y' in mergedBranches:
        superTC.loc[np.invert(superTC.HDM),"tc_y"] = geomDF_STC.loc[np.invert(geomDF_STC.HDM),"tc_y"]
    if 'tc_z' in mergedBranches:
        superTC.loc[np.invert(superTC.HDM),"tc_z"] = geomDF_STC.loc[np.invert(geomDF_STC.HDM),"tc_z"]

    superTC.reset_index(inplace=True)
    superTC.set_index(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'],inplace=True)

    superTC.reset_index(inplace=True)
    superTC.set_index(['entry'],inplace=True)
    superTC['tc_superTC'] = np.where(((superTC.HDM)), superTC.tc_superTC,
                                     superTC.tc_superTC.map(superTCMapFrom2x2To2x4)
                                 )

    superTCGroupFinal = superTC.reset_index().groupby(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'])
    superTCmax = superTC.reset_index().sort_values(val, ascending=False).drop_duplicates(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'])
    superTCmax.set_index(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'],inplace=True)

    superTCFinal = superTCGroupFinal.sum()

    if 'tc_cell' in mergedBranches:
        superTCFinal['tc_cell'] = superTCmax['tc_cell']
    if 'tc_eta' in mergedBranches:
        superTCFinal['tc_eta'] = superTCmax['tc_eta']
    if 'tc_phi' in mergedBranches:
        superTCFinal['tc_phi'] = superTCmax['tc_phi']
    if 'tc_x' in mergedBranches:
        superTCFinal['tc_x'] = superTCmax['tc_x']
    if 'tc_y' in mergedBranches:
        superTCFinal['tc_y'] = superTCmax['tc_y']
    if 'tc_z' in mergedBranches:
        superTCFinal['tc_z'] = superTCmax['tc_z']

    superTC = superTCFinal[mergedBranches+['HDM']]
    superTC.reset_index(inplace=True)
    superTC.set_index(['entry'],inplace=True)

    return superTC


def superTCMerging_STC16_LDM(fulldf,mergedBranches = ['tc_pt','tc_eta','tc_phi','tc_energy','tc_simenergy','tc_x','tc_y','tc_cell'],geomVersion='V9'):

    geomDF = getTCGeometry(geomVersion=geomVersion).reset_index().drop_duplicates(['tc_subdet','tc_zside','tc_layer','tc_wafer']).set_index(['tc_subdet','tc_zside','tc_layer','tc_wafer'])[['HDM']]

    fulldf.reset_index(inplace=True)
    fulldf.set_index(['tc_subdet','tc_zside','tc_layer','tc_wafer'],inplace=True)

    fulldf['HDM'] = geomDF.HDM
    fulldf.loc[5,'HDM']=True

    fulldf.reset_index(inplace=True)
    fulldf.set_index(['entry','subentry'])

    fulldf["tc_superTC"] = np.where(fulldf.tc_subdet==5,fulldf['tc_cell'].map(scintillatorMapToHGROC),
                                    np.where(fulldf.HDM, fulldf['tc_cell'].map(superTCMap2x2), fulldf['tc_cell'].map(superTCMap2x2)))

    columns = ['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC','HDM']+mergedBranches
    superTCGroup = fulldf.reset_index().loc[:,columns].groupby(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'])

#    geomDF = getSuperTCEqualShareGeometry(geomVersion=geomVersion, superTCMap=superTCMap)

    val = 'tc_energy'
    
    fullDFmax = fulldf.reset_index().sort_values(val, ascending=False).drop_duplicates(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'])

    fullDFmax.set_index(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'],inplace=True)

    superTC = superTCGroup.sum()
    superTC['HDM'] = superTCGroup.HDM.any()

    #use max energy TC for high density modules
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


    superTC = superTC[mergedBranches+['HDM']]
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
