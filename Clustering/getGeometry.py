import uproot
import numpy as np
from mapSuperTC import superTCMap2x2

def getTCGeometry(geomVersion="V9"):
    print "loading trigger cell geometry"
    fName = "root://cmseos.fnal.gov//store/user/dnoonan/HGCAL_Concentrator/triggerGeom%s.root"%geomVersion.upper()

    if not geomVersion.upper() in ['V8','V9-0','V9']:
        raise Exception("The geometry %s is not supported, please use one of the following: ['V8','V9-0','V9']"%geomVersion.upper())

    _tree = uproot.open(fName,xrootdsource=dict(chunkbytes=1024**3, limitbytes=1024**3))["hgcaltriggergeomtester/TreeTriggerCells"]

    geomDF = _tree.pandas.df(['subdet','zside','layer','wafer','triggercell','x','y','z'])
    geomDF['r'] = (geomDF.x**2 + geomDF.y**2)**.5
    geomDF['phi'] = np.arctan(geomDF.y / geomDF.x) + (geomDF.x<0)*((geomDF.y>0)*np.pi - (geomDF.y<0)*np.pi)
    geomDF['eta'] = np.arcsinh(geomDF.z/geomDF.r)

    #add 28 to the layer number for subdets 4 & 5, to match what's in the tc dataframe 
    if geomVersion.upper()=="V9":
        geomDF['layer'] = geomDF.layer + (geomDF.subdet>3)*28
    if geomVersion.upper()=="V9-0":
        geomDF['layer'] = geomDF.layer + (geomDF.subdet>3)*28
    if geomVersion.upper()=="V8":
        geomDF['layer'] = geomDF.layer + (geomDF.subdet==4)*28 + (geomDF.subdet==5)*40

    geomDF = geomDF[['subdet','zside','layer','wafer','triggercell','x','y','z','eta','phi']]
    geomDF.columns =['tc_subdet','tc_zside','tc_layer','tc_wafer','tc_cell','tc_x','tc_y','tc_z','tc_eta','tc_phi']
    geomDF.set_index(['tc_subdet','tc_zside','tc_layer','tc_wafer','tc_cell'],inplace=True)

    return geomDF

def getSuperTCGeometry(geomVersion="V9", superTCMap=superTCMap2x2):
    geomDF = getTCGeometry(geomVersion=geomVersion)
    geomDF.reset_index(inplace=True)
    geomDF['tc_superTC'] = np.where(geomDF.tc_subdet==5,geomDF.tc_cell,geomDF['tc_cell'].map(superTCMap))
    superTCGroup = geomDF.reset_index().groupby(['tc_subdet','tc_zside','tc_layer','tc_wafer','tc_superTC'])
    superTCGeom = superTCGroup.mean()
    return superTCGeom


if __name__=="__main__":
    geomDF = getSuperTCGeometry("V8")
    print geomDF
