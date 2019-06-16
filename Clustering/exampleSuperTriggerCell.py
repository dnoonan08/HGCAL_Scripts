import uproot
import numpy as np

from superTriggerCellGrouping import superTCMerging, superTCMerging_Fractional

from pyjet import cluster,DTYPE_EP,DTYPE_PTEPM

from mapSuperTC import superTCMap2x2, superTCMap4x4

print "Starting"

#fName = "root://cmseos.fnal.gov//store/user/dnoonan/HGCAL_Concentrator/L1THGCal_Ntuples/Electron_Particle_Gun/ntuple_etaphiTower_pid11_pt50_eta15-30_0PU_0.root"
fName = "root://cmseos.fnal.gov//store/user/dnoonan/HGCAL_Concentrator/L1THGCal_Ntuples/Electron_Particle_Gun_V9/ntuple_hgcalNtuples_pgun_pid11_pt50_eta15-30_200PU_1.root"

#use uproot (available from CMSSW after 10X) to load ntuples from root files into numpy arrays
_tree = uproot.open(fName,xrootdsource=dict(chunkbytes=1024**3, limitbytes=1024**3))["hgcalTriggerNtuplizer/HGCalTriggerNtuple"]

N = _tree.numentries

N = 2

print "File %s"%fName
fulldf = _tree.pandas.df(["tc_subdet","tc_zside","tc_layer","tc_wafer","tc_cell","tc_pt","tc_energy","tc_simenergy","tc_eta",'tc_x','tc_y','tc_z',"tc_mipPt","tc_phi"],entrystart=0, entrystop=N)
    
fulldfGen = _tree.pandas.df(["gen_pt","gen_eta","gen_phi","gen_pdgid","gen_status"],entrystart=0, entrystop=N)

fulldf = superTCMerging(fulldf, useMaxPtLocation=False, equalShare=True, geomVersion="V9", superTCMap = superTCMap2x2)
#fulldf = superTCMerging(fulldf, useMaxPtLocation=True, superTCMap = superTCMap2x2)
    

# for i_event in range(N):
#     print "     %i"%i_event

#     tc = fulldf.loc[i_event,['tc_pt','tc_eta','tc_phi']]

#     tc.columns=["pT","eta","phi"]
#     tc = tc.assign(mass=0.)

#     #go dataframe to np array, formatted for input into fastjet    
#     tcVectors = np.array(tc.to_records(index=False).astype([(u'pT', '<f8'), (u'eta', '<f8'), (u'phi', '<f8'), (u'mass', '<f8')]) )


#     clusterVals = cluster(tcVectors,R=0.4,algo="antikt")
#     _jets = clusterVals.inclusive_jets(ptmin=5)

#     for i,jet in enumerate(_jets):
#         print "Jet %i: pt=%.3f phi=%.3f eta=%.3f"%(i, jet.pt, jet.phi, jet.eta)

