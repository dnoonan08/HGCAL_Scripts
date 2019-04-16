import uproot
import numpy as np

import pandas as pd

from pyjet import cluster,DTYPE_EP,DTYPE_PTEPM


print "Starting"

fName = "root://cmseos.fnal.gov//store/user/dnoonan/HGCAL_Concentrator/L1THGCal_Ntuples/Particle_Gun/ntuple_etaphiTower_pid5_pt50_0PU_0.root"

#use uproot (available from CMSSW after 10X) to load ntuples from root files into numpy arrays
_tree = uproot.open(fName,xrootdsource=dict(chunkbytes=1024**3, limitbytes=1024**3))["hgcalTriggerNtuplizer/HGCalTriggerNtuple"]

print "File %s"%fName
fulldf = _tree.pandas.df(["tc_subdet","tc_zside","tc_layer","tc_wafer","tc_cell","tc_pt","tc_energy","tc_simenergy","tc_eta","tc_mipPt","tc_phi"])

N = _tree.numentries

N = 5

for i_event in range(N):
    print "     %i"%i_event

    #create df with subset of entries using only a single event
    tc = fulldf.loc[i_event,['tc_pt','tc_eta','tc_phi']]

    tc.columns=["pT","eta","phi"]
    tc = tc.assign(mass=0.)

    #go dataframe to np array, formatted for input into fastjet    
    tcVectors = np.array(tc.to_records(index=False).astype([(u'pT', '<f8'), (u'eta', '<f8'), (u'phi', '<f8'), (u'mass', '<f8')]) )

    clusterVals = cluster(tcVectors,R=0.4,algo="antikt")
    _jets = clusterVals.inclusive_jets(ptmin=5)

    for i,jet in enumerate(_jets):
        print "Jet %i: pt=%.3f phi=%.3f eta=%.3f"%(i, jet.pt, jet.phi, jet.eta)
