import numpy as np
import pandas as pd

import gc

###https://github.com/PFCal-dev/cmssw/blob/hgc-tpg-devel-CMSSW_11_0_0_pre2/L1Trigger/L1THGCal/data/links_mapping_decentralized_signaldriven_0.txt
linkAllocation = {1  : 2,
                  2  : 2,
                  3  : 2,
                  4  : 2,
                  5  : 2,
                  6  : 2,
                  7  : 4,
                  8  : 2,
                  9  : 5,
                  10 : 2,
                  11 : 4,
                  12 : 2,
                  13 : 3,
                  14 : 2,
                  15 : 2,
                  16 : 2,
                  17 : 2,
                  18 : 2,
                  19 : 2,
                  20 : 2,
                  21 : 2,
                  22 : 2,
                  23 : 2,
                  24 : 2,
                  25 : 2,
                  26 : 2,
                  27 : 2,
                  28 : 2,
                  29 : 2,
                  30 : 2,
                  31 : 2,
                  32 : 2,
                  33 : 2,
                  34 : 2,
                  35 : 2,
                  36 : 2,
                  37 : 2,
                  38 : 2,
                  39 : 2,
                  40 : 2,
                  41 : 2,
                  42 : 2,
                  43 : 2,
                  44 : 2,
                  45 : 2,
                  46 : 2,
                  47 : 2,
                  48 : 2,
                  49 : 2,
                  50 : 2,
                  51 : 2,
                  52 : 2}

###https://github.com/PFCal-dev/cmssw/blob/hgc-tpg-devel-CMSSW_11_0_0_pre2/L1Trigger/L1THGCal/python/hgcalConcentratorProducer_cfi.py#L33-L42
tcPerLink = [0, 1, 3, 6, 9, 14, 18, 23, 27, 32, 37, 41, 46]

tcAllocation = {x : tcPerLink[linkAllocation[x]] for x in linkAllocation}

#remap tc's to ECON-T for scintillators
scintillatorECONMap = {x: int(np.ceil(((x-1)%12+1)/6) + 2*(np.ceil(x/72)-1)) for x in range(1,145)}

def bestChoice_SignalAlloc(fulldf):

    fulldf['TC_Alloc'] = np.where(fulldf.tc_subdet==5,2,fulldf.tc_layer.map(tcAllocation))

    fulldf['ECON-T'] = np.where(fulldf.tc_subdet<5, 1, fulldf.tc_cell.map(scintillatorECONMap))

    fulldf.set_index(['tc_subdet','tc_zside','tc_layer','tc_wafer'],append=True,inplace=True)
    fulldf.reset_index(['subentry'],drop=True,inplace=True)

    ## sort by mipPt
    fulldf.sort_values(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','ECON-T',"tc_mipPt"], ascending=False,inplace=True)

    ## group by wafer, taking only the number of allocated 
    fulldf = fulldf.groupby(['entry','tc_subdet','tc_zside','tc_layer','tc_wafer','ECON-T']).head(fulldf.TC_Alloc)
    
    return fulldf
