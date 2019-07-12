###
# code based on simpleNewVBFreal.py from Uttiya Sarkar
###

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--vbf',help='Run on VBF samples', action='store_true')
parser.add_argument('--particlegun','--pgun',help='Run on Particle Gun samples', action='store_true',dest='particlegun')
parser.add_argument('--pid', default=1, type=int)
parser.add_argument('--pu', default=0, type=int)
parser.add_argument('--pt', default=50, type=int)

parser.add_argument('-N', "-n", default=-1, type=int)
parser.add_argument('--name', default="")
parser.add_argument('--job', default="1/1", help="parallelization jobs, job number over total number of jobs for splitting up files")

parser.add_argument('--verbose','-v',action='store_true',help='Print verbose output of matching')
parser.add_argument('--time',action='store_true',help='Print time information per event')
parser.add_argument('--minpt', default=20., type=float)

parser.add_argument('--genjet',action='store_true',help='Use genJet instead of gen partons')
parser.add_argument('--tcCut', default=-1., type=float,help='mipPT cut to apply to trigger cells')
parser.add_argument('--simEnergyCut', default=-1., type=float,help='simEnergy cut to apply to trigger cells')

parser.add_argument('--draw',action='store_true',help='Draw jet clusters')
parser.add_argument('--V8','--v8',action='store_true',help='Use V8 geometry')

parser.add_argument('--superTC',choices=['2x2','4x4'],default=None,help='Use super trigger cells, specify option of 2x2 or 4x4')
parser.add_argument('--superTCCenter',action='store_true',help='Use super trigger cells with location at center rather than max cell')

args = parser.parse_args()

print args

if args.draw:
    import drawJets

if not (args.vbf or args.particlegun):
    print "Specify vbf or particle gun"
    exit()

import sys
import uproot
import numpy as np
#import pandas as pd
from root_numpy import fill_hist

print "Using uproot version:",uproot.__version__
print "uproot path location:",uproot.__file__

import gc

from subprocess import Popen,PIPE

from pyjet import cluster,DTYPE_EP,DTYPE_PTEPM
#import histvbf as h
#from histIso import *

print "Starting"

#count number of gen particles, matched, and unmatched jets
nGenParticlesinHGCAL = 0
nmatchedJets = 0
nunmatchedJets = 0

#timers
import time
if args.time:
    start = time.clock()

skipFill = False

fileRage = []
fileName = ""

if args.V8:
    if args.particlegun:
        eosDir = "/store/user/dnoonan/HGCAL_Concentrator/L1THGCal_Ntuples/Particle_Gun/"
        if args.pid==11:
            eosDir = "/store/user/dnoonan/HGCAL_Concentrator/L1THGCal_Ntuples/Electron_Particle_Gun/"
        fileNameContent = "ntuple_etaphiTower_pid%i_pt%i_eta15-30_%iPU"%(args.pid, args.pt, args.pu)
        outputFileName = "particleGun_pid%i_pt%i_eta15-30_%iPU.root"%(args.pid, args.pt, args.pu)

    if args.vbf:
        eosDir = "/store/user/dnoonan/HGCAL_Concentrator/L1THGCal_Ntuples/VBF_Samples/"
        fileNameContent = "ntuple_hgcalNtuples_vbf_hmm_%iPU"%(args.pu)
        outputFileName = "VBF_%iPU.root"%(args.pu)
else:
    eosDir = "/store/user/dnoonan/HGCAL_Concentrator/L1THGCal_Ntuples/Electron_Particle_Gun_V9/"
    if args.particlegun:
        fileNameContent = "ntuple_hgcalNtuples_pgun_pid%i_pt%i_eta15-30_%iPU"%(args.pid, args.pt, args.pu)
        outputFileName = "particleGun_pid%i_pt%i_eta15-30_%iPU.root"%(args.pid, args.pt, args.pu)
    if args.vbf:
        fileNameContent = "ntuple_hgcalNtuples_vbf_hmm_%iPU"%(args.pu)
        outputFileName = "VBF_%iPU.root"%(args.pu)


if args.tcCut>-1:
    print "HERE"
    outputFileName = outputFileName.replace(".root","_tcMipPt_gt_%i.root"%int(args.tcCut))

if args.simEnergyCut>-1:
    outputFileName = outputFileName.replace(".root","_simEnergy_gt_%i.root"%int(args.simEnergy))

if args.V8:
    outputFileName = outputFileName.replace(".root","_geomV8.root")

if args.superTC:
    outputFileName = outputFileName.replace(".root","_superTC%s.root"%args.superTC)

    if args.superTCCenter:
        outputFileName = outputFileName.replace(".root","Centered.root")

if not args.name=="":
    outputFileName = outputFileName.replace(".root","_%s.root"%args.name)


fileList = []
# get list of files
cmd = "xrdfs root://cmseos.fnal.gov ls %s"%eosDir
dirContents,stderr = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE).communicate()
dirContentsList = dirContents.split("\n")
for fName in dirContentsList:
    if fileNameContent in fName:
        fileList.append("root://cmseos.fnal.gov/%s"%(fName))

if "/" in args.job:
    if not args.job=='1/1': 
        print "Running job",args.job

        jobN = int(args.job.split('/')[0])
        nJobs = int(args.job.split('/')[1])
        totalFiles = len(fileList)
        filesPerJob = int(round(totalFiles/nJobs+.4999))

        fileList = fileList[(jobN-1)*filesPerJob:jobN*filesPerJob]
        outputFileName = outputFileName.replace(".root","_%iof%i.root"%(jobN, nJobs))

from makeIsoTree import *

output = IsoOutputTree(outputFileName)

totalN = 0       

# if args.NFiles==-1:
#     nFiles = len(fileList)
# else:
#     nFiles = args.NFiles

for fName in fileList:

    if not args.N==-1 and (totalN > args.N):
        break

    print "File %s"%fName    
    try:
        _tree = uproot.open(fName,xrootdsource=dict(chunkbytes=1024**3, limitbytes=1024**3))["hgcalTriggerNtuplizer/HGCalTriggerNtuple"]
    except:
        print "---Unable to open file, skipping"
        continue

    if args.N==-1:
        N = _tree.numentries
    else:
        if totalN > args.N:
            break
        nRemaining = args.N - totalN
        if _tree.numentries < nRemaining:
            N = _tree.numentries
        else:
            N = nRemaining
        totalN += _tree.numentries


    branches = ["tc_pt","tc_energy","tc_eta","tc_mipPt","tc_simenergy","tc_phi"]
    if args.superTC:
        branches = ["tc_subdet","tc_zside","tc_layer","tc_wafer","tc_cell"] + branches


    fulldf = _tree.pandas.df(branches,entrystart=0,entrystop=N)

    
    if args.superTC:
        from mapSuperTC import superTCMap2x2, superTCMap4x4
        from superTriggerCellGrouping import superTCMerging
        if args.superTC=='2x2':
            superTCMap = superTCMap2x2
        if args.superTC=='4x4':
            superTCMap = superTCMap4x4
        geomVersion = "V9"
        if args.V8:
            geomVersion = "V8"
        useMaxPtLocation = not args.superTCCenter

        fulldf = superTCMerging(fulldf, mergedBranches=['tc_pt','tc_energy','tc_eta','tc_phi','tc_mipPt','tc_simenergy'], useMaxPtLocation=useMaxPtLocation, geomVersion=geomVersion, superTCMap = superTCMap)



    fulldf = fulldf[(fulldf.tc_mipPt > args.tcCut) & (fulldf.tc_simenergy > args.simEnergyCut)]
    gc.collect()

    fulldfGen = _tree.pandas.df(["gen_pt","gen_eta","gen_phi","gen_status","gen_pdgid","gen_energy"],entrystart=0,entrystop=N)
    fulldfGenJet = _tree.pandas.df(["genjet_pt","genjet_eta","genjet_phi","genjet_energy"],entrystart=0,entrystop=N)

    del _tree
    gc.collect()

    print "Loaded Tree"
    if args.time:
        print "Loaded tree ", time.clock()-start




    for i_event in range(N):
        if args.time: 
            start = time.clock()
            print '-'*20
            print "     %i"%i_event   

        # load into dictionary for creating pandas df
        tc = fulldf.loc[i_event,['tc_pt','tc_eta','tc_phi','tc_energy']]
        tc.columns=["pT","eta","phi","energy"]
        tc = tc.assign(mass=0)

        #go dataframe to np array, formatted for input into fastjet    
        # tcVectors = np.array(tc.to_records(index=False).astype([(u'pT', '<f8'), (u'eta', '<f8'), (u'phi', '<f8'), (u'energy', '<f8')]) ) 
        tcVectors = np.array(tc[["pT","eta","phi","mass"]].to_records(index=False).astype([(u'pT', '<f8'), (u'eta', '<f8'), (u'phi', '<f8'), (u'mass', '<f8')]) ) 

        if args.time:
            print "    Loaded in", time.clock()-start
        clusterVals = cluster(tcVectors,R=0.4,algo="antikt")
        _jets = clusterVals.inclusive_jets(ptmin=args.minpt)
        if args.time:
            print "    Clustered jets in", time.clock()-start
            print "      --- len =",len(tcVectors), (time.clock()-start)/len(tcVectors)
        del tcVectors
        del clusterVals
        gc.collect()

        genDF = fulldfGen.loc[i_event,['gen_pt','gen_eta','gen_phi','gen_status','gen_pdgid','gen_energy']]
        if args.vbf:
            genDF = genDF[(genDF.gen_status==23) & (abs(genDF.gen_pdgid)<6) & (abs(genDF.gen_eta)>1.5) & (abs(genDF.gen_eta)<3.)]
        else:
            genDF = genDF[((genDF.gen_status==1) | (genDF.gen_status==23)) & (abs(genDF.gen_pdgid)==args.pid) & (abs(genDF.gen_eta)>1.5) & (abs(genDF.gen_eta)<3.)]
        

        genJetDF = fulldfGenJet.loc[i_event,['genjet_pt','genjet_eta','genjet_phi','genjet_energy']]
        genJetDF.columns = ['gen_pt','gen_eta','gen_phi','gen_energy']
        genJetDF = genJetDF[(abs(genJetDF.gen_eta)>1.5) & (abs(genJetDF.gen_eta)<3.) & (genJetDF.gen_pt>10)]
        # genJetDF = genJetDF[(abs(genJetDF.gen_eta)>1.8) & (abs(genJetDF.gen_eta)<2.7) & (genJetDF.gen_pt>10)]
            
        output.nGen[0] = len(genDF)
        output.nGenJet[0] = len(genJetDF)
        for i in range(output.nGen[0]):
            output.genPt[i] = genDF.gen_pt.values[i]
            output.genEta[i] = genDF.gen_eta.values[i]
            output.genPhi[i] = genDF.gen_phi.values[i]
            output.genEn[i] = genDF.gen_energy.values[i]
            output.genPID[i] = genDF.gen_pdgid.values[i]

        for i in range(output.nGenJet[0]):
            output.genJetPt[i] = genJetDF.gen_pt.values[i]
            output.genJetEta[i] = genJetDF.gen_eta.values[i]
            output.genJetPhi[i] = genJetDF.gen_phi.values[i]
            output.genJetEn[i] = genJetDF.gen_energy.values[i]


        if args.time:
            print "    Got Gen in", time.clock()-start


        isGenJetMatched=[]
        isGenPartMatched=[]
        genjetVector = genJetDF.values

        for j in range(len(_jets)):
            output.jetGenMatch[j] = -1
            output.jetGenJetMatch[j] = -1
            output.jetMinGenDR[j] = 99
            output.jetMinGenJetDR[j] = 99

	genVector = genDF.values
        if args.verbose:
            print 'Gen Particles'
            print genDF
            print '------'

        if args.verbose:
            print 'Gen Jets Before Cleaning'
            print genJetDF
            print '------'
        #gen part to gen jet matching
        # remove gen jets which are not matched to a gen parton (within dR of at least 0.3)
        for k in range(len(genJetDF)-1,-1,-1):
            minDR = 999
            for i in range(len(genVector)):
                dR_gen_genjet=((genjetVector[k,1]-genVector[i,1])**2+(genjetVector[k,2]-genVector[i,2])**2)**0.5
                if dR_gen_genjet < minDR:
                    minDR = dR_gen_genjet
            if minDR>0.3:
                genJetDF.drop(genJetDF.index.values[k],inplace=True)
        genjetVector = genJetDF.values

        if args.verbose:
            print 'Gen Jets'
            print genJetDF
            print '------'



        if args.verbose:
            print 'Reco FastJets'
            for jet in _jets:
                print '   ', jet
            print '------'

        #print genjetVector
        if args.verbose:
            print 'Jet/GenJet Matching'

        for i in range(len(genjetVector)):
            foundMatch = False

            if args.verbose:
                print genjetVector[i][:4]
            for j,jet in enumerate(_jets):
                 dR_jet_genjet=((jet.eta-genjetVector[i,1])**2+(jet.phi-genjetVector[i,2])**2)**0.5

                 if dR_jet_genjet < output.jetMinGenJetDR[j]:
                     output.jetMinGenJetDR[j] = dR_jet_genjet

                 if args.verbose:
                     print "  --- %i %i %.5f\t%+.4f\t%+.4f\t%.5f"%(i, j, jet.pt, jet.eta, jet.phi, dR_jet_genjet)

                 if dR_jet_genjet<0.1:
                     output.jetGenJetMatch[j] = i
                     foundMatch = True
                     if args.verbose:
                         print 'Match'
                     break   

            isGenJetMatched.append(foundMatch)
            if args.verbose:
                if not foundMatch:
                    print "No Match Found"

        isGenMatched = []
         

#        genVector = genDF.values
#        if args.verbose:
#            print 'Gen Particles'

        if args.verbose:
            print 'Jet/GenPart Matching'
        for i in range(len(genVector)):
            foundMatch = False
            if args.verbose:
                print genVector[i][:4]

            for j,jet in enumerate(_jets):
                 dR_jet_gen=((jet.eta-genVector[i,1])**2+(jet.phi-genVector[i,2])**2)**0.5

                 if dR_jet_gen < output.jetMinGenDR[j]:
                     output.jetMinGenDR[j] = dR_jet_gen

                 if args.verbose:
                     print "  --- %i %i %.5f\t%+.4f\t%+.4f\t%.5f"%(i, j, jet.pt, jet.eta, jet.phi, dR_jet_gen)

                 if dR_jet_gen<0.1:
                     output.jetGenMatch[j] = i
                     foundMatch = True
                     break   

            isGenMatched.append(foundMatch)
            if args.verbose:
                if not foundMatch:
                    print "No Match Found"



        if args.time:
            print "    GenMatch in", time.clock()-start


        if args.draw:
            figureName = "Plots/JetAreas_%s"%(outputFileName.replace(".root","_%i.pdf"%i_event))
            drawJets.drawJets(jets = _jets,name = figureName ,genjet = genJetDF, genjet_matched=isGenJetMatched,genpart = genDF, genpart_matched=isGenMatched)
            # if args.genjet:
            #     drawJets.drawJets(jets = _jets,name = figureName ,genjet_pt = genjetDF.gen_pt.values, genjet_eta = genjetDF.gen_eta.values, genjet_phi = genjetDF.gen_phi.values, genjet_matched=isGenJetMatched)
            # else:
            #     drawJets.drawJets(jets = _jets,name = figureName ,genjet_pt = genDF.gen_pt.values, genjet_eta = genDF.gen_eta.values, genjet_phi = genDF.gen_phi.values, genjet_pid = genDF.gen_pdgid.values, genjet_matched=isGenMatched)


        output.nJet[0] = len(_jets)

        for j,jet in enumerate(_jets):

            tc['dEta'] = tc['eta']-jet.eta
            tc['dPhi'] = tc['phi']-jet.phi
            tc['dR'] = (tc.dEta**2 + tc.dPhi**2)**0.5
            
            cut = tc.dR<0.6
            tcThisJet = tc[cut]
            
            A = tcThisJet.loc[tcThisJet.dR<0.1,'pT'].sum()
            B = tcThisJet.loc[(tcThisJet.dR>0.1) & (tcThisJet.dR<0.2),'pT'].sum()
            C = tcThisJet.loc[(tcThisJet.dR>0.2) & (tcThisJet.dR<0.4),'pT'].sum()

            iso = (B-(3./12)*C)/(A-(1./12)*C)
            
            output.jetPt[j] = jet.pt
            output.jetEta[j] = jet.eta
            output.jetPhi[j] = jet.phi

            output.jetPtR01[j] = tcThisJet.loc[(tcThisJet.dR>=0.0) & (tcThisJet.dR<0.1),'pT'].sum()
            output.jetPtR02[j] = tcThisJet.loc[(tcThisJet.dR>=0.1) & (tcThisJet.dR<0.2),'pT'].sum()
            output.jetPtR03[j] = tcThisJet.loc[(tcThisJet.dR>=0.2) & (tcThisJet.dR<0.3),'pT'].sum()
            output.jetPtR04[j] = tcThisJet.loc[(tcThisJet.dR>=0.3) & (tcThisJet.dR<0.4),'pT'].sum()
            output.jetPtR05[j] = tcThisJet.loc[(tcThisJet.dR>=0.4) & (tcThisJet.dR<0.5),'pT'].sum()
            output.jetPtR06[j] = tcThisJet.loc[(tcThisJet.dR>=0.5) & (tcThisJet.dR<0.6),'pT'].sum()

            output.jetEnR01[j] = tcThisJet.loc[(tcThisJet.dR>=0.0) & (tcThisJet.dR<0.1),'energy'].sum()
            output.jetEnR02[j] = tcThisJet.loc[(tcThisJet.dR>=0.1) & (tcThisJet.dR<0.2),'energy'].sum()
            output.jetEnR03[j] = tcThisJet.loc[(tcThisJet.dR>=0.2) & (tcThisJet.dR<0.3),'energy'].sum()
            output.jetEnR04[j] = tcThisJet.loc[(tcThisJet.dR>=0.3) & (tcThisJet.dR<0.4),'energy'].sum()
            output.jetEnR05[j] = tcThisJet.loc[(tcThisJet.dR>=0.4) & (tcThisJet.dR<0.5),'energy'].sum()
            output.jetEnR06[j] = tcThisJet.loc[(tcThisJet.dR>=0.5) & (tcThisJet.dR<0.6),'energy'].sum()

            output.jetIso[j] = iso


        output.tree.Fill()

        if args.time:
            print "    process UnMatched in", time.clock()-start


        del tc
        del genVector
        del genjetVector

        gc.collect()

        
    
    del fulldf
    del fulldfGen
    gc.collect()

    sys.stdout.flush()


output.outFile.Write()
