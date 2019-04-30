###
# code based on simpleNewVBFreal.py from Uttiya Sarkar
###

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--vbf',help='Run on VBF samples', action='store_true')
parser.add_argument('--particlegun','--pgun',help='Run on Particle Gun samples', action='store_true',dest='particlegun')
parser.add_argument('--pid', default=1, type=int)
parser.add_argument('--pu', default=140, type=int)
parser.add_argument('--pt', default=50, type=int)

parser.add_argument('-N', "-n", default=-1, type=int)
#parser.add_argument('--NFiles', "--nFiles", default=-1, type=int)
parser.add_argument('--name', default="")
parser.add_argument('--job', default="1/1", help="parallelization jobs, job number over total number of jobs for splitting up files")

parser.add_argument('--verbose','-v',action='store_true',help='Print verbose output of matching')
parser.add_argument('--time',action='store_true',help='Print time information per event')
parser.add_argument('--minpt', default=20., type=float)

parser.add_argument('--genjet',action='store_true',help='Use genJet instead of gen partons')
parser.add_argument('--tcCut', default=-1., type=float,help='mipPT cut to apply to trigger cells')
parser.add_argument('--simEnergyCut', default=-1., type=float,help='simEnergy cut to apply to trigger cells')

parser.add_argument('--draw',action='store_true',help='Draw jet clusters')
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
from histIso import *

print "Starting"

#count number of gen particles, matched, and unmatched jets
nGenParticlesinHGCAL = 0
nmatchedJets = 0
nunmatchedJets = 0

#timers
if args.time:
    import time
    start = time.clock()

skipFill = False

fileRage = []
fileName = ""


if args.particlegun:
    eosDir = "/store/user/dnoonan/HGCAL_Concentrator/L1THGCal_Ntuples/Particle_Gun/"
    fileNameContent = "ntuple_etaphiTower_pid%i_pt%i_eta15-30_%iPU"%(args.pid, args.pt, args.pu)
    outputFileName = "particleGun_pid%i_pt%i_eta15-30_%iPU.root"%(args.pid, args.pt, args.pu)

if args.vbf:
    eosDir = "/store/user/dnoonan/HGCAL_Concentrator/L1THGCal_Ntuples/VBF_Samples/"
    fileNameContent = "ntuple_hgcalNtuples_vbf_hmm_%iPU"%(args.pu)
    outputFileName = "VBF_%iPU.root"%(args.pu)

if args.tcCut>-1:
    print "HERE"
    outputFileName = outputFileName.replace(".root","_tcMipPt_gt_%i.root"%int(args.tcCut))

if args.simEnergyCut>-1:
    outputFileName = outputFileName.replace(".root","_simEnergy_gt_%i.root"%int(args.simEnergy))

if args.genjet:
    outputFileName = outputFileName.replace(".root","_matchToGenJets.root")
else:
    outputFileName = outputFileName.replace(".root","_matchToGenPart.root")

if not args.name=="":
    outputFileName = outputFileName.replace(".root","_%s.root"%args.name)


fileList = []

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

h = MakeHist(outputFileName)

totalN = 0       
# if args.NFiles==-1:
#     nFiles = len(fileList)
# else:
#     nFiles = args.NFiles

for fName in fileList:

    if not args.N==-1 and (totalN > args.N):
        break

    _tree = uproot.open(fName,xrootdsource=dict(chunkbytes=1024**3, limitbytes=1024**3))["hgcalTriggerNtuplizer/HGCalTriggerNtuple"]

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

    print "File %s"%fName
    fulldf = _tree.pandas.df(["tc_pt","tc_energy","tc_eta","tc_mipPt","tc_simenergy","tc_phi"],entrystart=0,entrystop=N)
    fulldf = fulldf[(fulldf.tc_mipPt > args.tcCut) & (fulldf.tc_simenergy > args.simEnergyCut)]
    gc.collect()
    # if args.tcCut>-1:
    #     fulldf = fulldf[fulldf.tc_mipPt > args.tcCut]
    # if args.simEnergyCut>-1:
    #     fulldf = fulldf[fulldf.tc_simenergy > args.simEnergyCut]
    if not args.genjet:
        fulldfGen = _tree.pandas.df(["gen_pt","gen_eta","gen_phi","gen_status","gen_pdgid","gen_energy"],entrystart=0,entrystop=N)
    else:
        fulldfGen = _tree.pandas.df(["genjet_pt","genjet_eta","genjet_phi","genjet_energy"],entrystart=0,entrystop=N)

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

        if not args.genjet:
            genjetDF = fulldfGen.loc[i_event,['gen_pt','gen_eta','gen_phi','gen_status','gen_pdgid','gen_energy']]
            if args.vbf:
                genjetDF = genjetDF[(genjetDF.gen_status==23) & (abs(genjetDF.gen_pdgid)<6) & (abs(genjetDF.gen_eta)>1.5) & (abs(genjetDF.gen_eta)<3.)]
            else:
                genjetDF = genjetDF[(genjetDF.gen_status==23) & (abs(genjetDF.gen_pdgid)==args.pid) & (abs(genjetDF.gen_eta)>1.5) & (abs(genjetDF.gen_eta)<3.)]
        else:
            genjetDF = fulldfGen.loc[i_event,['genjet_pt','genjet_eta','genjet_phi','genjet_energy']]
            genjetDF.columns = ['gen_pt','gen_eta','gen_phi','gen_energy']
            genjetDF = genjetDF[(abs(genjetDF.gen_eta)>1.5) & (abs(genjetDF.gen_eta)<3.) & (genjetDF.gen_pt>10)]
            
            

        fill_hist(h.genJetPt   ,genjetDF.gen_pt.values)

        genjetVector = genjetDF.values
        
        nGenParticlesinHGCAL += len(genjetVector)
        if args.time:
            print "    Got Gen in", time.clock()-start

        matchedJets = []
        genMatch = []
        isGenMatched=[]

        for i in range(len(genjetVector)):
            foundMatch = False
            if args.verbose:
                print genjetVector[i][:4]
            for j in range(len(_jets)):
                 dRSq_jet_genjet=((_jets[j].eta-genjetVector[i,1])**2+(_jets[j].phi-genjetVector[i,2])**2)

                 if args.verbose:
                     print "  --- %i %i %.5f\t%+.4f\t%+.4f\t%.5f"%(i, j, _jets[j].pt, _jets[j].eta, _jets[j].phi, dRSq_jet_genjet)

                 if dRSq_jet_genjet<0.01:
                    matchedJets.append(j)
                    genMatch.append(i)
                    foundMatch = True
                    break   
            isGenMatched.append(foundMatch)
            if args.verbose:
                if not foundMatch:
                    print "No Match Found"

        if args.time:
            print "    GenMatch in", time.clock()-start
            print "       Matched %i jets"%len(matchedJets)


        if args.draw:
            figureName = "Plots/JetAreas_%s"%(outputFileName.replace(".root","_%i.pdf"%i_event))
            if args.genjet:
                drawJets.drawJets(jets = _jets,name = figureName ,genjet_pt = genjetDF.gen_pt.values, genjet_eta = genjetDF.gen_eta.values, genjet_phi = genjetDF.gen_phi.values, genjet_matched=isGenMatched)
            else:
                drawJets.drawJets(jets = _jets,name = figureName ,genjet_pt = genjetDF.gen_pt.values, genjet_eta = genjetDF.gen_eta.values, genjet_phi = genjetDF.gen_phi.values, genjet_pid = genjetDF.gen_pdgid.values, genjet_matched=isGenMatched)

        unmatchedJets = [x for x in range(len(_jets)) if x not in matchedJets]
        unmatchedGen = [x for x in range(len(genjetVector)) if x not in genMatch]

        for i in range(len(matchedJets)):
            j = matchedJets[i]

            tc['dEta'] = tc['eta']-_jets[j].eta
            tc['dPhi'] = tc['phi']-_jets[j].phi
            tc['dR'] = (tc.dEta**2 + tc.dPhi**2)**0.5
            
            cut = tc.dR<0.45
            tcThisJet = tc[cut]
            
            A = tcThisJet.loc[tcThisJet.dR<0.1,'pT'].sum()
            B = tcThisJet.loc[(tcThisJet.dR>0.1) & (tcThisJet.dR<0.2),'pT'].sum()
            C = tcThisJet.loc[(tcThisJet.dR>0.2) & (tcThisJet.dR<0.4),'pT'].sum()

            isoowen = (B-(3./12)*C)/(A-(1./12)*C)

            if skipFill: continue

            h.genJetPtMatched.Fill(genjetVector[genMatch[i],0])
            h.genVsClusterPt.Fill(genjetVector[genMatch[i],0], _jets[j].pt)
            h.genVsClusterPt_puSub.Fill(genjetVector[genMatch[i],0], _jets[j].pt - C/.75)

            h.jetPt.Fill(_jets[j].pt)
            h.jetPt_puSub.Fill(_jets[j].pt - C/.75)

            h.iso_PUsub_all.Fill(isoowen)
            h.iso_noPUsub_all.Fill(B/A)

            h.genPtVsClusterIso.Fill(genjetVector[genMatch[i],0],isoowen)
            h.clusterPtVsClusterIso.Fill(_jets[j].pt,isoowen)

            fill_hist(h.detadphiall   ,tcThisJet[['dEta','dPhi']].values,tcThisJet.pT.values)
            fill_hist(h.dRgall        ,tcThisJet[['dR','energy']].values)
            fill_hist(h.dRptall       ,tcThisJet[['dR','pT']].values)
            fill_hist(h.tcpt          ,tcThisJet.pT.values)
            fill_hist(h.tcenergy      ,tcThisJet.energy.values)
            fill_hist(h.DRe           ,tcThisJet.dR.values,tcThisJet.pT.values)
            fill_hist(h.DRp           ,tcThisJet.dR.values,tcThisJet.energy.values)
            
            if abs(_jets[j].eta)>1.5 and abs(_jets[j].eta)<1.9:
                h.iso_PUsub_eta15.Fill(isoowen)
                h.iso_noPUsub_eta15.Fill(B/A)
                fill_hist(h.detadphieta15        ,tcThisJet[['dEta','dPhi']].values, tcThisJet.pT.values)
                fill_hist(h.dReta15              ,tcThisJet[['dR','pT']].values)

            elif abs(_jets[j].eta)>1.9 and abs(_jets[j].eta)<2.3:
                h.iso_PUsub_eta20.Fill(isoowen)
                h.iso_noPUsub_eta20.Fill(B/A)
                fill_hist(h.detadphieta20        ,tcThisJet[['dEta','dPhi']].values, tcThisJet.pT.values)
                fill_hist(h.dReta20              ,tcThisJet[['dR','pT']].values)

            elif abs(_jets[j].eta)>2.3 and abs(_jets[j].eta)<2.7:
                h.iso_PUsub_eta25.Fill(isoowen)
                h.iso_noPUsub_eta25.Fill(B/A)
                fill_hist(h.detadphieta25        ,tcThisJet[['dEta','dPhi']].values, tcThisJet.pT.values)
                fill_hist(h.dReta25              ,tcThisJet[['dR','pT']].values)

            elif abs(_jets[j].eta)>2.7:
                h.iso_PUsub_eta30.Fill(isoowen)
                h.iso_noPUsub_eta30.Fill(B/A)
                fill_hist(h.detadphieta30        ,tcThisJet[['dEta','dPhi']].values, tcThisJet.pT.values)
                fill_hist(h.dReta30              ,tcThisJet[['dR','pT']].values)

        if args.time:
            print "    process Matched in", time.clock()-start
                       
        for j in unmatchedJets:
            tc['dEta'] = tc['eta']-_jets[j].eta
            tc['dPhi'] = tc['phi']-_jets[j].phi
            tc['dR'] = (tc.dEta**2 + tc.dPhi**2)**0.5

            cut = tc.dR<0.45
            tcThisJet = tc[cut]

            A = tcThisJet.loc[tcThisJet.dR<0.1,'pT'].sum()
            B = tcThisJet.loc[(tcThisJet.dR>0.1) & (tcThisJet.dR<0.2),'pT'].sum()
            C = tcThisJet.loc[(tcThisJet.dR>0.2) & (tcThisJet.dR<0.4),'pT'].sum()

            isoowenpu = (B-(3./12)*C)/(A-(1./12)*C)


            if skipFill: continue

            h.iso_PUsub_allpu.Fill(isoowenpu)
            h.iso_noPUsub_allpu.Fill(B/A)

            fill_hist(h.detadphiallpu ,tcThisJet[['dEta','dPhi']].values,tcThisJet.pT.values)
            fill_hist(h.dRgallpu      ,tcThisJet[['dR','energy']].values)
            fill_hist(h.dRptallpu     ,tcThisJet[['dR','pT']].values)
            fill_hist(h.tcptpu        ,tcThisJet.pT.values)
            fill_hist(h.tcenergypu    ,tcThisJet.energy.values)
            fill_hist(h.DRepu         ,tcThisJet.dR.values,tcThisJet.pT.values)
            fill_hist(h.DRppu         ,tcThisJet.dR.values,tcThisJet.energy.values)
            
            if abs(_jets[j].eta)>1.5 and abs(_jets[j].eta)<1.9:
                h.iso_PUsub_eta15pu.Fill(isoowenpu)
                h.iso_noPUsub_eta15pu.Fill(B/A)
                fill_hist(h.detadphieta15pu      ,tcThisJet[['dEta','dPhi']].values, tcThisJet.pT.values)
                fill_hist(h.dReta15pu            ,tcThisJet[['dR','pT']].values)
                
            elif abs(_jets[j].eta)>1.9 and abs(_jets[j].eta)<2.3:
                h.iso_PUsub_eta20pu.Fill(isoowenpu)
                h.iso_noPUsub_eta20pu.Fill(B/A)
                fill_hist(h.detadphieta20pu      ,tcThisJet[['dEta','dPhi']].values, tcThisJet.pT.values)
                fill_hist(h.dReta20pu            ,tcThisJet[['dR','pT']].values)

            elif abs(_jets[j].eta)>2.3 and abs(_jets[j].eta)<2.7:
                h.iso_PUsub_eta25pu.Fill(isoowenpu)
                h.iso_noPUsub_eta25pu.Fill(B/A)
                fill_hist(h.detadphieta25pu      ,tcThisJet[['dEta','dPhi']].values, tcThisJet.pT.values)
                fill_hist(h.dReta25pu            ,tcThisJet[['dR','pT']].values)

            elif abs(_jets[j].eta)>2.7:
                h.iso_PUsub_eta30pu.Fill(isoowenpu)
                h.iso_noPUsub_eta30pu.Fill(B/A)
                fill_hist(h.detadphieta30pu      ,tcThisJet[['dEta','dPhi']].values, tcThisJet.pT.values)
                fill_hist(h.dReta30pu            ,tcThisJet[['dR','pT']].values)

        if args.time:
            print "    process UnMatched in", time.clock()-start

	nmatchedJets=nmatchedJets + len(matchedJets)
	nunmatchedJets=nunmatchedJets + len(unmatchedJets)

        h.jetCount.Fill(0,len(genjetVector))
        h.jetCount.Fill(1,len(matchedJets))
        h.jetCount.Fill(2,len(unmatchedJets))

        del tc
        del genjetVector
        gc.collect()

        
    
    del fulldf
    del fulldfGen
    gc.collect()

    sys.stdout.flush()


print  '===nGenParticlesinHGCAL', nGenParticlesinHGCAL
print  '===matchedjets', nmatchedJets
print  '===unmatchedjets', nunmatchedJets

nmatchedJets = max(nmatchedJets,1)
nunmatchedJets = max(nunmatchedJets,1)

if not skipFill:
    h.dRgall.Scale(1.0/nmatchedJets)
    h.dRptall.Scale(1.0/nmatchedJets)
    h.iso_PUsub_all.Scale(1.0/nmatchedJets)
    h.iso_noPUsub_all.Scale(1.0/nmatchedJets)
    h.dReta15.Scale(1.0/nmatchedJets)
    h.iso_PUsub_eta15.Scale(1.0/nmatchedJets)
    h.iso_noPUsub_eta15.Scale(1.0/nmatchedJets)
    h.dReta20.Scale(1.0/nmatchedJets)
    h.iso_PUsub_eta20.Scale(1.0/nmatchedJets)
    h.iso_noPUsub_eta20.Scale(1.0/nmatchedJets)
    h.dReta25.Scale(1.0/nmatchedJets)
    h.iso_PUsub_eta25.Scale(1.0/nmatchedJets)
    h.iso_noPUsub_eta25.Scale(1.0/nmatchedJets)
    h.dReta30.Scale(1.0/nmatchedJets)
    h.iso_PUsub_eta30.Scale(1.0/nmatchedJets)
    h.iso_noPUsub_eta30.Scale(1.0/nmatchedJets)
    
    h.dRgallpu.Scale(1.0/nunmatchedJets)
    h.dRptallpu.Scale(1.0/nunmatchedJets)
    h.iso_PUsub_allpu.Scale(1.0/nunmatchedJets)
    h.iso_noPUsub_allpu.Scale(1.0/nunmatchedJets)
    h.dReta15pu.Scale(1.0/nunmatchedJets)
    h.iso_PUsub_eta15pu.Scale(1.0/nunmatchedJets)
    h.iso_noPUsub_eta15pu.Scale(1.0/nunmatchedJets)
    h.dReta20pu.Scale(1.0/nunmatchedJets)
    h.iso_PUsub_eta20pu.Scale(1.0/nunmatchedJets)
    h.iso_noPUsub_eta20pu.Scale(1.0/nunmatchedJets)
    h.dReta25pu.Scale(1.0/nunmatchedJets)
    h.iso_PUsub_eta25pu.Scale(1.0/nunmatchedJets)
    h.iso_noPUsub_eta25pu.Scale(1.0/nunmatchedJets)
    h.dReta30pu.Scale(1.0/nunmatchedJets)
    h.iso_PUsub_eta30pu.Scale(1.0/nunmatchedJets)
    h.iso_noPUsub_eta30pu.Scale(1.0/nunmatchedJets)
    
    h.hcalhisto.Write()
