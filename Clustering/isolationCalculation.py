
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
parser.add_argument('--NFiles', "--nFiles", default=-1, type=int)
parser.add_argument('--name', default="")
parser.add_argument('--verbose','-v',action='store_true',help='Print verbose output of matching')
parser.add_argument('--genjet',action='store_true',help='Use genJet instead of gen partons')
args = parser.parse_args()

print args

if not (args.vbf or args.particlegun):
    print "Specify vbf or particle gun"
    exit()

import uproot
import numpy as np
import pandas as pd
from root_numpy import fill_hist

from subprocess import Popen,PIPE

from pyjet import cluster,DTYPE_EP,DTYPE_PTEPM

#import histvbf as h
from histvbf import *

print "Starting"



nGenVBFinHGCAL = 0

nmatchedJets = 0
nunmatchedJets = 0

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

h = MakeHist(outputFileName)

fileList = []

if not args.name=="":
    outputFileName = outputFileName.replace(".root","_%s.root"%args.name)

cmd = "xrdfs root://cmseos.fnal.gov ls %s"%eosDir
dirContents,stderr = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE).communicate()
dirContentsList = dirContents.split("\n")
for fName in dirContentsList:
    if fileNameContent in fName:
        fileList.append("root://cmseos.fnal.gov/%s"%(fName))
        
if args.NFiles==-1:
    nFiles = len(fileList)
else:
    nFiles = args.NFiles

for fName in fileList[:nFiles]:

    _tree = uproot.open(fName,xrootdsource=dict(chunkbytes=1024**3, limitbytes=1024**3))["hgcalTriggerNtuplizer/HGCalTriggerNtuple"]
    if args.N==-1:
        N = _tree.numentries
    else:
        N = args.N

    print "File %s"%fName
    fulldf = _tree.pandas.df(["tc_pt","tc_energy","tc_eta","tc_mipPt","tc_phi"],entrystart=0,entrystop=N)
    if not args.genjet:
        fulldfGen = _tree.pandas.df(["gen_pt","gen_eta","gen_phi","gen_status","gen_pdgid","gen_energy"],entrystart=0,entrystop=N)
    else:
        fulldfGen = _tree.pandas.df(["genjet_pt","genjet_eta","genjet_phi","genjet_energy"],entrystart=0,entrystop=N)

    print "Loaded tree ", time.clock()-start
    print N

    for i_event in range(N):
        start = time.clock()
        print '-'*20
    	print "     %i"%i_event   
        # load into dictionary for creating pandas df
        tc = fulldf.loc[i_event,['tc_pt','tc_eta','tc_phi','tc_energy']]
        tc.columns=["pT","eta","phi","energy"]

    #go dataframe to np array, formatted for input into fastjet    
        tcVectors = np.array(tc.to_records(index=False).astype([(u'pT', '<f8'), (u'eta', '<f8'), (u'phi', '<f8'), (u'energy', '<f8')]) ) 
        print "    Loaded in", time.clock()-start
        clusterVals = cluster(tcVectors,R=0.4,algo="antikt")
        _jets = clusterVals.inclusive_jets(ptmin=20)
        print "    Clustered jets in", time.clock()-start
        print "      --- len =",len(tcVectors), (time.clock()-start)/len(tcVectors)

        if not args.genjet:
            genjetDF = fulldfGen.loc[i_event,['gen_pt','gen_eta','gen_phi','gen_status','gen_pdgid','gen_energy']]
            if args.vbf:
                genjetDF = genjetDF[(genjetDF.gen_status==23) & (abs(genjetDF.gen_pdgid)<6) & (abs(genjetDF.gen_eta)>1.5) & (abs(genjetDF.gen_eta)<3.)]
            else:
                genjetDF = genjetDF[(genjetDF.gen_status==23) & (abs(genjetDF.gen_pdgid)==args.pid) & (abs(genjetDF.gen_eta)>1.5) & (abs(genjetDF.gen_eta)<3.)]
        else:
            genjetDF = fulldfGen.loc[i_event,['genjet_pt','genjet_eta','genjet_phi','genjet_energy']]
            genjetDF.columns = ['gen_pt','gen_eta','gen_phi','gen_energy']
            

        fill_hist(h.genJetPt   ,genjetDF.gen_pt.values)

        genjetVector = genjetDF.values
        
        nGenVBFinHGCAL += len(genjetVector)

        print "    Got Gen in", time.clock()-start

        matchedJets = []
        genMatch = []


        for i in range(len(genjetVector)):
            if args.verbose:
                foundMatch = False
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
            if not foundMatch:
                print "No Match Found"


        print "    GenMatch in", time.clock()-start
        unmatchedJets = [x for x in range(len(_jets)) if x not in matchedJets]
        unmatchedGen = [x for x in range(len(genjetVector)) if x not in genMatch]

        for i in range(len(matchedJets)):
            j = matchedJets[i]

            tc['dEta'] = tc['eta']-_jets[j].eta
            tc['dPhi'] = tc['phi']-_jets[j].phi
            tc['dR'] = (tc.dEta**2 + tc.dPhi**2)**0.5
            
            cut = tc.dR<0.45
            tcThisJet = tc[cut] #.loc[dR2<0.45,:]
            # dEtaThisJet = dEta[dR<0.45].values
            # dPhiThisJet = dPhi[dR<0.45].values
            # dRThisJet = dR[dR<0.45].values
            # ptThisJet = tcThisJet.pT.values
	    # energyThisJet = tcThisJet.energy.values
            
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

            h.isoowenplotall.Fill(isoowen)

            fill_hist(h.detadphiall   ,tcThisJet[['dEta','dPhi']].values,tcThisJet.pT.values)
            fill_hist(h.dRgall        ,tcThisJet[['dR','energy']].values)
            fill_hist(h.dRptall       ,tcThisJet[['dR','pT']].values)
            fill_hist(h.tcpt          ,tcThisJet.pT.values)
            fill_hist(h.tcenergy      ,tcThisJet.energy.values)
            fill_hist(h.DRe           ,tcThisJet.dR.values,tcThisJet.pT.values)
            fill_hist(h.DRp           ,tcThisJet.dR.values,tcThisJet.energy.values)
            
            if abs(_jets[j].eta)>1.5 and abs(_jets[j].eta)<1.9:
                    h.isoowenploteta15.Fill(isoowen)
                    fill_hist(h.detadphieta15        ,tcThisJet[['dEta','dPhi']].values, tcThisJet.pT.values)
                    fill_hist(h.dReta15              ,tcThisJet[['dR','pT']].values)

            elif abs(_jets[j].eta)>1.9 and abs(_jets[j].eta)<2.3:
                    h.isoowenploteta20.Fill(isoowen)
                    fill_hist(h.detadphieta20        ,tcThisJet[['dEta','dPhi']].values, tcThisJet.pT.values)
                    fill_hist(h.dReta20              ,tcThisJet[['dR','pT']].values)

            elif abs(_jets[j].eta)>2.3 and abs(_jets[j].eta)<2.7:
                    h.isoowenploteta25.Fill(isoowen)
                    fill_hist(h.detadphieta25        ,tcThisJet[['dEta','dPhi']].values, tcThisJet.pT.values)
                    fill_hist(h.dReta25              ,tcThisJet[['dR','pT']].values)

            elif abs(_jets[j].eta)>2.7:
                    h.isoowenploteta30.Fill(isoowen)
                    fill_hist(h.detadphieta30        ,tcThisJet[['dEta','dPhi']].values, tcThisJet.pT.values)
                    fill_hist(h.dReta30              ,tcThisJet[['dR','pT']].values)


            # for i_tc in range(len(dEtaThisJet)):
            #     h.detadphiall.Fill(dEtaThisJet[i_tc],dPhiThisJet[i_tc],ptThisJet[i_tc])
            #     h.dRptall.Fill(dRThisJet[i_tc],ptThisJet[i_tc], ptThisJet[i_tc])
            #     h.dRgall.Fill(dRThisJet[i_tc],energyThisJet[i_tc],energyThisJet[i_tc])
	    #     h.tcpt.Fill(ptThisJet[i_tc],ptThisJet[i_tc])
	    #     h.tcenergy.Fill(energyThisJet[i_tc],energyThisJet[i_tc])
	    #     h.DRp.Fill(dRThisJet[i_tc],ptThisJet[i_tc])
	    #     h.DRe.Fill(dRThisJet[i_tc],energyThisJet[i_tc])
	    # if abs(genjetVector[0,1])>1.5 and abs(genjetVector[0,1])<1.9:
            #         h.isoowenploteta15.Fill(isoowen)
            #         for i in range(len(dEtaThisJet)):
            #            h.detadphieta15.Fill(dEtaThisJet[i_tc],dPhiThisJet[i_tc],ptThisJet[i_tc])
            #            h.dReta15.Fill(dRThisJet[i_tc],ptThisJet[i_tc])
            # elif abs(genjetVector[0,1])>1.9 and abs(genjetVector[0,1])<2.3:
            #         h.isoowenploteta20.Fill(isoowen)
            #         for i in range(len(dEtaThisJet)):
            #            h.detadphieta20.Fill(dEtaThisJet[i_tc],dPhiThisJet[i_tc],ptThisJet[i_tc])
            #            h.dReta20.Fill(dRThisJet[i_tc],ptThisJet[i_tc])
            # elif abs(genjetVector[0,1])>2.3 and abs(genjetVector[0,1])<2.7:
            #         h.isoowenploteta25.Fill(isoowen)
            #         for i in range(len(dEtaThisJet)):
            #            h.detadphieta25.Fill(dEtaThisJet[i_tc],dPhiThisJet[i_tc],ptThisJet[i_tc])
            #            h.dReta25.Fill(dRThisJet[i_tc],ptThisJet[i_tc])
            # elif abs(genjetVector[0,1])>2.7:
            #         h.isoowenploteta30.Fill(isoowen)
            #         for i in range(len(dEtaThisJet)):
            #            h.detadphieta30.Fill(dEtaThisJet[i_tc],dPhiThisJet[i_tc],ptThisJet[i_tc])
            #            h.dReta30.Fill(dRThisJet[i_tc],ptThisJet[i_tc])

        print "    process Matched in", time.clock()-start
                       
        for j in unmatchedJets:
            tc['dEta'] = tc['eta']-_jets[j].eta
            tc['dPhi'] = tc['phi']-_jets[j].phi
            tc['dR'] = (tc.dEta**2 + tc.dPhi**2)**0.5

            cut = tc.dR<0.45
            tcThisJet = tc[cut] #.loc[dR2<0.45,:]
#            dEtaThisJet2 = dEta2[cut].values
#            dPhiThisJet2 = dPhi2[cut].values
#            dRThisJet2 = dR2[cut].values
#            ptThisJet2 = tcThisJet2.pT.values
#	    energyThisJet2 = tcThisJet2.energy.values

            A = tcThisJet.loc[tcThisJet.dR<0.1,'pT'].sum()
            B = tcThisJet.loc[(tcThisJet.dR>0.1) & (tcThisJet.dR<0.2),'pT'].sum()
            C = tcThisJet.loc[(tcThisJet.dR>0.2) & (tcThisJet.dR<0.4),'pT'].sum()

            isoowenpu = (B-(3./12)*C)/(A-(1./12)*C)


            if skipFill: continue

            h.isoowenplotallpu.Fill(isoowenpu)

            fill_hist(h.detadphiallpu ,tcThisJet[['dEta','dPhi']].values,tcThisJet.pT.values)
            fill_hist(h.dRgallpu      ,tcThisJet[['dR','energy']].values)
            fill_hist(h.dRptallpu     ,tcThisJet[['dR','pT']].values)
            fill_hist(h.tcptpu        ,tcThisJet.pT.values)
            fill_hist(h.tcenergypu    ,tcThisJet.energy.values)
            fill_hist(h.DRepu         ,tcThisJet.dR.values,tcThisJet.pT.values)
            fill_hist(h.DRppu         ,tcThisJet.dR.values,tcThisJet.energy.values)
            
            if abs(_jets[j].eta)>1.5 and abs(_jets[j].eta)<1.9:
                    h.isoowenploteta15pu.Fill(isoowenpu)
                    fill_hist(h.detadphieta15pu      ,tcThisJet[['dEta','dPhi']].values, tcThisJet.pT.values)
                    fill_hist(h.dReta15pu            ,tcThisJet[['dR','pT']].values)

            elif abs(_jets[j].eta)>1.9 and abs(_jets[j].eta)<2.3:
                    h.isoowenploteta20pu.Fill(isoowenpu)
                    fill_hist(h.detadphieta20pu      ,tcThisJet[['dEta','dPhi']].values, tcThisJet.pT.values)
                    fill_hist(h.dReta20pu            ,tcThisJet[['dR','pT']].values)

            elif abs(_jets[j].eta)>2.3 and abs(_jets[j].eta)<2.7:
                    h.isoowenploteta25pu.Fill(isoowenpu)
                    fill_hist(h.detadphieta25pu      ,tcThisJet[['dEta','dPhi']].values, tcThisJet.pT.values)
                    fill_hist(h.dReta25pu            ,tcThisJet[['dR','pT']].values)

            elif abs(_jets[j].eta)>2.7:
                    h.isoowenploteta30pu.Fill(isoowenpu)
                    fill_hist(h.detadphieta30pu      ,tcThisJet[['dEta','dPhi']].values, tcThisJet.pT.values)
                    fill_hist(h.dReta30pu            ,tcThisJet[['dR','pT']].values)


        print "    process UnMatched in", time.clock()-start

	nmatchedJets=nmatchedJets + len(matchedJets)
	nunmatchedJets=nunmatchedJets + len(unmatchedJets)

print  '===nGenVBFinHGCAL', nGenVBFinHGCAL
print  '===matchedjets', nmatchedJets
print  '===unmatchedjets', nunmatchedJets

if not skipFill:
    h.dRgall.Scale(1.0/nmatchedJets)
    h.dRptall.Scale(1.0/nmatchedJets)
    h.isoowenplotall.Scale(1.0/nmatchedJets)
    h.dReta15.Scale(1.0/nmatchedJets)
    h.isoowenploteta15.Scale(1.0/nmatchedJets)
    h.dReta20.Scale(1.0/nmatchedJets)
    h.isoowenploteta20.Scale(1.0/nmatchedJets)
    h.dReta25.Scale(1.0/nmatchedJets)
    h.isoowenploteta25.Scale(1.0/nmatchedJets)
    h.dReta30.Scale(1.0/nmatchedJets)
    h.isoowenploteta30.Scale(1.0/nmatchedJets)
    
    h.dRgallpu.Scale(1.0/nunmatchedJets)
    h.dRptallpu.Scale(1.0/nunmatchedJets)
    h.isoowenplotallpu.Scale(1.0/nunmatchedJets)
    h.dReta15pu.Scale(1.0/nunmatchedJets)
    h.isoowenploteta15pu.Scale(1.0/nunmatchedJets)
    h.dReta20pu.Scale(1.0/nunmatchedJets)
    h.isoowenploteta20pu.Scale(1.0/nunmatchedJets)
    h.dReta25pu.Scale(1.0/nunmatchedJets)
    h.isoowenploteta25pu.Scale(1.0/nunmatchedJets)
    h.dReta30pu.Scale(1.0/nunmatchedJets)
    h.isoowenploteta30pu.Scale(1.0/nunmatchedJets)
    
    h.hcalhistoVBF.Write()
