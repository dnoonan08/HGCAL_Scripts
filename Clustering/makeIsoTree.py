from ROOT import TFile, TTree
from array import array
#import numpy as np

class IsoOutputTree():

    def __init__(self, _fName):
        self.outFile = TFile(_fName,'RECREATE')


        maxN = 10
        self.nGen = array('i',[0])
        self.genPt  = array('f',maxN*[0.])
        self.genEta = array('f',maxN*[0.])
        self.genPhi = array('f',maxN*[0.])
        self.genEn  = array('f',maxN*[0.])
        self.genPID = array('i',maxN*[0])
        self.nGenJet = array('i',[0])
        self.genJetPt  = array('f',maxN*[0.])
        self.genJetEta = array('f',maxN*[0.])
        self.genJetPhi = array('f',maxN*[0.])
        self.genJetEn  = array('f',maxN*[0.])        
        self.nJet = array('i',[0])
        maxN = 75
        self.jetPt    = array('f',maxN*[0.])
        self.jetEta   = array('f',maxN*[0.])
        self.jetPhi   = array('f',maxN*[0.])
        self.jetPtR01 = array('f',maxN*[0.])
        self.jetPtR02 = array('f',maxN*[0.])
        self.jetPtR04 = array('f',maxN*[0.])
        self.jetEnR01 = array('f',maxN*[0.])
        self.jetEnR02 = array('f',maxN*[0.])
        self.jetEnR04 = array('f',maxN*[0.])
        self.jetIso   = array('f',maxN*[0.])
        self.jetMinGenDR   = array('f',maxN*[0.])
        self.jetMinGenJetDR   = array('f',maxN*[0.])
        self.jetGenMatch    = array('i',maxN*[0])
        self.jetGenJetMatch = array('i',maxN*[0])
        
        self.tree = TTree("jetTree","jetTree")
        
        self.tree.Branch('nGen',self.nGen,'nGen/I')
        self.tree.Branch('genPt',self.genPt,'genPt[nGen]/F')
        self.tree.Branch('genEta',self.genEta,'genEta[nGen]/F')
        self.tree.Branch('genPhi',self.genPhi,'genPhi[nGen]/F')
        self.tree.Branch('genEn',self.genEn,'genEn[nGen]/F')
        self.tree.Branch('genPID',self.genPID,'genPID[nGen]/I')


        self.tree.Branch('nGenJet',self.nGenJet,'nGenJet/I')
        self.tree.Branch('genJetPt',self.genJetPt,'genJetPt[nGenJet]/F')
        self.tree.Branch('genJetEta',self.genJetEta,'genJetEta[nGenJet]/F')
        self.tree.Branch('genJetPhi',self.genJetPhi,'genJetPhi[nGenJet]/F')
        self.tree.Branch('genJetEn',self.genJetEn,'genJetEn[nGenJet]/F')

        self.tree.Branch('nJet',self.nJet,'nJet/I')
        self.tree.Branch('jetPt',self.jetPt,'jetPt[nJet]/F')
        self.tree.Branch('jetEta',self.jetEta,'jetEta[nJet]/F')
        self.tree.Branch('jetPhi',self.jetPhi,'jetPhi[nJet]/F')
        self.tree.Branch('jetPtR01',self.jetPtR01,'jetPtR01[nJet]/F')
        self.tree.Branch('jetPtR02',self.jetPtR02,'jetPtR02[nJet]/F')
        self.tree.Branch('jetPtR04',self.jetPtR04,'jetPtR04[nJet]/F')
        self.tree.Branch('jetEnR01',self.jetEnR01,'jetEnR01[nJet]/F')
        self.tree.Branch('jetEnR02',self.jetEnR02,'jetEnR02[nJet]/F')
        self.tree.Branch('jetEnR04',self.jetEnR04,'jetEnR04[nJet]/F')
        self.tree.Branch('jetIso',self.jetIso,'jetIso[nJet]/F')
        self.tree.Branch('jetMinGenDR',self.jetMinGenDR,'jetMinGenDR[nJet]/F')
        self.tree.Branch('jetMinGenJetDR',self.jetMinGenJetDR,'jetMinGenJetDR[nJet]/F')
        self.tree.Branch('jetGenMatch',self.jetGenMatch,'jetGenMatch[nJet]/I')
        self.tree.Branch('jetGenJetMatch',self.jetGenJetMatch,'jetGenJetMatch[nJet]/I')
