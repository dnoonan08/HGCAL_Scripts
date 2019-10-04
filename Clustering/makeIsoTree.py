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
        self.genMinDR  = array('f',maxN*[0.])
        self.genMinDPt  = array('f',maxN*[0.])
        self.genMinDEn  = array('f',maxN*[0.])

        self.nGen_Full = array('i',[0])
        self.genPt_Full  = array('f',maxN*[0.])
        self.genEta_Full = array('f',maxN*[0.])
        self.genPhi_Full = array('f',maxN*[0.])
        self.genEn_Full  = array('f',maxN*[0.])
        self.genPID_Full = array('i',maxN*[0])
        self.genMinDR_Full  = array('f',maxN*[0.])
        self.genMinDPt_Full  = array('f',maxN*[0.])
        self.genMinDEn_Full  = array('f',maxN*[0.])

        self.nGenJet = array('i',[0])
        self.genJetPt  = array('f',maxN*[0.])
        self.genJetEta = array('f',maxN*[0.])
        self.genJetPhi = array('f',maxN*[0.])
        self.genJetEn  = array('f',maxN*[0.])        


        maxN = 20
        self.nGenJet_Full = array('i',[0])
        self.genJetPt_Full  = array('f',maxN*[0.])
        self.genJetEta_Full = array('f',maxN*[0.])
        self.genJetPhi_Full = array('f',maxN*[0.])
        self.genJetEn_Full  = array('f',maxN*[0.])        
        self.genJetPartonMatch_Full  = array('i',maxN*[0])

        maxN = 75
        self.nJet = array('i',[0])
        self.jetPt    = array('f',maxN*[0.])
        self.jetEta   = array('f',maxN*[0.])
        self.jetPhi   = array('f',maxN*[0.])
        self.jetPtR01 = array('f',maxN*[0.])
        self.jetPtR02 = array('f',maxN*[0.])
        self.jetPtR03 = array('f',maxN*[0.])
        self.jetPtR04 = array('f',maxN*[0.])
        self.jetPtR05 = array('f',maxN*[0.])
        self.jetPtR06 = array('f',maxN*[0.])

        self.jetEnR01 = array('f',maxN*[0.])
        self.jetEnR02 = array('f',maxN*[0.])
        self.jetEnR03 = array('f',maxN*[0.])
        self.jetEnR04 = array('f',maxN*[0.])
        self.jetEnR05 = array('f',maxN*[0.])
        self.jetEnR06 = array('f',maxN*[0.])

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
        self.tree.Branch('genMinDR',self.genMinDR,'genMinDR[nGen]/F')
        self.tree.Branch('genMinDPt',self.genMinDPt,'genMinDPt[nGen]/F')
        self.tree.Branch('genMinDEn',self.genMinDEn,'genMinDEn[nGen]/F')


        self.tree.Branch('nGen_Full',self.nGen_Full,'nGen_Full/I')
        self.tree.Branch('genPt_Full',self.genPt_Full,'genPt_Full[nGen]/F')
        self.tree.Branch('genEta_Full',self.genEta_Full,'genEta_Full[nGen]/F')
        self.tree.Branch('genPhi_Full',self.genPhi_Full,'genPhi_Full[nGen]/F')
        self.tree.Branch('genEn_Full',self.genEn_Full,'genEn_Full[nGen]/F')
        self.tree.Branch('genPID_Full',self.genPID_Full,'genPID_Full[nGen]/I')
        self.tree.Branch('genMinDR_Full',self.genMinDR_Full,'genMinDR_Full[nGen]/F')
        self.tree.Branch('genMinDPt_Full',self.genMinDPt_Full,'genMinDPt_Full[nGen]/F')
        self.tree.Branch('genMinDEn_Full',self.genMinDEn_Full,'genMinDEn_Full[nGen]/F')



        self.tree.Branch('nGenJet',self.nGenJet,'nGenJet/I')
        self.tree.Branch('genJetPt',self.genJetPt,'genJetPt[nGenJet]/F')
        self.tree.Branch('genJetEta',self.genJetEta,'genJetEta[nGenJet]/F')
        self.tree.Branch('genJetPhi',self.genJetPhi,'genJetPhi[nGenJet]/F')
        self.tree.Branch('genJetEn',self.genJetEn,'genJetEn[nGenJet]/F')

        self.tree.Branch('nGenJet_Full',self.nGenJet_Full,'nGenJet_Full/I')
        self.tree.Branch('genJetPt_Full',self.genJetPt_Full,'genJetPt_Full[nGenJet]/F')
        self.tree.Branch('genJetEta_Full',self.genJetEta_Full,'genJetEta_Full[nGenJet]/F')
        self.tree.Branch('genJetPhi_Full',self.genJetPhi_Full,'genJetPhi_Full[nGenJet]/F')
        self.tree.Branch('genJetEn_Full',self.genJetEn_Full,'genJetEn_Full[nGenJet]/F')
        self.tree.Branch('genJetPartonMatch_Full',self.genJetPartonMatch_Full,'genJetPartonMatch_Full[nGenJet]/F')

        self.tree.Branch('nJet',self.nJet,'nJet/I')
        self.tree.Branch('jetPt',self.jetPt,'jetPt[nJet]/F')
        self.tree.Branch('jetEta',self.jetEta,'jetEta[nJet]/F')
        self.tree.Branch('jetPhi',self.jetPhi,'jetPhi[nJet]/F')
        self.tree.Branch('jetPtR01',self.jetPtR01,'jetPtR01[nJet]/F')
        self.tree.Branch('jetPtR02',self.jetPtR02,'jetPtR02[nJet]/F')
        self.tree.Branch('jetPtR03',self.jetPtR03,'jetPtR03[nJet]/F')
        self.tree.Branch('jetPtR04',self.jetPtR04,'jetPtR04[nJet]/F')
        self.tree.Branch('jetPtR05',self.jetPtR05,'jetPtR05[nJet]/F')
        self.tree.Branch('jetPtR06',self.jetPtR06,'jetPtR06[nJet]/F')
        self.tree.Branch('jetEnR01',self.jetEnR01,'jetEnR01[nJet]/F')
        self.tree.Branch('jetEnR02',self.jetEnR02,'jetEnR02[nJet]/F')
        self.tree.Branch('jetEnR03',self.jetEnR03,'jetEnR03[nJet]/F')
        self.tree.Branch('jetEnR04',self.jetEnR04,'jetEnR04[nJet]/F')
        self.tree.Branch('jetEnR05',self.jetEnR05,'jetEnR05[nJet]/F')
        self.tree.Branch('jetEnR06',self.jetEnR06,'jetEnR06[nJet]/F')
        self.tree.Branch('jetIso',self.jetIso,'jetIso[nJet]/F')
        self.tree.Branch('jetMinGenDR',self.jetMinGenDR,'jetMinGenDR[nJet]/F')
        self.tree.Branch('jetMinGenJetDR',self.jetMinGenJetDR,'jetMinGenJetDR[nJet]/F')
        self.tree.Branch('jetGenMatch',self.jetGenMatch,'jetGenMatch[nJet]/I')
        self.tree.Branch('jetGenJetMatch',self.jetGenJetMatch,'jetGenJetMatch[nJet]/I')
