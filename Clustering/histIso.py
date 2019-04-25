from ROOT import TFile, TH1F, TH2F, TProfile

class MakeHist():
    def __init__(self,_hName):
        self.hcalhisto = TFile (_hName,'RECREATE')

        self.genJetPt = TH1F("genJetPt","pT of gen parton", 100,0,200)
        self.genJetPtMatched = TH1F("genJetPtMatched","pT of matched gen parton", 100,0,200)
        self.genVsClusterPt = TH2F("genVsClusterPt","pT of gen parton vs pT of clustered jet", 100,0,200, 100,0,200)
        self.genVsClusterPt_puSub = TH2F("genVsClusterPt_puSub","pT of gen parton vs PU subtracted pT of clustered jet", 100,0,200, 100,0,200)

        self.jetCount = TH1F("jetCount","count gen, matched, and unmatched jets", 3,0,3)

        self.genPtVsClusterIso = TH2F("genPtVsClusterIso","pT of gen parton vs iso clustered jet", 100,0,200, 100,0,10)
        self.clusterPtVsClusterIso = TH2F("clusterPtVsClusterIso","pT of clustered jet vs iso clustered jet", 100,0,200, 100,0,10)

        self.jetPt = TH1F("jetPt","pt of clustered jet",100,0,200)
        self.jetPt_puSub = TH1F("jetPt_puSub","pt of clustered jet with PU subtraction",100,0,200)

        self.DRp = TH1F ('DRp','dR vs. tc_pt',25,0.0,0.5)
        self.DRe = TH1F ('DRe','dR vs. tc_energy',25,0.0,0.5)
        self.dRptall = TH2F ('dRptall','dR vs. tc_pt',25,0.0,0.5,50,0,100)
        self.dRgall = TH2F ('dRgall','dR vs. tc_energy',25,0.0,0.5,50,0,100)
        self.tcenergy = TH1F ('tcenergy','tc_energy',50,0,100)
        self.tcpt = TH1F ('tcpt','tc_pt',50,0,100)
        self.iso_PUsub_all = TH1F ('iso_PUsub_all','isoRatio R01',100,0,10)
        self.iso_noPUsub_all = TH1F ('iso_noPUsub_all','isoRatio R01',100,0,10)
        self.detadphiall = TH2F ('detadphiall','deta vs. dphi',50,-.5,.5,50,-.5,.5)

        self.dReta15 = TH2F ('dRgeta15','dRg',25,0.0,0.5,50,0,100)
        self.iso_noPUsub_eta15 = TH1F ('iso_noPUsub_eta15','isoRatio owen R01',100,0,10)
        self.iso_PUsub_eta15 = TH1F ('iso_PUsub_eta15','isoRatio owen R01',100,0,10)
        self.detadphieta15 = TH2F ('detadphieta15','deta vs. dphi',50,-.5,.5,50,-.5,.5)

        self.dReta20 = TH2F ('dRgeta20','dRg',25,0.0,0.5,50,0,100)
        self.iso_PUsub_eta20 = TH1F ('iso_PUsub_eta20','isoRatio owen R01',100,0,10)
        self.iso_noPUsub_eta20 = TH1F ('iso_noPUsub_eta20','isoRatio owen R01',100,0,10)
        self.detadphieta20 = TH2F ('detadphieta20','deta vs. dphi',50,-.5,.5,50,-.5,.5)

        self.dReta25 = TH2F ('dRgeta25','dRg',25,0.0,0.5,50,0,100)
        self.iso_PUsub_eta25 = TH1F ('iso_PUsub_eta25','isoRatio owen R01',100,0,10)
        self.iso_noPUsub_eta25 = TH1F ('iso_noPUsub_eta25','isoRatio owen R01',100,0,10)
        self.detadphieta25 = TH2F ('detadphieta25','deta vs. dphi',50,-.5,.5,50,-.5,.5)

        self.dReta30 = TH2F ('dRgeta30','dRg',25,0.0,0.5,50,0,100)
        self.iso_PUsub_eta30 = TH1F ('iso_PUsub_eta30','isoRatio owen R01',100,0,10)
        self.iso_noPUsub_eta30 = TH1F ('iso_noPUsub_eta30','isoRatio owen R01',100,0,10)
        self.detadphieta30 = TH2F ('detadphieta30','deta vs. dphi',50,-.5,.5,50,-.5,.5)

        self.DRppu = TH1F ('DRppu','dR vs. tc_pt PU',25,0.0,0.5)
        self.DRepu = TH1F ('DRepu','dR vs. tc_energy PU',25,0.0,0.5)
        self.dRptallpu = TH2F ('dRptallpu','dR vs. tc_pt PU',25,0.0,0.5,50,0,100)
        self.dRgallpu = TH2F ('dRgallpu','dR vs. tc_energy PU',25,0.0,0.5,50,0,100)
        self.tcenergypu = TH1F ('tcenergypu','tc_energy PU',50,0,100)
        self.tcptpu = TH1F ('tcptpu','tc_pt PU',50,0,100)
        self.iso_PUsub_allpu = TH1F ('iso_PUsub_allpu','isoRatio owen R01',100,0,10)
        self.iso_noPUsub_allpu = TH1F ('iso_noPUsub_allpu','isoRatio owen R01',100,0,10)
        self.detadphiallpu = TH2F ('detadphiallpu','deta vs. dphi',50,-.5,.5,50,-.5,.5)

        self.dReta15pu = TH2F ('dRgeta15pu','dRg',25,0.0,0.5,50,0,100)
        self.iso_PUsub_eta15pu = TH1F ('iso_PUsub_eta15pu','isoRatio owen R01',100,0,10)
        self.iso_noPUsub_eta15pu = TH1F ('iso_noPUsub_eta15pu','isoRatio owen R01',100,0,10)
        self.detadphieta15pu = TH2F ('detadphieta15pu','deta vs. dphi',50,-.5,.5,50,-.5,.5)

        self.dReta20pu = TH2F ('dRgeta20pu','dRg',25,0.0,0.5,50,0,100)
        self.iso_PUsub_eta20pu = TH1F ('iso_PUsub_eta20pu','isoRatio owen R01',100,0,10)
        self.iso_noPUsub_eta20pu = TH1F ('iso_noPUsub_eta20pu','isoRatio owen R01',100,0,10)
        self.detadphieta20pu = TH2F ('detadphieta20pu','deta vs. dphi',50,-.5,.5,50,-.5,.5)

        self.dReta25pu = TH2F ('dRgeta25pu','dRg',25,0.0,0.5,50,0,100)
        self.iso_PUsub_eta25pu = TH1F ('iso_PUsub_eta25pu','isoRatio owen R01',100,0,10)
        self.iso_noPUsub_eta25pu = TH1F ('iso_noPUsub_eta25pu','isoRatio owen R01',100,0,10)
        self.detadphieta25pu = TH2F ('detadphieta25pu','deta vs. dphi',50,-.5,.5,50,-.5,.5)

        self.dReta30pu = TH2F ('dRgeta30pu','dRg',25,0.0,0.5,50,0,100)
        self.iso_PUsub_eta30pu = TH1F ('iso_PUsub_eta30pu','isoRatio owen R01',100,0,10)
        self.iso_noPUsub_eta30pu = TH1F ('iso_noPUsub_eta30pu','isoRatio owen R01',100,0,10)
        self.detadphieta30pu = TH2F ('detadphieta30pu','deta vs. dphi',50,-.5,.5,50,-.5,.5)
        
        NBins=100
        binMin=0.0
        binMax=0.5
        self.hArea = TProfile('hArea','hArea',NBins,binMin,binMax)
        for j in range(NBins):
            binwidth=(binMax-binMin)/NBins
            lowbin=(j)*binwidth
            highbin=lowbin+binwidth
            area=3.14159265*(highbin**2-lowbin**2)
            self.hArea.Fill(j*binwidth,area)
    
