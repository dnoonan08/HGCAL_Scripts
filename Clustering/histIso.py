from ROOT import TFile, TH1F, TH2F, TProfile


class MakeHist():
    def __init__(self,_hName):
        self.hcalhisto = TFile (_hName,'RECREATE')

        self.genJetPt = TH1F("genJetPt","pT of gen parton", 100,0,200)
        self.genJetPtMatched = TH1F("genJetPtMatched","pT of matched gen parton", 100,0,200)
        self.genVsClusterPt = TH2F("genVsClusterPt","pT of gen parton vs pT of clustered jet", 100,0,200, 100,0,200)
        self.genVsClusterPt_puSub = TH2F("genVsClusterPt_puSub","pT of gen parton vs PU subtracted pT of clustered jet", 100,0,200, 100,0,200)

        self.jetPt = TH1F("jetPt","pt of clustered jet",100,0,200)
        self.jetPt_puSub = TH1F("jetPt_puSub","pt of clustered jet with PU subtraction",100,0,200)

        self.DRp = TH1F ('DRp','dR vs. tc_pt',100,0.0,0.5)
        self.DRe = TH1F ('DRe','dR vs. tc_energy',100,0.0,0.5)
        self.dRptall = TH2F ('dRptall','dR vs. tc_pt',100,0.0,0.5,50,0,100)
        self.dRgall = TH2F ('dRgall','dR vs. tc_energy',100,0.0,0.5,50,0,50)
        self.tcenergy = TH1F ('tcenergy','tc_energy',50,0,100)
        self.tcpt = TH1F ('tcpt','tc_pt',50,0,100)
        self.isoowenplotall = TH1F ('isoowenplotall','isoRatio owen R01',500,0,10)
        self.detadphiall = TH2F ('detadphiall','deta vs. dphi',500,-1.0,1.0,500,-1.0,1.0)

        self.dReta15 = TH2F ('dRgeta15','dRg',100,0.0,0.5,50,0,100)
        self.isoowenploteta15 = TH1F ('isoowenploteta15','isoRatio owen R01',500,0,10)
        self.detadphieta15 = TH2F ('detadphieta15','deta vs. dphi',500,-1.0,1.0,500,-1.0,1.0)

        self.dReta20 = TH2F ('dRgeta20','dRg',100,0.0,0.5,50,0,100)
        self.isoowenploteta20 = TH1F ('isoowenploteta20','isoRatio owen R01',500,0,10)
        self.detadphieta20 = TH2F ('detadphieta20','deta vs. dphi',500,-1.0,1.0,500,-1.0,1.0)

        self.dReta25 = TH2F ('dRgeta25','dRg',100,0.0,0.5,50,0,100)
        self.isoowenploteta25 = TH1F ('isoowenploteta25','isoRatio owen R01',500,0,10)
        self.detadphieta25 = TH2F ('detadphieta25','deta vs. dphi',500,-1.0,1.0,500,-1.0,1.0)

        self.dReta30 = TH2F ('dRgeta30','dRg',100,0.0,0.5,50,0,100)
        self.isoowenploteta30 = TH1F ('isoowenploteta30','isoRatio owen R01',500,0,10)
        self.detadphieta30 = TH2F ('detadphieta30','deta vs. dphi',500,-1.0,1.0,500,-1.0,1.0)

        self.DRppu = TH1F ('DRppu','dR vs. tc_pt PU',100,0.0,0.5)
        self.DRepu = TH1F ('DRepu','dR vs. tc_energy PU',100,0.0,0.5)
        self.dRptallpu = TH2F ('dRptallpu','dR vs. tc_pt PU',100,0.0,0.5,50,0,100)
        self.dRgallpu = TH2F ('dRgallpu','dR vs. tc_energy PU',100,0.0,0.5,50,0,100)
        self.tcenergypu = TH1F ('tcenergypu','tc_energy PU',50,0,100)
        self.tcptpu = TH1F ('tcptpu','tc_pt PU',50,0,100)
        self.isoowenplotallpu = TH1F ('isoowenplotallpu','isoRatio owen R01',500,0,10)
        self.detadphiallpu = TH2F ('detadphiallpu','deta vs. dphi',500,-1.0,1.0,500,-1.0,1.0)

        self.dReta15pu = TH2F ('dRgeta15pu','dRg',100,0.0,0.5,50,0,100)
        self.isoowenploteta15pu = TH1F ('isoowenploteta15pu','isoRatio owen R01',500,0,10)
        self.detadphieta15pu = TH2F ('detadphieta15pu','deta vs. dphi',500,-1.0,1.0,500,-1.0,1.0)

        self.dReta20pu = TH2F ('dRgeta20pu','dRg',100,0.0,0.5,50,0,100)
        self.isoowenploteta20pu = TH1F ('isoowenploteta20pu','isoRatio owen R01',500,0,10)
        self.detadphieta20pu = TH2F ('detadphieta20pu','deta vs. dphi',500,-1.0,1.0,500,-1.0,1.0)

        self.dReta25pu = TH2F ('dRgeta25pu','dRg',100,0.0,0.5,50,0,100)
        self.isoowenploteta25pu = TH1F ('isoowenploteta25pu','isoRatio owen R01',500,0,10)
        self.detadphieta25pu = TH2F ('detadphieta25pu','deta vs. dphi',500,-1.0,1.0,500,-1.0,1.0)

        self.dReta30pu = TH2F ('dRgeta30pu','dRg',100,0.0,0.5,50,0,100)
        self.isoowenploteta30pu = TH1F ('isoowenploteta30pu','isoRatio owen R01',500,0,10)
        self.detadphieta30pu = TH2F ('detadphieta30pu','deta vs. dphi',500,-1.0,1.0,500,-1.0,1.0)
        
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
    
