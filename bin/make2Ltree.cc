#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TSystem.h"
#include "TF1.h"
#include "TRandom3.h"

#include <string>
#include <vector>

#include "HeavyIonsAnalysis/topskim/include/tnp_weight.h"
#include "HeavyIonsAnalysis/topskim/include/tnp_electrons.h"
#include "HeavyIonsAnalysis/topskim/include/ForestHiTree.h"
#include "HeavyIonsAnalysis/topskim/include/ForestHLTObject.h"
#include "HeavyIonsAnalysis/topskim/include/ForestLeptons.h"
#include "HeavyIonsAnalysis/topskim/include/ForestSkim.h"
#include "HeavyIonsAnalysis/topskim/include/ForestPFCands.h"
#include "HeavyIonsAnalysis/topskim/include/ForestJets.h"
#include "HeavyIonsAnalysis/topskim/include/JetCorrector.h"
#include "HeavyIonsAnalysis/topskim/include/JetUncertainty.h"
#include "HeavyIonsAnalysis/topskim/include/LumiRun.h"
#include "HeavyIonsAnalysis/topskim/include/HistTool.h"
#include "HeavyIonsAnalysis/topskim/include/PFAnalysis.h"
#include "HeavyIonsAnalysis/topskim/include/LeptonSummary.h"
#include "HeavyIonsAnalysis/topskim/include/ForestGen.h"
#include "HeavyIonsAnalysis/topskim/include/ElectronId.h"
#include "HeavyIonsAnalysis/topskim/include/BtagUncertaintyComputer.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "../scripts/functions.cc"


const bool isDebug = true;
const float lepPtCut  = 15.;
const float muEtaCut = 2.4;
const float eleEtaCut = 2.1;
//see https://indico.cern.ch/event/803679/contributions/3342407/attachments/1808912/2953435/egm-minipog-190308.pdf
const float eeScaleShift = 6.8182E-2/5.9097E-2;
const int firstEEScaleShiftRun = 327402; 
const float barrelEndcapEta[2]={1.4442,1.5660};
const float hem1516Eta[2]={-3.0,-1.9};
const float hem1516Phi[2]={-1.6,-0.9};
const float csvWPList[2]={0.81,0.91};
int csvWP = 0;
// runs with HLT issues: any path using tracking was disabled => includes L3 muon paths. 
// total lumi in these runs is 30.103/ub (total lumi available in golden json unblinde period is 425.349/ub)
std::vector<int> badMuonTriggerRuns={326482,326483,326500,326520,326527,326528,326530,326532,326533,326534,326535,326546,326548,326549,326550,326568,326569,326571};

using namespace std;
using namespace fastjet;
using namespace TMVA;

//dispersion of rapidities of the final state
std::vector<float> getRapidityMoments(std::vector<TLorentzVector> & coll){
  std::vector<float> mom(3,0.);
  for(size_t i=0; i<coll.size(); i++) {
    TLorentzVector pi(coll[i]);
    mom[0]+=pi.Rapidity();
    mom[1]+=pow(pi.Rapidity(),2);
    for(size_t j=i; j<coll.size(); j++) {
      TLorentzVector pj(coll[j]);
      float dy=fabs(pj.Rapidity()-pi.Rapidity());
      mom[2]=max(dy,mom[2]);
    }
  }
  
  mom[0]=mom[0]/float(mom.size());
  mom[1]=sqrt(mom[1]/float(mom.size())-pow(mom[0],2));
  return mom;
}

//Madgraph BR(W->lnu) correction
float getMadgraphBRWlCorrection(int nl) {
  if(nl==0) return 1.0236380625;
  if(nl==1) return 0.987973875;
  if(nl==2) return 0.95355225;
  return 1.0;
}

int getRhoIndex(float eta,std::vector<Double_t> *etaMin=NULL, std::vector<Double_t> *etaMax=NULL){
  if(etaMin==NULL || etaMax==NULL) return -1;
  for(size_t i=0; i<etaMin->size(); i++) {
    if(eta>=etaMin->at(i) && eta<etaMax->at(i)) return i;
  }
  return -1;
}

TRandom3 * smearRand = new TRandom3(42);

float calibratedPt(float pt, float eta, float cen, bool isMC) {

  float newpt;
  float scale = 0.;
  //https://twiki.cern.ch/twiki/pub/CMS/HiEgamma2019/electron-aa-scalesmear-190918.pdf
  if (isMC){
    if (fabs(eta) < 1.45) {
      if      (cen < 10.) scale = 0.971;
      else if (cen < 30.) scale = 0.990;
      else                scale = 1.003;
    }else {
      if      (cen < 10.) scale = 0.910;
      else if (cen < 30.) scale = 0.951;
      else                scale = 0.988;
    }
  } else {
    if (fabs(eta) < 1.45) {
      if      (cen < 10.) scale = 0.989;
      else if (cen < 30.) scale = 1.006;
      else                scale = 1.016;
    }else {
      if      (cen < 10.) scale = 0.967;
      else if (cen < 30.) scale = 1.018;
      else                scale = 1.054;
    }
  }

  newpt = pt * scale;

  float smear = 1.;

  if (isMC){
    if (fabs(eta) < 1.45) {
      if      (cen < 10.) smear = 1.19113;
      else if (cen < 30.) smear = 1.28262;
      else                smear = 2.20549;
    } else {
      if      (cen < 10.) smear = 3.13988;
      else if (cen < 30.) smear = 3.15340;
      else                smear = 3.17194;
    }

    newpt = newpt * smearRand->Gaus(1., smear / 91.1876);
  }

  return newpt;

}

// index, ntks in svtx, m svtx, csv
typedef std::tuple<int,int,float,float,TLorentzVector,int,int> BtagInfo_t;
static bool orderByBtagInfo(const BtagInfo_t &a, const BtagInfo_t &b)
{
  //int ntks_a(std::get<1>(a)), ntks_b(std::get<1>(b));
  //if(ntks_a>ntks_b) return true;

  float csv_a(std::get<3>(a)), csv_b(std::get<3>(b));
  if(csv_a>csv_b) return true;
  return false;
}

// btag efficiencies from the AN
float btagEfficiencies(int flavor, float cenbin){
  if(csvWP==0){
    if      (fabs(flavor) == 5) return (cenbin <= 30 ? 0.556 : 0.683); // bs
    else if (fabs(flavor) == 0) return (cenbin <= 30 ? 0.057 : 0.042); // unmatched
    else                        return (cenbin <= 30 ? 0.023 : 0.008); // udsg
  }
  else{
    if      (fabs(flavor) == 5) return (cenbin <= 30 ? 0.385 : 0.546); // bs
    else if (fabs(flavor) == 0) return (cenbin <= 30 ? 0.013 : 0.008); // unmatched
    else                        return (cenbin <= 30 ? 0.005 : 0.001); // udsg
  }
}

//
static bool orderByPt(const LeptonSummary &a, const LeptonSummary &b)
{
  float pt_a(a.p4.Pt()), pt_b(b.p4.Pt());
  if(pt_a>pt_b) return true;
  return false;
}

//parameterizations from Glauber MC
Double_t findNcoll(int hiBin) {
  const int nbins = 200;
  const Double_t Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
  return Ncoll[hiBin];
};

//
TF1 *getRBW(float m,float g) {
  //define the relativistic Breit-Wigner function
  TF1 *bwigner=new TF1("bwigner",
                       "[0]*([1]*[2]*sqrt([1]*[1]*([1]*[1]+[2]*[2]))/sqrt([1]*[1]+sqrt([1]*[1]*([1]*[1]+[2]*[2]))))/(TMath::Power(x*x-[1]*[1],2)+TMath::Power([1]*[2],2))",
                       max(float(0.),m-50*g),m+50*g);

  bwigner->SetParName(0,"N");
  bwigner->SetParameter(0,1.0);
  bwigner->SetParName(1,"m_{0}");
  bwigner->FixParameter(1,m);
  bwigner->SetParName(2,"#Gamma_{t}");
  bwigner->FixParameter(2,g);

  return bwigner;
}

//
float weightBW(TF1 *bwigner,std::vector<float> obsm,float g,float m,float gini,float mini) {

  //
  if(bwigner==NULL || obsm.size()!=2 || obsm[0]<150 || obsm[1]<150) return 1.;

  bwigner->FixParameter(1,mini);
  bwigner->FixParameter(2,gini);
  float nini=bwigner->Integral(max(m-50*g,float(0.)),m+50*g);
      
  bwigner->FixParameter(1,m);
  bwigner->FixParameter(2,g);
  float n=bwigner->Integral(max(m-50*g,float(0.)),m+50*g);

  float wgt(1.0);
  for(auto obsm_i : obsm){
    bwigner->FixParameter(1,mini);
    bwigner->FixParameter(2,gini);
    float vini=bwigner->Eval(obsm_i);

    bwigner->FixParameter(1,m);
    bwigner->FixParameter(2,g);
    float v=bwigner->Eval(obsm_i);

    wgt *= (v/n) / (vini/nini);
  }

  return wgt;
}



//
int main(int argc, char* argv[])
{
  TF1 * quenchingModel = new TF1("quenchingModel", "[0]/(TMath::Sqrt(2.*TMath::Pi())*0.73*x)*TMath::Exp(-1.*TMath::Power(TMath::Log(x/[0])+1.5,2)/2./0.73/0.73)", 0., 50.);
  TF1 * centralityModel = new TF1("centralityModel", "gaus");
  centralityModel->SetParameter(0,  1.090);
  centralityModel->SetParameter(1, -0.144);
  centralityModel->SetParameter(2,  0.442);

  bool blind(false);
  TString inURL,outURL;
  bool isMC(false),isPP(false),isAMCATNLO(false),isSkim(false);
  int maxEvents(-1);
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--in")!=string::npos && i+1<argc)         { inURL=TString(argv[i+1]); i++;}
    else if(arg.find("--out")!=string::npos && i+1<argc)   { outURL=TString(argv[i+1]); i++;}
    else if(arg.find("--max")!=string::npos && i+1<argc)   { sscanf(argv[i+1],"%d",&maxEvents); }
    else if(arg.find("--csvWP")!=string::npos && i+1<argc) { sscanf(argv[i+1],"%d",&csvWP); }
    else if(arg.find("--mc")!=string::npos)                { isMC=true;  }
    else if(arg.find("--pp")!=string::npos)                { isPP=true;  }
    else if(arg.find("--amcatnlo")!=string::npos)          { isAMCATNLO=true;  }
    else if(arg.find("--skim")!=string::npos)              { isSkim=true;  }
  }
  
  bool isSingleMuPD( !isMC && inURL.Contains("SkimMuons"));
  bool isSingleElePD( !isMC && inURL.Contains("SkimElectrons"));
  bool isMuSkimedMCPD( isMC && inURL.Contains("HINPbPbAutumn18DR_skims") && inURL.Contains("Muons"));
  bool isEleSkimedMCPD( isMC && inURL.Contains("HINPbPbAutumn18DR_skims") && inURL.Contains("Electrons"));
  LumiRun lumiTool;
  ElectronEfficiencyWrapper eleEff("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/data", false);

  //read expected trigger efficiencies
  TString trigEffURL("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/data/trigeff_mc.root");
  gSystem->ExpandPathName(trigEffURL);
  TFile *fIn=TFile::Open(trigEffURL);
  TGraphAsymmErrors *e_mctrigeff=(TGraphAsymmErrors *)fIn->Get("e_pt_trigeff");
  //TGraphAsymmErrors *m_mctrigeff=(TGraphAsymmErrors *)fIn->Get("m_eta_trigeff");
  fIn->Close();

  //read isolation efficiencies
  TString isosfURL("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/data/isolation_sf.root");
  gSystem->ExpandPathName(isosfURL);
  fIn=TFile::Open(isosfURL);
  std::map<TString, TH2 *> isoEffSFs;
  TString isomaps[]={"cen_169","periph_169","cen_121","periph_121"};
  for(size_t i=0; i<sizeof(isomaps)/sizeof(TString); i++){
    TString key(isomaps[i]);
    isoEffSFs[key]=(TH2*)fIn->Get("sfiso2eff_"+key);
    isoEffSFs[key]->SetDirectory(0);
  }
  fIn->Close();
  
  // initialize the JEC and associated unc files
  std::vector<std::string> FilesData;
  TString DATA_L2RelativeURL("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/data/Autumn18_HI_V6_DATA_L2Relative_AK4PF.txt");
  gSystem->ExpandPathName(DATA_L2RelativeURL);
  TString DATA_L2L3ResidualURL("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/data/Autumn18_HI_V6_DATA_L2L3Residual_AK4PF.txt");
  gSystem->ExpandPathName(DATA_L2L3ResidualURL);
  FilesData.push_back(DATA_L2RelativeURL.Data());
  FilesData.push_back(DATA_L2L3ResidualURL.Data());
  
  JetCorrector JECData(FilesData);
  TString JEUDataURL("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/data/Autumn18_HI_V6_DATA_Uncertainty_AK4PF.txt");
  gSystem->ExpandPathName(DATA_L2RelativeURL);
  JetUncertainty JEUData(DATA_L2RelativeURL.Data());
  
  std::vector<std::string> FilesMC;
  TString FilesMCURL("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/data/Autumn18_HI_V6_MC_L2Relative_AK4PF.txt");
  gSystem->ExpandPathName(FilesMCURL);
  
  FilesMC.push_back(FilesMCURL.Data());
  
  JetCorrector JECMC(FilesMC);
  TString JECMCURL("${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/data/Autumn18_HI_V6_MC_Uncertainty_AK4PF.txt");
  gSystem->ExpandPathName(JECMCURL);
  JetUncertainty JEUMC(JECMCURL.Data());
  
  
  if(isPP)
    cout << "Treating as a pp collision file" << endl;
  if(isMC)
    cout << "Treating as a MC file" << endl;
  if(isAMCATNLO)
    cout << "This is an amc@NLO file (for the ME weights)" << endl;
  else
    cout << "Will assume that this is a powheg file!" << endl;

  
  std::vector<size_t> meIdxList={0,1,2,3,4,6,8}; //qcd weights
  if(!isMC) 
    meIdxList.clear();
  else {
    if(isPP){
      for(size_t i=10; i<=111; i++) meIdxList.push_back(i); //hessian NNPDF3.1
      meIdxList.push_back(116); meIdxList.push_back(117);   //alphaS variation +/-0.001
    }else if(isAMCATNLO){    
      meIdxList[0]=1080;
      for(size_t i=1081; i<=1176; i++) meIdxList.push_back(i); //EPPS16nlo_CT14nlo_Pb208
      for(size_t i=1177; i<=1209; i++) meIdxList.push_back(i); //nCTEQ15FullNuc_208_82
    }else{    
      meIdxList[0]=864;
      for(size_t i=865; i<=960; i++) meIdxList.push_back(i); //EPPS16nlo_CT14nlo_Pb208
      for(size_t i=961; i<=993; i++) meIdxList.push_back(i); //nCTEQ15FullNuc_208_82
    }
    cout << "Will store " <<  meIdxList.size() << " ME weights" << endl;
  }

  //for breit-wigner reweithing
  TF1 *rbwigner( isMC ? getRBW(172.5,1.31) : NULL );

  //book some histograms
  HistTool ht;
  if(isMC){
    ht.addHist("fidcounter",  new TH2F("fidcounter", ";Fiducial counter;Events",5,0,5,meIdxList.size(),0,meIdxList.size()));
    ht.get2dPlots()["fidcounter"]->GetXaxis()->SetBinLabel(1,"all");
    ht.get2dPlots()["fidcounter"]->GetXaxis()->SetBinLabel(2,"=2l");
    ht.get2dPlots()["fidcounter"]->GetXaxis()->SetBinLabel(3,"=2l fid");
    ht.get2dPlots()["fidcounter"]->GetXaxis()->SetBinLabel(4,"=2l,#geq1b fid");
    ht.get2dPlots()["fidcounter"]->GetXaxis()->SetBinLabel(5,"=2l,#geq2b fid");
  }else{
    ht.addHist("ratevsrun",lumiTool.getLumiMonitor());
  }

  //generic histograms
  ht.addHist("br",  new TH1F("br",    ";BR;Events",6,0,6));
  ht.addHist("trig_pt",  new TH1F("trig_pt",    ";Lepton transverse momentum [GeV];Events",20,20,200));
  ht.addHist("trig_eta", new TH1F("trig_eta",   ";Lepton pseudo-rapidity;Events",20,0,2.5));
  for(int i=0; i<2; i++) {
    TString pf(Form("l%d",i+1));
    ht.addHist(pf+"pt",             new TH1F(pf+"pt",            ";Lepton transverse momentum [GeV];Events",20,20,200));
    ht.addHist(pf+"eta",            new TH1F(pf+"eta",           ";Lepton pseudo-rapidity;Events",20,0,2.5));

    for(int j=0; j<3; j++) {
      TString comp("ch");
      if(j==1) comp="pho";
      if(j==2) comp="nh";
      ht.addHist(pf+comp+"iso",      new TH1F(pf+comp+"iso",      ";PF "+comp+" isolation;Leptons",50,0,250));
      ht.addHist(pf+comp+"reliso",   new TH1F(pf+comp+"reliso",   ";Relative PF "+comp+" isolation;Leptons",50,0,2.0));
      ht.addHist(pf+comp+"isovscen", new TH2F(pf+comp+"isovscen", ";Centrality bin;PF "+comp+" isolation [GeV];Leptons",10,0,100,50,0,100));
      ht.addHist(pf+comp+"isovsrho", new TH2F(pf+comp+"isovsrho", ";#rho_{"+comp+"};PF "+comp+" isolation [GeV];Leptons",10,0,100,20,0,100));
    }    
  }

  //electron specific
  ht.addHist("esihih",       new TH1F("esihih",      ";#sigma(i#etai#eta);Electrons",       50,0,0.06));
  ht.addHist("edetaseedvtx", new TH1F("edetaseedvtx",    ";#Delta#eta(vtx);Electrons",          50,0,0.015));
  ht.addHist("edphivtx",     new TH1F("edphivtx",    ";#Delta#phi(vtx) [rad];Electrons",    50,0,0.015));
  ht.addHist("ehoebc",       new TH1F("ehoebc"    ,    ";h/e;Electrons",                      50,0,0.25));
  ht.addHist("eempinv",      new TH1F("eempinv",     ";|1/E-1/p| [1/GeV];Electrons",        50,0,0.05));
  ht.addHist("ed0",          new TH1F("ed0",         ";d_{0} [cm];Electrons",               50,0,0.05));
  ht.addHist("edz",          new TH1F("edz",         ";d_{z} [cm];Electrons",               50,0,0.05));
  ht.addHist("emll",         new TH1F("emll",        ";Di-electron invariant mass [GeV];Events",  40,20,200));

  //muon specific
  ht.addHist("mmusta",     new TH1F("mmusta",      ";Muon stations;Muons",            15,0,15));   
  ht.addHist("mtrklay",    new TH1F("mtrklay",     ";Tracker layers;Muons",           25,0,25));
  ht.addHist("mchi2ndf",   new TH1F("mchi2ndf",    ";#chi^2/ndf;Muons",               50,0,15));
  ht.addHist("mmuhits",    new TH1F("mmuhits",     ";Muon hits;Muons",                25,0,25));
  ht.addHist("mpxhits",    new TH1F("mpxhits",     ";Pixel hits;Muons",               15,0,15));
  ht.addHist("md0",        new TH1F("md0",         ";d_{0} [cm];Muons",               50,0,0.5));
  ht.addHist("mdz",        new TH1F("mdz",         ";d_{z} [cm];Muons",               50,0,1.0));
  ht.addHist("mmll",       new TH1F("mmll",        ";Di-muon invariant mass [GeV];Events",  40,20,200));
 
  ht.addHist("mll",      new TH1F("mll",      ";Dilepton invariant mass [GeV];Events",40,0,200));
  ht.addHist("ptll",     new TH1F("ptll",     ";Dilepton transverse momentum [GeV];Events",25,0,200));
  ht.addHist("ptsum",    new TH1F("ptsum",    ";p_{T}(l)+p_{T}(l') [GeV];Events",25,0,200));
  ht.addHist("acopl",    new TH1F("acopl" ,   ";1-#Delta#phi(l,l')/#pi;Events",20,0,1.0));
  ht.addHist("detall",   new TH1F("detall",   ";#Delta#eta(l,l');Events",20,0,4));
  ht.addHist("drll",     new TH1F("drll"  ,   ";#DeltaR(l,l');Events",20,0,2*TMath::Pi()));

  ht.addHist("pfrapavg",      new TH1F("pfrapavg",     ";Average rapidity;Events",25,0,2.5));
  ht.addHist("pfraprms",      new TH1F("pfraprms",     ";#sigma(rapidity);Events",50,0,2.5));
  ht.addHist("pfrapmaxspan",  new TH1F("pfrapmaxspan", ";Maximum rapidity span;Events",25,0,5));
  ht.addHist("pfht",          new TH1F("pfht",         ";H_{T} [GeV];Events",25,0,500));
  ht.addHist("pfmht",         new TH1F("pfmht",        ";Missing H_{T} [GeV];Events",25,0,200));
  ht.addHist("npfjets",    new TH1F("npfjets",    ";Jet multiplicity;Events",8,0,8));
  ht.addHist("npfbjets",   new TH1F("npfbjets",   ";b-jet multiplicity;Events",5,0,5));    
  ht.addHist("npfsvtx",    new TH1F("npfsvtx",    ";Secondary vertex multiplicity;Events",5,0,5));

  ht.addHist("jetptprequench" ,    new TH1F("jetptprequench" , ";jets pT pre quenching" ,25,0,100.));
  ht.addHist("jetptpostquench",    new TH1F("jetptpostquench", ";jets pT post quenching",25,0,100.));
  ht.addHist("jetquenchloss"  ,    new TH1F("jetquenchloss"  , ";pT loss in quenching"  ,10,0,30.));

  ht.addHist("jetptprequenchB" ,    new TH1F("jetptprequenchB" , ";matched b: jets pT pre quenching" ,25,0,100.));
  ht.addHist("jetptpostquenchB",    new TH1F("jetptpostquenchB", ";matched b: jets pT post quenching",25,0,100.));
  ht.addHist("jetquenchlossB"  ,    new TH1F("jetquenchlossB"  , ";matched b: pT loss in quenching"  ,10,0,30.));

  for(size_t j=1; j<=2; j++){
    TString ppf(j==1 ? "1" : "2");
    ht.addHist("pf"+ppf+"jbalance",    new TH1F("pf"+ppf+"jbalance", ";R = p_{T}(j)/p_{T}(ll);Events",50,0,3));
    ht.addHist("pf"+ppf+"jpt",         new TH1F("pf"+ppf+"jpt",      ";Jet transverse momentum [GeV];Events",30,00,300));
    ht.addHist("pf"+ppf+"jeta",        new TH1F("pf"+ppf+"jeta",     ";Jet pseudo-rapidity;Events",20,0,2.5));
    ht.addHist("pf"+ppf+"jetavsphi",   new TH2F("pf"+ppf+"jetavsphi", ";Jet pseudo-rapidity;Jet azimuthal angle [rad];Events",100,-2.5,2.5,100,-TMath::Pi(),TMath::Pi()));
    ht.addHist("pf"+ppf+"jsvtxm",      new TH1F("pf"+ppf+"jsvtxm",   ";Secondary vertex mass;Events",25,0,6));
    ht.addHist("pf"+ppf+"jsvtxntk",    new TH1F("pf"+ppf+"jsvtxntk", ";Secondary vertex track multiplicity;Events",5,0,5));
    ht.addHist("pf"+ppf+"jcsv",        new TH1F("pf"+ppf+"jcsv",     ";CSVv2;Events",25,0,1));
  }
  ht.addHist("jptvsjptquench",  new TH2F("jptvsjptquench", ";Reconstructed jet p_{T} [GeV];Quenched jet p_{T} [GeV];Jets",50,0,100,50,0,100) );

  // Initialize the btagging SF stuff
  BTagSFUtil * myBTagUtil = new BTagSFUtil(42);
  TRandom3 * rand = new TRandom3(2);
  cout << __LINE__ << endl;
  //Get Tree info
  char* read = new char[100];
  // TChain *hiInfoTree_p  = new TChain("HiForest/HiForestInfo");
  TChain *hiInfoTree_p  = new TChain("HiForestInfo/HiForest");
  hiInfoTree_p->Add(inURL);
  TBranch *branch = hiInfoTree_p->GetBranch("GlobalTag");
  branch->SetAddress((void*)read);
  hiInfoTree_p->GetEntry(0);
  string GT(read, 0, 100);

  //Get global event filters 
  TChain *globalTree_p     = new TChain("skimanalysis/HltTree");
  globalTree_p->Add(inURL);
  ForestSkim fForestSkim(globalTree_p);

  //configure leptons
  // TString lepTreeName("ggHiNtuplizerGED/EventTree");
  TString lepTreeName("ggHiNtuplizer/EventTree");
  if(isPP) lepTreeName="ggHiNtuplizer/EventTree";
  if(GT.find("75X_mcRun2")!=string::npos) lepTreeName="ggHiNtuplizer/EventTree";
  TChain *lepTree_p     = new TChain(lepTreeName);
  lepTree_p->Add(inURL);
  ForestLeptons fForestLep(lepTree_p);
  ForestGen fForestGen(lepTree_p);

  //configure PF cands
  // TChain *pfCandTree_p  = new TChain("pfcandAnalyzer/pfTree");
  TChain *pfCandTree_p  = new TChain("particleFlowAnalyser/pftree");
  pfCandTree_p->Add(inURL);
  ForestPFCands fForestPF(pfCandTree_p);

  //configure jets
  // TChain *jetTree_p     = new TChain(isPP ? "ak4PFJetAnalyzer/t" : "akFlowPuCs4PFJetAnalyzer/t");
  TChain *jetTree_p     = new TChain(isPP ? "ak4PFJetAnalyzer/t" : "akCs4PFJetAnalyzer/t");
  jetTree_p->Add(inURL);
  ForestJets fForestJets(jetTree_p);

  //global variables
  TChain *hiTree_p      = new TChain("hiEvtAnalyzer/HiTree");
  hiTree_p->Add(inURL);
  HiTree fForestTree(hiTree_p);

  //trigger
  TChain *hltTree_p     = new TChain("hltanalysis/HltTree");
  hltTree_p->Add(inURL);
  int etrig(0),mtrig(0);
  TString muTrigName(""),eTrigName("");
  if(isPP){
    muTrigName="HLT_HIL3Mu20_v";
    if( !hltTree_p->FindBranch(muTrigName+"1") ) muTrigName="HLT_HIL3Mu15ForPPRef_v";
    if( GT.find("75X_mcRun2")!=string::npos ) muTrigName="HLT_HIL2Mu15ForPPRef_v";
    hltTree_p->SetBranchStatus(muTrigName+"1",1);
    hltTree_p->SetBranchAddress(muTrigName+"1",&mtrig);
    eTrigName="HLT_HIEle20_WPLoose_Gsf_v";
    if( !hltTree_p->FindBranch(eTrigName+"1") ) eTrigName="HLT_HISinglePhoton20_Eta3p1ForPPRef_v";
    if( GT.find("75X_mcRun2")!=string::npos ) eTrigName="HLT_HISinglePhoton40_Eta3p1ForPPRef_v";
    hltTree_p->SetBranchStatus(eTrigName+"1",1);
    hltTree_p->SetBranchAddress(eTrigName+"1",&etrig);
  }else{
    // muTrigName="HLT_HIL3Mu12_v";
    muTrigName="HLT_HIL3SingleMu12_v";
    eTrigName="HLT_HIEle20Gsf_v";
    if(isSingleMuPD || isMuSkimedMCPD) mtrig = 1;
    if(isSingleElePD || isEleSkimedMCPD) etrig = 1;
    if(isMC and !isMuSkimedMCPD and !isEleSkimedMCPD){
      hltTree_p->SetBranchStatus(muTrigName+"1",1);
      hltTree_p->SetBranchAddress(muTrigName+"1",&mtrig);
      hltTree_p->SetBranchStatus(eTrigName+"1",1);
      hltTree_p->SetBranchAddress(eTrigName+"1",&etrig);
    }
  }

  //trigger objects
  cout << "Using " << muTrigName << " " << eTrigName << " as triggers" << endl;
  TChain *muHLTObj_p =NULL, *eleHLTObj_p=NULL;
  ForestHLTObject *muHLTObjs=NULL, *eleHLTObjs=NULL;
  if(!isPP){
    muHLTObj_p=new TChain("hltobject/"+muTrigName);
    muHLTObj_p->Add(inURL);
    muHLTObjs=new ForestHLTObject(muHLTObj_p);
    eleHLTObj_p = new TChain("hltobject/"+eTrigName);
    eleHLTObj_p->Add(inURL);
    eleHLTObj_p->AddFriend(muHLTObj_p);
    eleHLTObjs=new ForestHLTObject(eleHLTObj_p);
  }

  // TChain *rhoTree_p = new TChain("hiFJRhoAnalyzerFinerBins/t");
  // rhoTree_p->Add(inURL);
  // std::vector<Double_t> *t_rho=0,*t_rhom=0,*t_etaMin=0,*t_etaMax=0;
  // if(rhoTree_p){
  //   rhoTree_p->SetBranchAddress("rho", &t_rho);
  //   rhoTree_p->SetBranchAddress("etaMin", &t_etaMin);
  //   rhoTree_p->SetBranchAddress("etaMax", &t_etaMax);
  //   rhoTree_p->SetBranchAddress("rhom", &t_rhom);
  // }else{
  //   std::cout << "[WARN] Can't find rho tree hiFJRhoAnalyzerFinerBins/t" << std::endl;
  // }
  
  // =============================================================
  // marc here make the output tree

  TTree * outTree = new TTree("tree", "tree with 2lepton selection and combined collections");

  // event and trigger variables
  Int_t  t_run, t_lumi, t_etrig, t_mtrig, t_isData;
  Long_t t_event;
  Float_t t_vx, t_vy, t_vz, t_weight, t_weight_BRW(1.0), t_cenbin, t_ncollWgt, t_trigSF, t_trigSFUnc;
  std::vector<Float_t> t_meWeights;
  outTree->Branch("run"   , &t_run  , "run/I");
  outTree->Branch("lumi"  , &t_lumi , "lumi/I");
  outTree->Branch("event" , &t_event, "event/L");
  outTree->Branch("isData", &t_isData, "isData/I");

  outTree->Branch("vx", &t_vx, "vx/F");
  outTree->Branch("vy", &t_vy, "vy/F");
  outTree->Branch("vz", &t_vz, "vz/F");

  outTree->Branch("weightBRW", &t_weight_BRW, "weightBRW/F");
  outTree->Branch("weight", &t_weight, "weight/F");
  outTree->Branch("meWeights", &t_meWeights);


  // centrality and different flavors of rho
  outTree->Branch("cenbin", &t_cenbin, "cenbin/F");
  outTree->Branch("ncollWgt", &t_ncollWgt, "ncollWgt/F");
  
  // outTree->Branch("rho",    &t_rho);
  // outTree->Branch("rhom",   &t_rhom);

  Float_t t_globalrho;
  outTree->Branch("globalrho", &t_globalrho, "globalrho/F");

  outTree->Branch("etrig" , &t_etrig , "etrig/I");
  outTree->Branch("mtrig" , &t_mtrig , "mtrig/I");

  outTree->Branch("trigSF" ,   &t_trigSF ,     "trigSF/F");
  outTree->Branch("trigSFUnc" , &t_trigSFUnc , "trigSFUnc/F");

  // variables per lepton, including iso
  Int_t t_nlep, t_lep_ind1, t_lep_ind2;
  std::vector<Float_t> t_lep_pt, t_lep_calpt, t_lep_eta, t_lep_phi, t_lep_d0, t_lep_dz, t_lep_d0err, t_lep_phiso, t_lep_chiso, t_lep_nhiso, t_lep_rho, t_lep_isofull, t_lep_miniiso,t_lep_isofull20,t_lep_isofull25,t_lep_isofull30;
  std::vector<Bool_t> t_lep_matched,t_lep_taufeeddown,t_lep_trigmatch;
  std::vector<Int_t  > t_lep_pdgId, t_lep_charge,t_lep_idflags;
  std::vector<Float_t> t_lepSF,t_lepSFUnc,t_lepIsoSF,t_lepIsoSFUnc;
  outTree->Branch("nlep"       , &t_nlep      , "nlep/I"            );
  outTree->Branch("lep_ind1"   , &t_lep_ind1  , "lep_ind1/I");
  outTree->Branch("lep_ind2"   , &t_lep_ind2  , "lep_ind2/I");
  outTree->Branch("lep_pt"     , &t_lep_pt);
  outTree->Branch("lep_calpt"     , &t_lep_calpt);
  outTree->Branch("lep_eta"    , &t_lep_eta);
  outTree->Branch("lep_phi"    , &t_lep_phi);
  outTree->Branch("lep_idflags", &t_lep_idflags);
  outTree->Branch("lep_d0"    , &t_lep_d0);
  outTree->Branch("lep_d0err"  , &t_lep_d0err);
  outTree->Branch("lep_dz"     , &t_lep_dz);
  outTree->Branch("lep_phiso"  , &t_lep_phiso);
  outTree->Branch("lep_chiso"  , &t_lep_chiso);
  outTree->Branch("lep_nhiso"  , &t_lep_nhiso);
  outTree->Branch("lep_rho"    , &t_lep_rho);
  outTree->Branch("lep_pdgId"  , &t_lep_pdgId);
  outTree->Branch("lep_charge" , &t_lep_charge);
  outTree->Branch("lep_isofull", &t_lep_isofull);
  outTree->Branch("lep_isofull20", &t_lep_isofull20);
  outTree->Branch("lep_isofull25", &t_lep_isofull25);
  outTree->Branch("lep_isofull30", &t_lep_isofull30);
  outTree->Branch("lep_miniiso", &t_lep_miniiso);
  outTree->Branch("lep_matched", &t_lep_matched);
  outTree->Branch("lep_taufeeddown", &t_lep_taufeeddown);
  outTree->Branch("lep_trigmatch", &t_lep_trigmatch);
  outTree->Branch("lepSF",      &t_lepSF);
  outTree->Branch("lepSFUnc",      &t_lepSFUnc);
  outTree->Branch("lepIsoSF",      &t_lepIsoSF);
  outTree->Branch("lepIsoSFUnc",      &t_lepIsoSFUnc);

  // variables from dilepton system
  Float_t t_zpt(-1),t_llpt, t_llpt_raw, t_lleta, t_llphi, t_llm, t_llm_raw, t_dphi, t_deta, t_sumeta;
  outTree->Branch("zpt"    , &t_zpt   , "zpt/F");
  outTree->Branch("llpt"   , &t_llpt   , "llpt/F");
  outTree->Branch("llpt_raw", &t_llpt_raw   , "llpt_raw/F");
  outTree->Branch("lleta"  , &t_lleta  , "lleta/F");
  outTree->Branch("llphi"  , &t_llphi  , "llphi/F");
  outTree->Branch("llm"    , &t_llm    , "llm/F");
  outTree->Branch("llm_raw", &t_llm_raw, "llm_raw/F");
  outTree->Branch("dphi"   , &t_dphi   , "dphi/F");
  outTree->Branch("deta"   , &t_deta   , "deta/F");
  outTree->Branch("sumeta" , &t_sumeta , "sumeta/F");

  // variables per bjet (jets ordered by csvv2)
  Int_t t_nbjet;
  Bool_t t_bjet_leadPassTight;
  std::vector<Float_t> t_bjet_pt, t_bjet_eta, t_bjet_phi, t_bjet_mass, t_bjet_csvv2;
  outTree->Branch("nbjet"      , &t_nbjet      , "nbjet/I"            );
  outTree->Branch("bjet_leadPassTight"    , &t_bjet_leadPassTight);
  outTree->Branch("bjet_pt"    , &t_bjet_pt    );
  outTree->Branch("bjet_eta"   , &t_bjet_eta   );
  outTree->Branch("bjet_phi"   , &t_bjet_phi   );
  outTree->Branch("bjet_mass"  , &t_bjet_mass  );
  outTree->Branch("bjet_csvv2" , &t_bjet_csvv2 );

  Int_t t_nbjet_sel, t_nbjet_sel_jecup, t_nbjet_sel_jecdn, t_nbjet_sel_jerup, t_nbjet_sel_jerdn, t_nbjet_sel_bup, t_nbjet_sel_bdn, t_nbjet_sel_udsgup, t_nbjet_sel_udsgdn, t_nbjet_sel_quenchup, t_nbjet_sel_quenchdn;
  outTree->Branch("nbjet_sel"       , &t_nbjet_sel       , "nbjet_sel/I"       );
  outTree->Branch("nbjet_sel_jecup" , &t_nbjet_sel_jecup , "nbjet_sel_jecup/I" );
  outTree->Branch("nbjet_sel_jecdn" , &t_nbjet_sel_jecdn , "nbjet_sel_jecdn/I" );
  outTree->Branch("nbjet_sel_jerup" , &t_nbjet_sel_jerup , "nbjet_sel_jerup/I" );
  outTree->Branch("nbjet_sel_jerdn" , &t_nbjet_sel_jerdn , "nbjet_sel_jerdn/I" );
  outTree->Branch("nbjet_sel_bup"   , &t_nbjet_sel_bup   , "nbjet_sel_bup/I"   );
  outTree->Branch("nbjet_sel_bdn"   , &t_nbjet_sel_bdn   , "nbjet_sel_bdn/I"   );
  outTree->Branch("nbjet_sel_udsgup", &t_nbjet_sel_udsgup, "nbjet_sel_udsgup/I");
  outTree->Branch("nbjet_sel_udsgdn", &t_nbjet_sel_udsgdn, "nbjet_sel_udsgdn/I");
  outTree->Branch("nbjet_sel_quenchup", &t_nbjet_sel_quenchup, "nbjet_sel_quenchup/I");
  outTree->Branch("nbjet_sel_quenchdn", &t_nbjet_sel_quenchdn, "nbjet_sel_quenchdn/I");

  std::vector<Float_t> t_bjet_matchpt, t_bjet_matcheta, t_bjet_matchphi, t_bjet_matchmass;
  outTree->Branch("bjet_genpt"    , &t_bjet_matchpt    );
  outTree->Branch("bjet_geneta"   , &t_bjet_matcheta   );
  outTree->Branch("bjet_genphi"   , &t_bjet_matchphi   );
  outTree->Branch("bjet_genmass"  , &t_bjet_matchmass  );

  std::vector<Int_t> t_bjet_flavor, t_bjet_flavorForB;
  outTree->Branch("bjet_flavor"  , &t_bjet_flavor  );
  outTree->Branch("bjet_flavorB"  , &t_bjet_flavorForB  );

  // constructed variables like ht and stuff
  Float_t t_ht, t_mht, t_apt, t_dphilll2;
  outTree->Branch("ht"     , &t_ht     , "ht/F");
  outTree->Branch("mht"    , &t_mht    , "mht/F");
  outTree->Branch("apt"    , &t_apt    , "apt/F");
  outTree->Branch("dphilll2"    , &t_dphilll2    , "dphilll2/F");

  // Float_t t_bdt, t_bdt_rarity, t_fisher2;
  // outTree->Branch("bdt", &t_bdt, "bdt/F");
  // outTree->Branch("bdtrarity", &t_bdt_rarity, "bdtrarity/F");
  // outTree->Branch("fisher2", &t_fisher2, "fisher2/F");
  // // =============================================================
  // //
  // TMVA::Tools::Instance();
  // TMVA::Reader *reader        = new TMVA::Reader( "!Color:!Silent" );
  // TMVA::Reader *readerFisher2 = new TMVA::Reader( "!Color:!Silent" );

  // // make new variables because i'm too lazy to think right now
  // Float_t bdt_l1pt, bdt_apt, bdt_abslleta, bdt_dphilll2, bdt_sumabseta;//, bdt_flavor;

  // // these must have the same name as in the training. and the same ordeeeeeeer
  // // copy directly from the script that runs the training:
  // //dataloader.AddVariable('lep_pt[0]'  , 'p_{T}^{lep1}'     , 'GeV' , 'F')
  // //dataloader.AddVariable('apt'        , 'A_{pt}'           , ''    , 'F')
  // //dataloader.AddVariable('llpt'       , 'p_{T}^{ll}'       , 'GeV' , 'F')
  // //dataloader.AddVariable('abs(lleta)' , '|#eta^{ll}|'      , ''    , 'F')
  // //dataloader.AddVariable('dphi'       , '|#Delta #phi|'    , 'rad' , 'F')
  // //dataloader.AddVariable('abs(lep_eta[0])+abs(lep_eta[1])' , '#sum |#eta_{i}|', ''    , 'F')
  // //dataloader.AddVariable('abs(lep_pdgId[0]*lep_pdgId[1])'  , 'flavor', ''    , 'F')
  // //
  // reader->AddVariable("lep_pt[0]"  , &bdt_l1pt    );
  // reader->AddVariable("apt"        , &bdt_apt     );
  // reader->AddVariable("llpt"       , &t_llpt      );
  // reader->AddVariable("abs(lleta)" , &bdt_abslleta);
  // reader->AddVariable("abs(dphi)"  , &t_dphi      );
  // reader->AddVariable("abs(lep_eta[0])+abs(lep_eta[1])", &bdt_sumabseta);

  // // for the fisher just take these two
  // //dataloader.AddVariable('llpt'       , 'p_{T}^{ll}'       , 'GeV' , 'F')
  // //dataloader.AddVariable('dphi'       , '|#Delta #phi|'    , 'rad' , 'F')
  // readerFisher2->AddVariable("llpt", &t_llpt);
  // readerFisher2->AddVariable("abs(dphi)", &t_dphi);

  // TString methodName       ("BDTG method");
  // TString methodNameFisher2("Fisher method");
  // // hard coded path for now ...
  // TString weightFile("/afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_9_4_6_patch1/src/HeavyIonsAnalysis/topskim/scripts/training_dy/weights/TMVAClassification_BDTG.weights.xml");
  // reader->BookMVA( methodName, weightFile);

  // TString weightFileFisher2("/afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_9_4_6_patch1/src/HeavyIonsAnalysis/topskim/scripts/training_dy_fisher2/weights/TMVAClassification_Fisher.weights.xml");
  // readerFisher2->BookMVA( methodNameFisher2, weightFileFisher2);

    
  Double_t wgtSum(0);
  std::vector<Double_t> allWgtSum;
  int nEntries = (int)lepTree_p->GetEntries();  
  int entryDiv = ((int)(nEntries/20));    
  cout << inURL << " has " << nEntries << " events to process" << endl;
  if(maxEvents>0) { 
    nEntries=TMath::Min(nEntries,maxEvents); 
    cout << "Number of events to process limited to " << nEntries << endl;
  }

  //get ncoll weighting norm factor
  double ncollWgtNorm(1.0);
  if(isMC && !isPP){
    double ncollSum(0.);
    for(int entry = 0; entry < nEntries; entry++){
      hiTree_p->GetEntry(entry);
      ncollSum+=findNcoll(fForestTree.hiBin);
    }
    if(ncollSum>0) ncollWgtNorm=double(nEntries)/ncollSum;
  }

  //loop over events
  for(int entry = 0; entry < nEntries; entry++){
    
    if(entryDiv!=0)if(entry%entryDiv == 0) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
    globalTree_p->GetEntry(entry);
    lepTree_p->GetEntry(entry);
    pfCandTree_p->GetEntry(entry);
    jetTree_p->GetEntry(entry);    
    hltTree_p->GetEntry(entry);
    hiTree_p->GetEntry(entry);
    // if(rhoTree_p) rhoTree_p->GetEntry(entry);
    if(muHLTObj_p ) muHLTObj_p->GetEntry(entry);
    if(eleHLTObj_p)  eleHLTObj_p->GetEntry(entry);
    
    //gen level analysis
    float evWgt(1.0),topPtWgt(1.0),topMassUpWgt(1.0),topMassDnWgt(1.0);
    int genDileptonCat(1.);
    std::vector<TLorentzVector> genLeptons, genZLeptons, genBjets;
    std::vector<bool> genTauLeptons;
    std::vector<int> genLeptonIds;
    bool isGenDilepton(false),isLeptonFiducial(false),is1bFiducial(false),is2bFiducial(false);    
    if(isMC) {

      //gen level selection
      int nlFromTopW(0);
      TLorentzVector topP4(0,0,0,0),antitopP4(0,0,0,0);
      for(size_t i=0; i<fForestGen.mcPID->size(); i++) {
        int pid=fForestGen.mcPID->at(i);
        //int sta=fForestGen.mcStatus->at(i);
        int mom_pid=fForestGen.mcMomPID->at(i);
        int gmom_pid=fForestGen.mcGMomPID->at(i);
        
        TLorentzVector p4(0,0,0,0);
        p4.SetPtEtaPhiM( fForestGen.mcPt->at(i), fForestGen.mcEta->at(i), fForestGen.mcPhi->at(i), fForestGen.mcMass->at(i) );

        if(pid==6)  topP4=p4;
        if(pid==-6) antitopP4=p4;
        if( abs(pid)<6  && abs(mom_pid)==6 ) {          
          if(p4.Pt()>30 && fabs(p4.Eta())<2.5) genBjets.push_back(p4);          
        }
        
        bool isFromTop(abs(mom_pid)==6 || abs(gmom_pid)==6 );
        bool isFromZ( abs(mom_pid)==23 || abs(gmom_pid)==23 );
        bool isFromW( abs(mom_pid)==24 || abs(gmom_pid)==24 );
        bool isTauFeedDown( abs(mom_pid)==15 );
        
        //count W leptonic decays
        if( abs(pid)==11 || abs(pid)==13 )
          {
            if(isFromTop && isFromW) nlFromTopW++;
          }
        if(abs(pid)==16) //use tau neutrino here
          {
            if(isFromW && isTauFeedDown) nlFromTopW++;
          }

        //charged leptons
        if(abs(pid)==11 || abs(pid)==13) {

          //leptons from t->W->l or W->tau->l
          if(isFromW && (isTauFeedDown || isFromTop ) ) {
            genLeptons.push_back(p4);
            genTauLeptons.push_back(isTauFeedDown);
            genLeptonIds.push_back(pid);
            genDileptonCat *= abs(pid);
          }

          //leptons from Z->ll or Z->tt->ll
          if(isFromZ || isTauFeedDown){
            genZLeptons.push_back(p4);
          }
        }

        //neutrinos
        if(abs(pid)==12 || abs(pid)==14 || abs(pid)==16){
          if(isFromZ || isTauFeedDown) {
            genZLeptons.push_back(p4);
          }
        }
      }
      t_weight_BRW=getMadgraphBRWlCorrection(nlFromTopW);
      ht.fill("br",2*nlFromTopW,1);
      ht.fill("br",2*nlFromTopW+1,t_weight_BRW);

      TLorentzVector gendil(0,0,0,0);
      for(auto &l:genZLeptons) gendil += l;
      t_zpt=gendil.Pt();
      
      topPtWgt = TMath::Exp(0.199-0.00166*topP4.Pt());
      topPtWgt *= TMath::Exp(0.199-0.00166*antitopP4.Pt());
      topPtWgt = TMath::Sqrt(topPtWgt);
      std::vector<float> obsm={float(topP4.M()),float(antitopP4.M())};
      topMassDnWgt = weightBW(rbwigner,obsm,171.5,1.28,172.5,1.31);
      topMassUpWgt = weightBW(rbwigner,obsm,173.5,1.34,172.5,1.31);
      
      isGenDilepton=(genLeptons.size()==2);      
      isLeptonFiducial=(isGenDilepton && 
                        genLeptons[0].Pt()>lepPtCut && fabs(genLeptons[0].Eta())<muEtaCut && 
                        genLeptons[1].Pt()>lepPtCut && fabs(genLeptons[1].Eta())<muEtaCut);  
      
      //further cuts for electrons (EE-EB transition, HEM15/16 transition)
      if(isLeptonFiducial){
        for(size_t igl=0; igl<2; igl++){
          if(abs(genLeptonIds[igl])!=11) continue;
          float eta(genLeptons[igl].Eta());
          float phi(genLeptons[igl].Phi());
          if(fabs(eta) > barrelEndcapEta[0] && fabs(eta) < barrelEndcapEta[1])  {
            isLeptonFiducial=false;
            break;
          }
          if(fabs(eta)>eleEtaCut) {
            isLeptonFiducial=false;
            break;
          }
          if(eta>hem1516Eta[0] && eta<hem1516Eta[1] && phi>hem1516Phi[0] && phi<hem1516Phi[1]){
            isLeptonFiducial=false;
            break;
          }          
        }
      }

      is1bFiducial=(isLeptonFiducial && genBjets.size()>0);
      is2bFiducial=(isLeptonFiducial && genBjets.size()>1);
      
      //event weights and fiducial counters   
      if(isMC) {
        evWgt=fForestTree.ttbar_w->size()==0 ? 1. : fForestTree.ttbar_w->at(meIdxList[0]);
        if(allWgtSum.size()==0) allWgtSum.resize(meIdxList.size(),0.);
        for(size_t i=0; i<meIdxList.size(); i++) {
          Double_t iwgt(fForestTree.ttbar_w->size()<i  || fForestTree.ttbar_w->size() == 0 ? 1. : fForestTree.ttbar_w->at(meIdxList[i]));
          allWgtSum[i]+=iwgt;
          ht.fill2D("fidcounter",0,i,iwgt,"gen");
          if(isGenDilepton)    ht.fill2D("fidcounter",1,i,iwgt,"gen");
          if(isLeptonFiducial) ht.fill2D("fidcounter",2,i,iwgt,"gen");
          if(is1bFiducial)     ht.fill2D("fidcounter",3,i,iwgt,"gen");
          if(is2bFiducial)     ht.fill2D("fidcounter",4,i,iwgt,"gen");
        }
      }
    }
        
    wgtSum += evWgt;    
    float plotWgt(evWgt);
    
    //apply global filters
    if(!isMC && GT.find("103X")!=string::npos){
      if(TMath::Abs(fForestTree.vz) > 20) continue;
      if(!fForestSkim.phfCoincFilter2Th4) continue;
      if(!fForestSkim.pclusterCompatibilityFilter) continue;
      if(!fForestSkim.pprimaryVertexFilter) continue;
    }
    else if(!isMC && GT.find("75X")!=string::npos){
      if(TMath::Abs(fForestTree.vz) > 15) continue;
      if(!fForestSkim.phfCoincFilter) continue;
      if(!fForestSkim.HBHENoiseFilterResult) continue;
      if(!fForestSkim.pcollisionEventSelection) continue;
      if(!fForestSkim.pprimaryVertexFilter) continue;
    }

    //build jets from different PF candidate collections  
    std::cout << "checking PF cand\t" << fForestPF.nPF << std::endl; 
    SlimmedPFCollection_t pfColl;
    for(int ipf=0; ipf<fForestPF.nPF; ipf++) {
      int id(abs(fForestPF.pfId->at(ipf)));
      float mass(0.13957);  //pions
      if(id==4) mass=0.;    //photons
      if(id>=5) mass=0.497; //K0L
      pfColl.push_back( getSlimmedPF( id, fForestPF.pfPt->at(ipf),fForestPF.pfEta->at(ipf),fForestPF.pfPhi->at(ipf),mass) );
    }

    Float_t globalrho = getRho(pfColl,{1,2,3,4,5,6},-1.,5.);

    //monitor trigger and centrality
    float cenBin(0.),ncoll(1.);
    bool isCentralEvent(false);
    if(!isPP){
      isCentralEvent=(fForestTree.hiBin<30);
      cenBin=0.5*fForestTree.hiBin;
      ncoll=findNcoll(fForestTree.hiBin)*ncollWgtNorm;
    }
    if(!isMC){
      Int_t runBin=lumiTool.getRunBin(fForestTree.run);
      Float_t lumi=lumiTool.getLumi(fForestTree.run);
      if(lumi>0.){
        if(etrig>0) ht.fill("ratevsrun",runBin,1./lumi,"e");
        if(mtrig>0) ht.fill("ratevsrun",runBin,1./lumi,"m");
      }
    }

    //the selected leptons
    std::vector<LeptonSummary> selLeptons;
    
    //select muons
    std::vector<LeptonSummary> noIdMu;
    std::vector<TLorentzVector> muHLTP4;
    if(muHLTObjs) muHLTP4=muHLTObjs->getHLTObjectsP4();
    cout << __LINE__ << endl;
    for(unsigned int muIter = 0; muIter < fForestLep.muPt->size(); ++muIter) {
    cout << __LINE__ << endl;
      //kinematics selection
      TLorentzVector p4(0,0,0,0);
      float rawpt(fForestLep.muPt->at(muIter));
      p4.SetPtEtaPhiM(rawpt,fForestLep.muEta->at(muIter),fForestLep.muPhi->at(muIter),0.1057);

      //no specific calibration for muons:this is just left for symmetry with what is done for electrons
      float calpt=p4.Pt(); 
      p4 *= calpt/rawpt;  

      if(TMath::Abs(p4.Eta()) > muEtaCut) continue;
      if(p4.Pt() < lepPtCut) continue;

      bool isTrigMatch(false);
      for(auto hp4: muHLTP4) {
        if(hp4.DeltaR(p4)>0.1) continue;
        isTrigMatch=true;
        break;
      }
      
      LeptonSummary l(13,p4);
      l.rawpt = p4.Pt()*(rawpt/calpt); // calpt == pt for muons
      l.isTrigMatch=isTrigMatch;
      l.charge  = fForestLep.muCharge->at(muIter);
      l.chiso   = fForestLep.muPFChIso->at(muIter);
      l.nhiso   = fForestLep.muPFNeuIso->at(muIter);
      l.phoiso  = fForestLep.muPFPhoIso->at(muIter);
      l.isofull = l.chiso+l.nhiso+l.phoiso;
      // int   tmp_rhoind  = getRhoIndex(p4.Eta(),t_etaMin,t_etaMax);
      // l.rho = isPP ? globalrho : t_rho->at(tmp_rhoind);
      
      l.isofullR=getIsolationFull( pfColl, l.p4);
      l.miniiso = getMiniIsolation( pfColl ,l.p4, l.id);
      l.d0      = fForestLep.muD0   ->at(muIter);
      l.d0err   = 0.; //fForestLep.muD0Err->at(muIter); // no d0err for muons!!!
      l.dz      = fForestLep.muDz   ->at(muIter);
      l.origIdx = muIter;
      l.isMatched=false;
      l.isTauFeedDown=false;
      for(size_t ig=0;ig<genLeptons.size(); ig++) {
        if(genLeptons[ig].DeltaR(l.p4)<0.1) continue;
        l.isMatched=true;
        l.isTauFeedDown=genTauLeptons[ig];
      }
    
      noIdMu.push_back(l);

      //id (Tight muon requirements)
      int type=fForestLep.muType->at(muIter);
      bool isGlobal( ((type>>1)&0x1) );
      if(!isGlobal) continue;
      bool isPF( ((type>>5)&0x1) );
      if(!isPF) continue;
      bool isGlobalMuonPromptTight(fForestLep.muChi2NDF->at(muIter)<10. && fForestLep.muMuonHits->at(muIter)>0);
      if(!isGlobalMuonPromptTight) continue;
      if(fForestLep.muStations->at(muIter)<=1) continue;
      if(fForestLep.muTrkLayers->at(muIter) <= 5) continue;
      if(fForestLep.muPixelHits->at(muIter) == 0) continue;
      if(TMath::Abs(fForestLep.muInnerD0->at(muIter)) >=0.2 ) continue;
      if(TMath::Abs(fForestLep.muInnerDz->at(muIter)) >=0.5) continue;

      l.idFlags=1;

      //selected a good muon
      selLeptons.push_back(l);
    }
    std::sort(noIdMu.begin(),noIdMu.end(),orderByPt);
        
    //monitor muon id variables
    if(noIdMu.size()>1) {
      TLorentzVector p4[2] = {noIdMu[0].p4,      noIdMu[1].p4};
      int midx[2]          = {noIdMu[0].origIdx, noIdMu[1].origIdx};
      int charge(noIdMu[0].charge*noIdMu[1].charge);
 
      TString cat("zmmctrl");
      if(charge>0) cat="ss"+cat;
      float mmm((p4[0]+p4[1]).M());
      if(mmm<20) continue;

      ht.fill("mmll",  mmm,         plotWgt,cat);
      if( fabs(mmm-91)<15 && isSingleMuPD && mtrig>0 ) {
        for(size_t i=0; i<2; i++) {
          ht.fill("mmusta",    fForestLep.muStations->at(midx[i]),            plotWgt,cat);
          ht.fill("mtrklay",   fForestLep.muTrkLayers->at(midx[i]),           plotWgt,cat);
          ht.fill("mchi2ndf",  fForestLep.muChi2NDF->at(midx[i]),             plotWgt,cat);
          ht.fill("mmuhits",   fForestLep.muMuonHits->at(midx[i]),            plotWgt,cat);
          ht.fill("mpxhits",   fForestLep.muPixelHits->at(midx[i]),           plotWgt,cat);
          ht.fill("md0",       TMath::Abs(fForestLep.muInnerD0->at(midx[i])), plotWgt,cat);
          ht.fill("mdz",       TMath::Abs(fForestLep.muInnerDz->at(midx[i])), plotWgt,cat);          
        }
      }
    }


    //select electrons
    //cf. https://twiki.cern.ch/twiki/pub/CMS/HiHighPt2019/HIN_electrons2018_followUp.pdf
    std::vector<LeptonSummary> noIdEle;
    std::vector<TLorentzVector> eleHLTP4;
    if(eleHLTObjs) eleHLTP4=eleHLTObjs->getHLTObjectsP4() ;
    for(unsigned int eleIter = 0; eleIter < fForestLep.elePt->size(); ++eleIter) {

      //kinematics selection
      TLorentzVector p4(0,0,0,0);
      float rawpt(fForestLep.elePt->at(eleIter));
      p4.SetPtEtaPhiM(rawpt,fForestLep.eleEta->at(eleIter),fForestLep.elePhi->at(eleIter),0.000511);

      //deprecated
      //apply ad-hoc shift for endcap electrons if needed, i.e., PromptReco'18
      //if(!isMC && fForestTree.run<=firstEEScaleShiftRun && TMath::Abs(p4.Eta())>=barrelEndcapEta[1] && GT.find("fixEcalADCToGeV")==string::npos && GT.find("75X")==string::npos)
      //  p4 *=eeScaleShift;         
      float calpt=calibratedPt(rawpt, p4.Eta(), cenBin, isMC);
      p4 *= calpt/rawpt;
      if(TMath::Abs(p4.Eta()) > eleEtaCut) continue;
      if(TMath::Abs(p4.Eta()) > barrelEndcapEta[0] && TMath::Abs(p4.Eta()) < barrelEndcapEta[1] ) continue;
      if(p4.Pt() < lepPtCut) continue;
      bool isTrigMatch(false);
      for(auto hp4: eleHLTP4) {
        if(hp4.DeltaR(p4)>0.2) continue;
        isTrigMatch=true;
        break;
      }

      LeptonSummary l(11,p4);
      l.rawpt = p4.Pt()*(rawpt/calpt);
      l.isTrigMatch=isTrigMatch;
      l.charge  = fForestLep.eleCharge->at(eleIter);
      if(GT.find("75X_mcRun2")==string::npos) {
	l.chiso   = fForestLep.elePFChIso03->at(eleIter);
	l.nhiso   = fForestLep.elePFNeuIso03->at(eleIter);
	l.phoiso  = fForestLep.elePFPhoIso03->at(eleIter);
      }else{
        l.chiso   = fForestLep.elePFChIso->at(eleIter);
        l.nhiso   = fForestLep.elePFNeuIso->at(eleIter);
        l.phoiso  = fForestLep.elePFPhoIso->at(eleIter);
      }
      l.isofull = l.chiso+l.nhiso+l.phoiso;
      // int   tmp_rhoind  = getRhoIndex(p4.Eta(),t_etaMin,t_etaMax);
      // l.rho = isPP ? globalrho : t_rho->at(tmp_rhoind);

      l.isofullR= getIsolationFull( pfColl, l.p4);
      l.miniiso = getMiniIsolation( pfColl ,l.p4, l.id);
      l.d0      = fForestLep.eleD0   ->at(eleIter);
      l.d0err   = fForestLep.eleD0Err->at(eleIter);
      l.dz      = fForestLep.eleDz   ->at(eleIter);
      l.origIdx=eleIter;
      l.isMatched=false;
      l.isTauFeedDown=false;
      for(size_t ig=0;ig<genLeptons.size(); ig++) {
        if(genLeptons[ig].DeltaR(l.p4)<0.1) continue;
        l.isMatched=true;
        l.isTauFeedDown=genTauLeptons[ig];
      }

      noIdEle.push_back(l);
      
      l.idFlags=getElectronId(TMath::Abs(fForestLep.eleSCEta->at(eleIter))< barrelEndcapEta[0],
                                 fForestLep.eleSigmaIEtaIEta->at(eleIter),
                                 fForestLep.eledEtaSeedAtVtx->at(eleIter),
                                 fForestLep.eledPhiAtVtx->at(eleIter),
                                 fForestLep.eleHoverEBc->at(eleIter),
                                 fForestLep.eleEoverPInv->at(eleIter),
                                 fForestLep.eleIP3D->at(eleIter),
                                 fForestLep.eleMissHits->at(eleIter),
                                 isCentralEvent);
      
      //id'ed electron
      if(!isLooseElectron(l.idFlags)) continue;
      selLeptons.push_back(l);
    }
    std::sort(noIdEle.begin(),noIdEle.end(),orderByPt);       

    //monitor electron id variables
    if(noIdEle.size()>1) {
      TLorentzVector p4[2] = {noIdEle[0].p4,noIdEle[1].p4};
      int eidx[2]          = {noIdEle[0].origIdx,noIdEle[1].origIdx};
      int charge(noIdEle[0].charge*noIdEle[1].charge);
      float mee((p4[0]+p4[1]).M());
      if(mee<20) continue;
      
      TString basecat("zeectrl");
      if(charge>0) basecat="ss"+basecat;      
      TString cat(basecat);
      if(fabs(noIdEle[0].p4.Eta())>=barrelEndcapEta[1] && fabs(noIdEle[1].p4.Eta())>=barrelEndcapEta[1])
        cat +="EE";
      else if(fabs(noIdEle[0].p4.Eta())>=barrelEndcapEta[1] || fabs(noIdEle[1].p4.Eta())>=barrelEndcapEta[1])
        cat+="EB";
      else
        cat+="BB";
      ht.fill("emll",  mee,         plotWgt,cat);

      if( fabs(mee-91)<15 && isSingleElePD && etrig>0) {
        for(size_t i=0; i<2; i++) {
          cat=basecat;
          cat += (fabs(p4[i].Eta())>=barrelEndcapEta[1] ? "EE" : "EB");
          ht.fill("esihih",  fForestLep.eleSigmaIEtaIEta->at(eidx[i]),         plotWgt,cat);
          ht.fill("edetaseedvtx", TMath::Abs(fForestLep.eledEtaSeedAtVtx->at(eidx[i])), plotWgt,cat);
          ht.fill("edphivtx", TMath::Abs(fForestLep.eledPhiAtVtx->at(eidx[i])), plotWgt,cat);
          ht.fill("ehoebc",     fForestLep.eleHoverEBc->at(eidx[i]),                plotWgt,cat);
          ht.fill("eempinv",  fForestLep.eleEoverPInv->at(eidx[i]),             plotWgt,cat);
          ht.fill("e3dip",    TMath::Abs(fForestLep.eleIP3D->at(eidx[i])),        plotWgt,cat);
        }
      }
    }

    //sort selected electrons by pt
    std::sort(selLeptons.begin(),selLeptons.end(),orderByPt);

    //monitor trigger efficiency
    if(selLeptons.size()>=2){
      for(size_t ilep=0; ilep<2; ilep++){
        if(!selLeptons[ilep].isMatched) continue;
        TString cat( abs(selLeptons[ilep].id)==11 ? "e" : "m");
        float pt(selLeptons[ilep].p4.Pt()), abseta(fabs(selLeptons[ilep].p4.Eta()));
        ht.fill("trig_pt",  pt,     ncoll, cat);          
        ht.fill("trig_eta", abseta, ncoll, cat);          
        if(!selLeptons[ilep].isTrigMatch) continue;
        cat+="match";
        ht.fill("trig_pt",  pt,     ncoll, cat);          
        ht.fill("trig_eta", abseta, ncoll, cat);          
      }
    }

    //require at least two leptons matched to trigger objects
    if(selLeptons.size()<2) continue;
    bool hasOneTrigMatchLepton(selLeptons[0].isTrigMatch || selLeptons[1].isTrigMatch);
    if( !hasOneTrigMatchLepton ) continue;
    
    //apply trigger preselection & duplicate event removal
    //in skim mode assume that the hltobject matched to offline will give the HLT trigger bit
    //this is a hack when hlt tree is missing...
    if(isSkim) {
      if( (abs(selLeptons[0].id)==11 && selLeptons[0].isTrigMatch) ||
          (abs(selLeptons[1].id)==11 && selLeptons[1].isTrigMatch) ) etrig=true;
      if( (abs(selLeptons[0].id)==13 && selLeptons[0].isTrigMatch) ||
          (abs(selLeptons[1].id)==13 && selLeptons[1].isTrigMatch) ) mtrig=true;
    }

    int trig=etrig+mtrig;
    if(trig==0) continue;

    if(isSingleMuPD || isMuSkimedMCPD) {
      if(std::find(badMuonTriggerRuns.begin(), badMuonTriggerRuns.end(), fForestTree.run) != badMuonTriggerRuns.end() and !isMuSkimedMCPD) continue;
      if(mtrig==0) continue;
      if(etrig!=0) continue;
      if(selLeptons.size()>=2)
	if ( (abs(selLeptons[0].id)==11 and selLeptons[0].isTrigMatch==1) and (abs(selLeptons[1].id)==11 and selLeptons[1].isTrigMatch==1) ) continue;
    }
    if(isSingleElePD || isEleSkimedMCPD) {
      if(etrig==0) continue;
      if(selLeptons.size()>=2)
	if ( (abs(selLeptons[0].id)==13 and selLeptons[0].isTrigMatch==1) or (abs(selLeptons[1].id)==13 and selLeptons[1].isTrigMatch==1) ) continue;
    }

    //dilepton selection
    TLorentzVector ll(selLeptons[0].p4+selLeptons[1].p4);
    TLorentzVector ll_raw(selLeptons[0].p4*(selLeptons[0].rawpt/selLeptons[0].p4.Pt())+selLeptons[1].p4*(selLeptons[1].rawpt/selLeptons[1].p4.Pt()));
    t_llpt=ll.Pt();
    t_llpt_raw=ll_raw.Pt();
    t_lleta=ll.Eta();
    t_llphi=ll.Phi();
    t_llm=ll.M();
    t_llm_raw=ll_raw.M();
    t_dphi=TMath::Abs(selLeptons[0].p4.DeltaPhi(selLeptons[1].p4));
    t_deta=fabs(selLeptons[0].p4.Eta()-selLeptons[1].p4.Eta());
    t_sumeta=selLeptons[0].p4.Eta()+selLeptons[1].p4.Eta();
    int dilCode(selLeptons[0].id*selLeptons[1].id);
    TString dilCat("mm");
    if(dilCode==11*13) dilCat="em";
    if(dilCode==11*11) dilCat="ee";

    //ee and mm events should come from the appropriate primary dataset
    if(isSingleMuPD || isSingleMuPD || isMuSkimedMCPD || isEleSkimedMCPD) {
      if(dilCode==11*11 && !isSingleElePD && !isEleSkimedMCPD) continue;
      if(dilCode==13*13 && !isSingleMuPD && !isMuSkimedMCPD) continue;
    }

    if(blind) {
      bool isZ( dilCode!=11*13 && fabs(t_llm-91)<15);
      int charge(selLeptons[0].charge*selLeptons[1].charge);
      if(!isMC && !isZ && charge<0 && fForestTree.run>=326887) continue;
    }      
              
    //analyze jets
    std::vector<BtagInfo_t> pfJetsIdx;
    std::vector<TLorentzVector> pfJetsP4;
    int npfjets(0),npfbjets(0); 

    // initialize all the counters
    t_nbjet_sel = 0; 
    t_nbjet_sel_jecup    = 0; t_nbjet_sel_jecdn    = 0;
    t_nbjet_sel_jerup    = 0; t_nbjet_sel_jerdn    = 0;
    t_nbjet_sel_bup      = 0; t_nbjet_sel_bdn      = 0;
    t_nbjet_sel_udsgup   = 0; t_nbjet_sel_udsgdn   = 0;
    t_nbjet_sel_quenchup = 0; t_nbjet_sel_quenchdn = 0;

	
    for(int jetIter = 0; jetIter < fForestJets.nref; jetIter++){

      //at least two tracks
      if(fForestJets.trackN[jetIter]<2) continue;

      TLorentzVector jp4(0,0,0,0);
      if(isMC){
	JECMC.SetJetPT(fForestJets.rawpt[jetIter]);
	JECMC.SetJetEta(fForestJets.jteta[jetIter]);
	JECMC.SetJetPhi(fForestJets.jtphi[jetIter]);
	jp4.SetPtEtaPhiM( JECMC.GetCorrectedPT(),fForestJets.jteta[jetIter],fForestJets.jtphi[jetIter],fForestJets.jtm[jetIter]);
	JEUMC.SetJetPT(JECMC.GetCorrectedPT());
	JEUMC.SetJetEta(fForestJets.jteta[jetIter]);
	JEUMC.SetJetPhi(fForestJets.jtphi[jetIter]);
      }
      else{
	JECData.SetJetPT(fForestJets.rawpt[jetIter]);
	JECData.SetJetEta(fForestJets.jteta[jetIter]);
	JECData.SetJetPhi(fForestJets.jtphi[jetIter]);
	jp4.SetPtEtaPhiM( JECData.GetCorrectedPT(),fForestJets.jteta[jetIter],fForestJets.jtphi[jetIter],fForestJets.jtm[jetIter]);
	JEUData.SetJetPT(JECMC.GetCorrectedPT());
	JEUData.SetJetEta(fForestJets.jteta[jetIter]);
	JEUData.SetJetPhi(fForestJets.jtphi[jetIter]);
      }

      float csvVal=fForestJets.discr_csvV2[jetIter];
      int nsvtxTk=fForestJets.svtxntrk[jetIter];
      float msvtx=fForestJets.svtxm[jetIter];

      if(jp4.Pt()<20.) continue; // smaller pT cut here to avoid the full loop
      if(fabs(jp4.Eta())>2.0) continue;
      bool isBTagged(csvVal>csvWPList[csvWP]);      

      // simple matching to the closest jet in dR. require at least dR < 0.3
      TLorentzVector matchjp4(0,0,0,0);
      int refFlavor(0),refFlavorForB(0);
      if (isMC){
        //std::vector<TLorentzVector> matchedJets;
        for (int genjetIter = 0; genjetIter < fForestJets.ngen; genjetIter++){
          if (jetIter == fForestJets.genmatchindex[genjetIter]) {
            matchjp4.SetPtEtaPhiM( fForestJets.genpt[genjetIter],fForestJets.geneta[genjetIter],fForestJets.genphi[genjetIter],fForestJets.genm[genjetIter]);
          }
        }  
        refFlavor=fForestJets.refparton_flavor[jetIter];
        refFlavorForB=fForestJets.refparton_flavorForB[jetIter];
      }

      //cross clean wrt to leptons
      if(jp4.DeltaR(selLeptons[0].p4)<0.4 || jp4.DeltaR(selLeptons[1].p4)<0.4) continue;
      
      pfJetsIdx.push_back(std::make_tuple(pfJetsP4.size(),nsvtxTk,msvtx,csvVal,matchjp4,refFlavor,refFlavorForB));
      pfJetsP4.push_back(jp4);
      npfjets++;
      npfbjets += isBTagged;
      
      
      if (jp4.Pt() > 30. && isBTagged) t_nbjet_sel       += 1;

      if (isMC){

        if (jp4.Pt() > 30 * (1 + JEUMC.GetUncertainty().first) && isBTagged) t_nbjet_sel_jecup += 1;
        if (jp4.Pt() > 30 * (1 - JEUMC.GetUncertainty().second) && isBTagged) t_nbjet_sel_jecdn += 1;
        float cjer(0.);
        if ( abs(refFlavorForB) ) cjer = 1. + (1.2 -1.) * (jp4.Pt() - matchjp4.Pt()) / jp4.Pt(); // hard coded 1.2
        else cjer = rand->Gaus(1., 0.2);

        if (jp4.Pt()*cjer > 30. && isBTagged) t_nbjet_sel_jerup += 1;
        if (jp4.Pt()/cjer > 30. && isBTagged) t_nbjet_sel_jerdn += 1;

        quenchingModel->SetParameter(0, 50.); // this sets the omega_c parameter. if we want to make this centrality dependent
        float tmp_quench_loss = quenchingModel->GetRandom();
        tmp_quench_loss = TMath::Abs(TMath::Sin(jp4.Theta())*tmp_quench_loss); // make it only on the transverse part...
        // make it centrality dependent
        float centralitySuppression = centralityModel->Eval(cenBin/100.);
        //
        // std::cout << "energy loss " << tmp_quench_loss << " due to quenching " << std::endl;
        // std::cout << "jet has eta " << jp4.Eta() << std::endl;
        // std::cout << "pT loss     " << TMath::Abs(TMath::Sin(jp4.Theta())*tmp_quench_loss) << " due to quenching " << std::endl;
        // std::cout << "pT loss centrality dependen    " << TMath::Abs(TMath::Sin(jp4.Theta())*tmp_quench_loss)*centralitySuppression << " (cen/supp) " << cenBin << " / " << centralitySuppression << std::endl << std::endl;

        if (jp4.Pt()                  > 30. && isBTagged) t_nbjet_sel_quenchup += 1;
        if (jp4.Pt()-tmp_quench_loss*centralitySuppression  > 30. && isBTagged) t_nbjet_sel_quenchdn += 1;
        ht.fill("jetptprequench" ,  jp4.Pt()                                      ,  plotWgt);
        ht.fill("jetquenchloss"  ,  tmp_quench_loss*centralitySuppression         ,  plotWgt);
        if (jp4.Pt()-tmp_quench_loss*centralitySuppression > 20.)  ht.fill("jetptpostquench",  jp4.Pt()-tmp_quench_loss*centralitySuppression,  plotWgt);

        if (abs(refFlavorForB) == 5){
            ht.fill("jetptprequenchB" ,  jp4.Pt()                                      ,  plotWgt);
            ht.fill("jetquenchlossB"  ,  tmp_quench_loss*centralitySuppression         ,  plotWgt);
            if (jp4.Pt()-tmp_quench_loss*centralitySuppression > 20.)  ht.fill("jetptpostquenchB",  jp4.Pt()-tmp_quench_loss*centralitySuppression,  plotWgt);
        }

        bool isBTaggedNew(0);
        float tmp_btageff = btagEfficiencies(refFlavorForB, cenBin);

        if (jp4.Pt() < 30.) continue;
        t_nbjet_sel_bup    = t_nbjet_sel; t_nbjet_sel_bdn    = t_nbjet_sel;
        t_nbjet_sel_udsgup = t_nbjet_sel; t_nbjet_sel_udsgdn = t_nbjet_sel;

        if ( abs(refFlavorForB) == 5){

          isBTaggedNew = isBTagged; 
          myBTagUtil->modifyBTagsWithSF(isBTaggedNew, 1.05, tmp_btageff );
          if (!isBTagged &&  isBTaggedNew) t_nbjet_sel_bup = t_nbjet_sel+1;
          if ( isBTagged && !isBTaggedNew) t_nbjet_sel_bup = t_nbjet_sel-1;

          isBTaggedNew = isBTagged; 
          myBTagUtil->modifyBTagsWithSF(isBTaggedNew, 0.95, tmp_btageff );
          if (!isBTagged &&  isBTaggedNew) t_nbjet_sel_bdn = t_nbjet_sel+1;
          if ( isBTagged && !isBTaggedNew) t_nbjet_sel_bdn = t_nbjet_sel-1;
        }

        else {

          isBTaggedNew = isBTagged; 
          myBTagUtil->modifyBTagsWithSF(isBTaggedNew, 1.15, tmp_btageff );
          if (!isBTagged &&  isBTaggedNew) t_nbjet_sel_udsgup = t_nbjet_sel+1;
          if ( isBTagged && !isBTaggedNew) t_nbjet_sel_udsgup = t_nbjet_sel-1;

          isBTaggedNew = isBTagged; 
          myBTagUtil->modifyBTagsWithSF(isBTaggedNew, 0.85, tmp_btageff );
          if (!isBTagged &&  isBTaggedNew) t_nbjet_sel_udsgdn = t_nbjet_sel+1;
          if ( isBTagged && !isBTaggedNew) t_nbjet_sel_udsgdn = t_nbjet_sel-1;

        }
      }
    }

    // set the bjet variations to the nominal one for data
    if (!isMC){
      t_nbjet_sel_jecup    = t_nbjet_sel; t_nbjet_sel_jecdn    = t_nbjet_sel;
      t_nbjet_sel_jerup    = t_nbjet_sel; t_nbjet_sel_jerdn    = t_nbjet_sel;
      t_nbjet_sel_bup      = t_nbjet_sel; t_nbjet_sel_bdn      = t_nbjet_sel;
      t_nbjet_sel_udsgup   = t_nbjet_sel; t_nbjet_sel_udsgdn   = t_nbjet_sel;
      t_nbjet_sel_quenchup = t_nbjet_sel; t_nbjet_sel_quenchdn = t_nbjet_sel;
    }

    std::sort(pfJetsIdx.begin(),       pfJetsIdx.end(),      orderByBtagInfo);

    //for gen fill again fiducial counters
    if(isMC) {      
      
      bool isMatchedDilepton(abs(genDileptonCat)==abs(dilCode));
      if( (genDileptonCat==11*11 && etrig==0) || (genDileptonCat==13*13 && mtrig==0)) 
        isMatchedDilepton=false;

      std::vector<TString> fidCats;
      fidCats.push_back( isMatchedDilepton   ? "lep"    : "fakelep" );
      if(npfbjets>0) {
        fidCats.push_back( isMatchedDilepton && is1bFiducial ? "lep1b" : "fakelep1b" );
        if(npfbjets>1) {
          fidCats.push_back( isMatchedDilepton && is2bFiducial ? "lep2b" : "fakelep2b" );
        }
      }    
      
      for(size_t i=0; i<meIdxList.size(); i++) {
        Double_t iwgt(fForestTree.ttbar_w->size()<i || fForestTree.ttbar_w->size() == 0 ? 1. : fForestTree.ttbar_w->at(meIdxList[i]));
        ht.fill2D("fidcounter",0,i,iwgt,fidCats);
        if(isGenDilepton)    ht.fill2D("fidcounter",1,i,iwgt,fidCats);
        if(isLeptonFiducial) ht.fill2D("fidcounter",2,i,iwgt,fidCats);
        if(is1bFiducial)     ht.fill2D("fidcounter",3,i,iwgt,fidCats);
        if(is2bFiducial)     ht.fill2D("fidcounter",4,i,iwgt,fidCats);
      }
    }


    //define categories for pre-selection control histograms
    std::vector<TString> categs;
    categs.push_back(dilCat);
    
    std::vector<TString> addCategs;
    
    //monitor after run where EE scale shift changed
    if(!isPP){
      addCategs.clear();
      TString pf( fForestTree.run>=firstEEScaleShiftRun ? "after" : "before" );
      for(auto c : categs) {
        addCategs.push_back(c); addCategs.push_back(c+pf); 
      }
      categs=addCategs;
    }

    //monitor according to the b-tagging category
    addCategs.clear();
    TString pfbcat(Form("%dpfb",min(npfbjets,2)));
    for(auto c : categs) { 
      addCategs.push_back(c); 
      if(npfbjets==0) addCategs.push_back(c+"0pfb"); 
      if(npfbjets>0) addCategs.push_back(c+"geq1pfb"); 
    }
    categs=addCategs;

        
    //fill histograms
    for(int i=0; i<2; i++) {
      TString pf(Form("l%d",i+1));
      float pt(selLeptons[i].p4.Pt());
      ht.fill(pf+"pt",             pt,                            plotWgt, categs);
      ht.fill(pf+"eta",            fabs(selLeptons[i].p4.Eta()),  plotWgt, categs);
      ht.fill(pf+"chiso",          selLeptons[i].chiso,           plotWgt, categs);
      ht.fill(pf+"phoiso",         selLeptons[i].phoiso,          plotWgt, categs);
      ht.fill(pf+"nhiso",          selLeptons[i].nhiso,           plotWgt, categs);

    }

    ht.fill( "acopl",     1-fabs(t_dphi)/TMath::Pi(),                   plotWgt, categs);
    ht.fill( "detall",    t_deta,                                       plotWgt, categs);
    ht.fill( "drll",      selLeptons[0].p4.DeltaR(selLeptons[1].p4),    plotWgt, categs);
    ht.fill( "mll",       t_llm,                                        plotWgt, categs);
    ht.fill( "ptll",      t_llpt,                                       plotWgt, categs);
    ht.fill( "ptsum",     selLeptons[0].p4.Pt()+selLeptons[1].p4.Pt(),  plotWgt, categs);

    //PF jets
    ht.fill( "npfjets",   npfjets,   plotWgt, categs);
    ht.fill( "npfbjets",  npfbjets,  plotWgt, categs);
    std::vector<TLorentzVector> pfFinalState;
    pfFinalState.push_back(selLeptons[0].p4);
    pfFinalState.push_back(selLeptons[1].p4);
    for(size_t ij=0; ij<min(pfJetsIdx.size(),size_t(2)); ij++) {     
      int idx(std::get<0>(pfJetsIdx[ij]));
      int ntks(std::get<1>(pfJetsIdx[ij]));
      float svm(std::get<2>(pfJetsIdx[ij]));
      float csv(std::get<3>(pfJetsIdx[ij]));
      TLorentzVector p4=pfJetsP4[idx];
      if(csv>csvWPList[csvWP]) pfFinalState.push_back(p4);
      TString ppf(ij==0 ? "1" : "2");
      ht.fill( "pf"+ppf+"jbalance", p4.Pt()/ll.Pt(), plotWgt, categs);
      ht.fill( "pf"+ppf+"jpt",      p4.Pt(),         plotWgt, categs);
      ht.fill( "pf"+ppf+"jeta",     fabs(p4.Eta()),  plotWgt, categs);
      ht.fill( "pf"+ppf+"jsvtxm",   ntks,            plotWgt, categs);
      ht.fill( "pf"+ppf+"jsvtxntk", svm,             plotWgt, categs);
      ht.fill( "pf"+ppf+"jcsv",     csv,             plotWgt, categs);
      ht.fill2D( "pf"+ppf+"jetavsphi",   p4.Eta(),p4.Phi(),   plotWgt, categs);


      quenchingModel->SetParameter(0, 50.); // this sets the omega_c parameter. if we want to make this centrality dependent
      float tmp_quench_loss = quenchingModel->GetRandom();
      tmp_quench_loss = TMath::Abs(TMath::Sin(p4.Theta())*tmp_quench_loss); // make it only on the transverse part...
      // make it centrality dependent
      float centralitySuppression = centralityModel->Eval(cenBin/100.);
      float quenchedPt=p4.Pt()-tmp_quench_loss*centralitySuppression;
      ht.fill2D("jptvsjptquench", p4.Pt(), quenchedPt, plotWgt,categs);
    }



    
    std::vector<float> rapMoments=getRapidityMoments(pfFinalState);
    ht.fill( "pfrapavg",     rapMoments[0], plotWgt, categs);
    ht.fill( "pfraprms",     rapMoments[1], plotWgt, categs);
    ht.fill( "pfrapmaxspan", rapMoments[2], plotWgt, categs);
    float pfht(0.);
    TLorentzVector vis(0,0,0,0);
    for(auto p : pfFinalState) { vis+=p; pfht+=p.Pt(); }
    ht.fill( "pfht",         pfht, plotWgt, categs);
    ht.fill( "pfmht",        vis.Pt(), plotWgt, categs);

    // for tree filling set all the proper variables
    t_run    = fForestTree.run;
    t_lumi   = fForestTree.lumi;
    t_event  = fForestTree.evt;
    t_vx     = fForestTree.vx;
    t_vy     = fForestTree.vy;
    t_vz     = fForestTree.vz;
    t_weight = plotWgt;
    t_meWeights.clear();
    if(isMC){
      t_meWeights.push_back(topPtWgt);
      t_meWeights.push_back(1./topPtWgt);
      t_meWeights.push_back(topMassUpWgt);
      t_meWeights.push_back(topMassDnWgt);
      if(fForestTree.ttbar_w->size()>0) {
        float nomWgt=fForestTree.ttbar_w->at(meIdxList[0]);
        for(auto meIdx : meIdxList){
          if(meIdx< fForestTree.ttbar_w->size()){
            t_meWeights.push_back( fForestTree.ttbar_w->at(meIdx)/nomWgt );
          }
        }
      }
      else {t_meWeights.push_back(1.0); } 
    }
        
    //centrality
    t_cenbin   = cenBin;
    t_ncollWgt = ncoll;

    t_globalrho = globalrho;
    t_etrig  = etrig;
    t_mtrig  = mtrig;

    //trigger scale factor
    t_trigSF    = 1.0;
    t_trigSFUnc = 0.;

    //get expected trigger efficiencies and measured scale factors
    std::vector<std::pair<float,float>  > ltrigEff, ltrigSF;
    for(size_t ilep=0; ilep<2; ilep++){
      float pt(selLeptons[ilep].p4.Pt()),eta(selLeptons[ilep].p4.Eta()),abseta(fabs(eta));

      if(abs(selLeptons[ilep].id)==11){
        ltrigEff.push_back(  std::pair<float,float>(e_mctrigeff->Eval(pt),0.0) );
        ltrigSF.push_back( eleEff.eval(pt, abseta<barrelEndcapEta[0], cenBin, true, false) ); //HLT (L1 is unity by definition in this trigger menu)
      }else{

        ltrigEff.push_back(  std::pair<float,float>(tnp_weight_trig_pbpb(pt,eta,cenBin),0.0) );
        float mutrigSF=tnp_weight_trig_pbpb(pt,eta,cenBin,0);
        float deltaTnp=max(fabs(mutrigSF-tnp_weight_trig_pbpb(pt,eta,cenBin,-1)),fabs(mutrigSF-tnp_weight_trig_pbpb(pt,eta,cenBin,-2)));
        float deltaStat=max(fabs(mutrigSF-tnp_weight_trig_pbpb(pt,eta,cenBin,1)),fabs(mutrigSF-tnp_weight_trig_pbpb(pt,eta,cenBin,2)));
        float deltaSF=sqrt(deltaTnp*deltaTnp+deltaStat*deltaStat);
        ltrigSF.push_back(  std::pair<float,float>(mutrigSF,deltaSF)  );        
      }
    }

    //trigeff = e1*e2 +e1*(1-e2)+e2*(1-e1), the rest is scale factor and error propagation
    t_trigSF  = (ltrigSF[0].first*ltrigEff[0].first+ltrigSF[1].first*ltrigEff[1].first-ltrigSF[0].first*ltrigSF[1].first*ltrigEff[0].first*ltrigEff[1].first);
    t_trigSF /= (           
      ltrigEff[0].first+                 ltrigEff[1].first-                                  ltrigEff[0].first*ltrigEff[1].first);

    t_trigSFUnc  = pow( ltrigSF[0].second*(ltrigEff[0].first-ltrigSF[1].first*ltrigEff[0].first*ltrigEff[1].first), 2 );
    t_trigSFUnc += pow( ltrigSF[1].second*(ltrigEff[1].first-ltrigSF[0].first*ltrigEff[0].first*ltrigEff[1].first), 2 );
    t_trigSFUnc  = sqrt(t_trigSFUnc);

    // fill the leptons ordered by pt
    t_lep_pt    .clear();
    t_lep_calpt .clear();
    t_lep_eta   .clear();
    t_lep_phi   .clear();
    t_lep_pdgId .clear();
    t_lep_idflags.clear();
    t_lep_d0 .clear();
    t_lep_d0err .clear();
    t_lep_dz  .clear();
    t_lep_charge.clear();
    t_lep_chiso.clear();    
    t_lep_phiso.clear();
    t_lep_nhiso.clear();
    t_lep_rho.clear();    
    t_lep_isofull.clear();
    t_lep_isofull20.clear();
    t_lep_isofull25.clear();
    t_lep_isofull30.clear();
    t_lep_miniiso.clear();
    t_lep_matched.clear();
    t_lep_taufeeddown.clear();
    t_lep_trigmatch.clear();
    t_lepSF.clear();
    t_lepSFUnc.clear();
    t_lepIsoSF.clear();
    t_lepIsoSFUnc.clear();
    t_nlep = selLeptons.size();
    t_lep_ind1 = -1;
    t_lep_ind2 = -1;
    for (int ilep = 0; ilep < t_nlep; ++ilep){
      t_lep_pt    .push_back( selLeptons[ilep].rawpt  );
      t_lep_calpt .push_back( selLeptons[ilep].p4.Pt()  );
      t_lep_eta   .push_back( selLeptons[ilep].p4.Eta() );
      t_lep_phi   .push_back( selLeptons[ilep].p4.Phi() );
      t_lep_idflags.push_back(selLeptons[ilep].idFlags);
      t_lep_d0    .push_back( selLeptons[ilep].d0 );
      t_lep_d0err .push_back( selLeptons[ilep].d0err);
      t_lep_dz    .push_back( selLeptons[ilep].dz  );
      t_lep_chiso .push_back( selLeptons[ilep].chiso );
      t_lep_phiso .push_back( selLeptons[ilep].phoiso );
      t_lep_nhiso .push_back( selLeptons[ilep].nhiso );
      t_lep_rho   .push_back( selLeptons[ilep].rho );
      t_lep_pdgId .push_back( selLeptons[ilep].id );
      t_lep_charge.push_back( selLeptons[ilep].charge );
      t_lep_isofull.push_back( selLeptons[ilep].isofull );
      t_lep_isofull20.push_back( selLeptons[ilep].isofullR[0] );
      t_lep_isofull25.push_back( selLeptons[ilep].isofullR[1] );
      t_lep_isofull30.push_back( selLeptons[ilep].isofullR[2] );
      t_lep_miniiso.push_back( selLeptons[ilep].miniiso );
      t_lep_matched.push_back( selLeptons[ilep].isMatched );
      t_lep_taufeeddown.push_back( selLeptons[ilep].isTauFeedDown );
      t_lep_trigmatch.push_back( selLeptons[ilep].isTrigMatch );
      
      //reco/tracking+id scale factors
      float sfVal(1.0),sfValUnc(0.0);
      if(abs(selLeptons[ilep].id)==13) {
	//ID
        sfVal=tnp_weight_muid_pbpb( selLeptons[ilep].p4.Eta(), 0 );                         //central value
        sfValUnc += pow(fabs(tnp_weight_muid_pbpb( selLeptons[ilep].p4.Eta(),+1)-sfVal),2); //stat, +1 sigma
        sfValUnc += pow(fabs(tnp_weight_muid_pbpb( selLeptons[ilep].p4.Eta(),-1)-sfVal),2); //stat, -1 sigma
	sfValUnc += pow(fabs(tnp_weight_muid_pbpb( selLeptons[ilep].p4.Eta(),+2)-sfVal),2); //syst, +1 sigma                                                                                         
	sfValUnc += pow(fabs(tnp_weight_muid_pbpb( selLeptons[ilep].p4.Eta(),-2)-sfVal),2); //syst, -1 sigma
	sfValUnc += pow(0.01,2);                                                            //identification centrality dependence
	//Tracking
	sfVal*=tnp_weight_glbtrk_pbpb( selLeptons[ilep].p4.Eta(), cenBin,0 );                      //central value                                                                                
	sfValUnc += pow(fabs(tnp_weight_glbtrk_pbpb( selLeptons[ilep].p4.Eta(),cenBin,+1)-sfVal),2); //stat, +1 sigma                                                                                 
	sfValUnc += pow(fabs(tnp_weight_glbtrk_pbpb( selLeptons[ilep].p4.Eta(),cenBin,-1)-sfVal),2); //stat, -1 sigma                                                                                    
	sfValUnc += pow(fabs(tnp_weight_glbtrk_pbpb( selLeptons[ilep].p4.Eta(),cenBin,+2)-sfVal),2); //syst, +1 sigma                                                                                    
	sfValUnc += pow(fabs(tnp_weight_glbtrk_pbpb( selLeptons[ilep].p4.Eta(),cenBin,-2)-sfVal),2); //syst, -1 sigma                                                                                    
	sfValUnc = sqrt(sfValUnc);
      }else {
        std::pair<float,float > eleIDsf=eleEff.eval(selLeptons[ilep].p4.Pt(), fabs(selLeptons[ilep].p4.Eta())<barrelEndcapEta[0], cenBin, false, false); //ID
        sfVal=eleIDsf.first;
        sfValUnc=eleIDsf.second;
        std::pair<float,float > eleRECOsf=eleEff.eval(selLeptons[ilep].p4.Pt(), fabs(selLeptons[ilep].p4.Eta())<barrelEndcapEta[0], cenBin, false, true); //RECO                                     
        sfValUnc = sqrt(pow(sfValUnc/sfVal,2)+pow(eleRECOsf.second/eleRECOsf.first,2));
        sfVal*=eleRECOsf.first;
      }    
      t_lepSF.push_back(sfVal);
      t_lepSFUnc.push_back(sfValUnc);

      TString isoKey(cenBin<30 ? "cen" : "periph");
      isoKey+=abs(selLeptons[ilep].id)==13 ? "_169" : "_121";
      Int_t xbin=isoEffSFs[isoKey]->GetXaxis()->FindBin( min(selLeptons[ilep].p4.Pt(), isoEffSFs[isoKey]->GetXaxis()->GetXmax()-0.01) );
      Int_t ybin=isoEffSFs[isoKey]->GetYaxis()->FindBin( min(fabs(selLeptons[ilep].p4.Eta()), isoEffSFs[isoKey]->GetYaxis()->GetXmax()) );
      if (isoEffSFs[isoKey]->GetBinContent(xbin,ybin) != 0.){
        t_lepIsoSF   .push_back( isoEffSFs[isoKey]->GetBinContent(xbin,ybin) );
        t_lepIsoSFUnc.push_back( isoEffSFs[isoKey]->GetBinError  (xbin,ybin) );
      } else{
        t_lepIsoSF   .push_back( 1.00);
        t_lepIsoSFUnc.push_back( 0.05);
      }

      //isolation-based indices
      bool isIso(true);
      if(abs(selLeptons[ilep].id)==13) {
        float rho=selLeptons[ilep].rho;
        float ue=0.00102*pow(rho+12.6255,2)+0.18535*(rho+12.6255);
        float iso=(selLeptons[ilep].isofull-ue)/selLeptons[ilep].p4.Pt();
        if(iso>0.12) isIso=false;
      }else {
        float rho=selLeptons[ilep].rho;
        float ue=0.000817*pow(rho+14.696,2)+0.201661*(rho+14.696);
        float iso=(selLeptons[ilep].isofull-ue)/selLeptons[ilep].p4.Pt();
        if(iso>0.) isIso=false;
      }
      if(isIso && t_lep_ind1 < 0)                    t_lep_ind1 = ilep;
      if(isIso && t_lep_ind1 > -1 && t_lep_ind2 < 0) t_lep_ind2 = ilep;
    }
    
    // fill the jets ordered by b-tag
    t_bjet_leadPassTight=false;
    t_bjet_pt   .clear();
    t_bjet_eta  .clear();
    t_bjet_phi  .clear();
    t_bjet_mass .clear();
    t_bjet_csvv2.clear();
    t_bjet_matchpt  .clear();
    t_bjet_matcheta .clear();
    t_bjet_matchphi .clear();
    t_bjet_matchmass.clear();
    t_bjet_flavor.clear();
    t_bjet_flavorForB.clear();
    t_nbjet = pfJetsIdx.size();
    for (int ij = 0; ij < t_nbjet; ij++) {
      int idx = std::get<0>(pfJetsIdx[ij]);
      t_bjet_pt   .push_back( pfJetsP4[idx].Pt()  );
      t_bjet_eta  .push_back( pfJetsP4[idx].Eta() );
      t_bjet_phi  .push_back( pfJetsP4[idx].Phi() );
      t_bjet_mass .push_back( pfJetsP4[idx].M()   );
      t_bjet_csvv2.push_back( std::get<3>(pfJetsIdx[ij])   );      
      t_bjet_matchpt  .push_back( std::get<4>(pfJetsIdx[ij]).Pt());
      t_bjet_matcheta .push_back( std::get<4>(pfJetsIdx[ij]).Eta());
      t_bjet_matchphi .push_back( std::get<4>(pfJetsIdx[ij]).Phi());
      t_bjet_matchmass.push_back( std::get<4>(pfJetsIdx[ij]).M());
      t_bjet_flavor.push_back( std::get<5>(pfJetsIdx[ij]) );
      t_bjet_flavorForB.push_back( std::get<6>(pfJetsIdx[ij]) );
    }


    t_ht  = pfht;
    t_mht = vis.Pt();
    
    // // now set the 4 variables that we added for the tmva reader for the bdt evaluation
    // bdt_l1pt      = t_lep_calpt[0];
    // bdt_apt       = (t_lep_calpt[0]-t_lep_calpt[1])/(t_lep_calpt[0]+t_lep_calpt[1]);
    // bdt_abslleta  = fabs(t_lleta);
    // bdt_dphilll2  = fabs(dphi_2(t_lep_calpt[0],t_lep_eta[0],t_lep_phi[0],t_lep_calpt[1],t_lep_eta[1],t_lep_phi[1],2)); // this function is in functions.cc in scripts/
    // bdt_sumabseta = fabs(t_lep_eta[0])+fabs(t_lep_eta[1]);
    // //bdt_flavor    = abs(t_lep_pdgId[0]*t_lep_pdgId[1]); //abs should be fine here, it's an int
    // t_apt         = bdt_apt;
    // t_dphilll2    = bdt_dphilll2;
    // t_bdt         = reader->EvaluateMVA( methodName );
    // t_bdt_rarity  = reader->GetRarity  ( methodName );
    // t_fisher2     = readerFisher2->EvaluateMVA( methodNameFisher2 );
    t_bjet_leadPassTight = (t_bjet_csvv2.size()>0 && t_bjet_csvv2[0]>csvWPList[1]);
    t_isData = !isMC;
    
    outTree->Fill();
  }

  //save histos to file  
  if(outURL!=""){
    TFile *fOut=TFile::Open(outURL,"RECREATE");
    fOut->cd();

    outTree->Write();

    //store the weight sum for posterior normalization
    TH1D *wgtH=new TH1D("wgtsum","wgtsum",1,0,1);
    wgtH->SetBinContent(1,wgtSum);
    wgtH->SetDirectory(fOut);
    wgtH->Write();

    TH1D *allwgtH=new TH1D("allwgtsum","allwgtsum",allWgtSum.size(),0,allWgtSum.size());
    for(size_t i=0; i<allWgtSum.size(); i++)
      allwgtH->SetBinContent(i+1,allWgtSum[i]);
    allwgtH->SetDirectory(fOut);
    allwgtH->Write();
    for (auto& it : ht.getPlots())  { 
      if(it.second->GetEntries()==0) continue;
      it.second->SetDirectory(fOut); it.second->Write(); 
    }
    for (auto& it : ht.get2dPlots())  { 
      if(it.second->GetEntries()==0) continue;
      it.second->SetDirectory(fOut); it.second->Write(); 
    }
    fOut->Close();
  }

  return 0;
}
