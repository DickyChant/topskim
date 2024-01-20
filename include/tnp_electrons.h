#ifndef tnp_electrons_h
#define tnp_electrons_h

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TObjString.h"
#include <map>

/**
   @short reads scale factors for ID/HLT for electrons
 */
class ElectronEfficiencyWrapper
{

 public:

  /**
     @short CTOR - just parses all the files from the directory given
   */
  ElectronEfficiencyWrapper(TString url,bool useOldId=true) 
    {
      TString baseF("ScaleFactors_PbPb_LooseWP");
      TString regs[2]={"EB","EE"};
      TString centr[2]={"0_30","30_100"};
      TString pfix[3]={useOldId ? "_AlpaFixedDataOnly_BWResCBErfExp_preliminary_v2" : "" ,"_HLT","_RECO"};
      for(size_t i=0; i<2; i++) {
        for(size_t j=0; j<2; j++) {
          for(size_t k=0; k<3; k++) {
            TString path(Form("%s/%s_%s_Centr_%s%s.root",
                              url.Data(),
                              baseF.Data(),
                              regs[i].Data(),
                              centr[j].Data(),
                              pfix[k].Data()));
            gSystem->ExpandPathName(path);
            TFile *f=TFile::Open(path);
            TString key(regs[i]+centr[j]+pfix[k]);
	    TObjArray *tx = key.Tokenize("_");
	    for (Int_t t = 0; t < tx->GetEntries(); t++) {
	      if (t==0) {key.Clear(); key.Append(((TObjString *)(tx->At(t)))->String()+"_");}
	      else if (t<=2) if (t==2 and !pfix[k].Contains("HLT") and !pfix[k].Contains("RECO")) continue; else key.Append(((TObjString *)(tx->At(t)))->String()+"_");
	      else continue; 
	    }
	    sfs_[key]=(TGraphAsymmErrors *)f->Get("g_scalefactors");
            f->Close();
          }
        }
      }
    }

  /**
     @short returns the SF and the uncertainty
   */
  std::pair<float,float> eval(float pt, bool isEB, int cenbin, bool hlt, bool reco) {

    std::pair<float,float> sfVal(1.0,0.0);

    //build the key to the map
    TString reg(isEB ? "EB" : "EE");
    TString cen(cenbin<30 ? "0_30" : "30_100");
    TString pfix("_");
    if(hlt) pfix.Append("HLT_");
    if(reco) pfix.Append("RECO_");
    TString key(reg+cen+pfix);
    //check key exists
    if(sfs_.find(key)==sfs_.end()) {
      std::cout << "Unable to find " << key << " in electron SFs map..." << std::endl;
      return sfVal;
    }

    //find closest value in x (prefer to spline interpolation)
    Double_t x,y,mindiff(9.e+6);
    for(int i=0; i<sfs_[key]->GetN(); i++){
      sfs_[key]->GetPoint(i,x,y);
      if(fabs(x-pt)>mindiff) continue;
      mindiff=fabs(x-pt);
      sfVal.first=float(y);
      sfVal.second=sfs_[key]->GetErrorY(i);
    }

    return sfVal;
  }

 private:
  std::map< TString, TGraphAsymmErrors *> sfs_;

};



#endif
