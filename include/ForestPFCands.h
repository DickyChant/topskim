#ifndef ForestPFCands_h
#define ForestPFCands_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>

class ForestPFCands {
public :
 ForestPFCands(TChain *t): nPF(0),
    pfId(0),
    pfPt(0),
    pfEta(0),
    pfPhi(0),
    pfM(0) // ,
    // There seems to be no Track info in Run3 HIForest PF candidates
    // trkAlgo(0),
    // trkPtError(0),
    // trkNHit(0),
    // trkChi2(0),
    // trkNdof(0)
      {     
        std::cout << "Setting up pfCand branches." << std::endl;
        t->SetBranchAddress("nPF", &nPF); 
        t->SetBranchAddress("pfId", &pfId);
        t->SetBranchAddress("pfPt", &pfPt);
        t->SetBranchAddress("pfEta", &pfEta);
        t->SetBranchAddress("pfPhi", &pfPhi);
        t->SetBranchAddress("pfM", &pfM);
        // t->SetBranchAddress("trkAlgo", &trkAlgo);
        // t->SetBranchAddress("trkPtError", &trkPtError);
        // t->SetBranchAddress("trkNHit", &trkNHit);
        // t->SetBranchAddress("trkChi2", &trkChi2);
        // t->SetBranchAddress("trkNdof", &trkNdof);
      }

  ~ForestPFCands() {}

   int nPF;
   std::vector<int>     *pfId;
   std::vector<float>   *pfPt;
   std::vector<float>   *pfEta;
   std::vector<float>   *pfPhi;
   std::vector<float>   *pfM;
  //  std::vector<int>     *trkAlgo;
  //  std::vector<float>   *trkPtError;
  //  std::vector<float>   *trkNHit;
  //  std::vector<float>   *trkChi2;
  //  std::vector<float>   *trkNdof;
};



#endif 
