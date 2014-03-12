//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 20 17:14:55 2010 by ROOT version 5.22/00d
// from TTree ttree/ttree
// found on file: clustersTree.root
//////////////////////////////////////////////////////////

#ifndef BackplaneCorrectionAnalysisTreeMaker_h
#define BackplaneCorrectionAnalysisTreeMaker_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <map>
#include <vector>

#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>

#include "../../../interface/TreeEvent.h"
#include "../../VoltageSteps/Code/VoltageStepsMaker.h"

class BackplaneCorrectionAnalysisTreeMaker {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   TreeEvent		*event;
   TBranch			*branch;

   BackplaneCorrectionAnalysisTreeMaker(TTree *tree=0, std::string file="Run.txt", int subdet=1, bool timestamp=false,
    int angle = 0, int width = 0, int onstrip = 0 );
   virtual ~BackplaneCorrectionAnalysisTreeMaker();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
   
   VoltageStepsMaker VSmaker;   

   std::map<ULong64_t , std::vector<TH1F*> > HistSoN_TIB;
   std::vector< std::map<ULong64_t , std::vector<TH1F*> > > HistSoN_TIB_onstrip;
   std::map<ULong64_t , std::vector<TH1F*> > HistSoN_TOB;
   std::map<ULong64_t , std::vector<TH1F*> > HistSoN_TID;
   std::vector< std::map<ULong64_t , std::vector<TH1F*> > > HistSoN_TID_onstrip;
   std::map<ULong64_t , std::vector<TH1F*> > HistSoN_TEC;
   
   std::map<ULong64_t , std::vector<TProfile*> > ProfVsAngle_TIB;
   std::map<ULong64_t , std::vector<TProfile*> > ProfDeltaxVsAngle_TIB;
   std::map<ULong64_t , std::vector<TProfile*> > ProfVsAngle_TOB;
   std::map<ULong64_t , std::vector<TProfile*> > ProfVsAngle_TID;
   std::map<ULong64_t , std::vector<TProfile*> > ProfVsAngle_TEC;

   unsigned int t_monitor_start;
   unsigned int t_monitor_end;
   std::map<ULong64_t, TProfile* > Monitors_TIB;
   std::map<ULong64_t, TProfile* > Monitors_TOB;
   std::map<ULong64_t, TProfile* > Monitors_TID;
   std::map<ULong64_t, TProfile* > Monitors_TEC;
   
   void FillHitInfo(std::map<ULong64_t , std::vector<TH1F*> > &, std::map< ULong64_t , std::vector<TProfile*> > &, 
                    TreeHit *, int , bool , std::vector< TH1F* > , std::vector< TH2F* >);
   void FillHitInfo(std::map<ULong64_t , std::vector<TH1F*> > &, std::map< ULong64_t , std::vector<TProfile*> > &, 
                    TreeFullHit *, int , bool , std::vector< TH1F* >, std::vector< TH2F* > );
   void FillHitInfo(std::map<ULong64_t , std::vector<TH1F*> > &, std::map< ULong64_t , std::vector<TProfile*> > &, 
                    TreeHit *, int , bool , std::vector< TH1F* > , std::vector< TH2F* >, float, float , float );

   void FillHistos(std::map< ULong64_t , std::vector<TH1F*> > &, std::map< ULong64_t , std::vector<TProfile*> > &,
                   std::vector<TreeHit> &, int , bool , std::vector< TH1F* >, std::vector< TH2F* >, std::map<ULong64_t, TProfile* > );
   void FillHistos(std::map< ULong64_t , std::vector<TH1F*> > &, std::map< ULong64_t , std::vector<TProfile*> > &,
                   std::vector<TreeFullHit> &, int , bool , std::vector< TH1F* >, std::vector< TH2F* >, std::map<ULong64_t, TProfile* > );
   
   void FitHistos(std::map< ULong64_t , std::vector<TH1F*> > &, std::string , std::vector< TH1F* >, std::vector< TH2F* >,std::map<ULong64_t, TProfile* >);
   void FitProfiles(std::map< ULong64_t , std::vector<TProfile*> > &, std::string);

   enum SubDet { All, TIB, TOB, TID, TEC};

   const int part;
   const std::string stepsfile;
   const bool usetimestamp;
   const int use_angle;
   const int use_width;
   int use_onstrip;

};

#endif

#ifdef BackplaneCorrectionAnalysisTreeMaker_cxx
BackplaneCorrectionAnalysisTreeMaker::BackplaneCorrectionAnalysisTreeMaker(TTree *tree, std::string file,
 int subdet, bool timestamp, int angle, int width, int onstrip):stepsfile(file),part(subdet),usetimestamp(timestamp),
 use_angle(angle), use_width(width), use_onstrip(onstrip)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

   event = 0;

   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("clustersTree.root");
      if (!f) {
         f = new TFile("clustersTree.root");
         f->cd("clustersTree.root:/demo");
      }
      tree = (TTree*)gDirectory->Get("ttree");

   }
   Init(tree);
   
   for(int i=0; i<6; i++){
     std::map<ULong64_t , std::vector<TH1F*> > histomap;
     HistSoN_TIB_onstrip.push_back( histomap );
     HistSoN_TID_onstrip.push_back( histomap );
   }
   
   t_monitor_start = 0;
   t_monitor_end = 0;
}

BackplaneCorrectionAnalysisTreeMaker::~BackplaneCorrectionAnalysisTreeMaker()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BackplaneCorrectionAnalysisTreeMaker::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BackplaneCorrectionAnalysisTreeMaker::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void BackplaneCorrectionAnalysisTreeMaker::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   //fChain->SetMakeClass(1);
   
   fChain->SetBranchAddress("event", &event);

   Notify();
}

Bool_t BackplaneCorrectionAnalysisTreeMaker::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BackplaneCorrectionAnalysisTreeMaker::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BackplaneCorrectionAnalysisTreeMaker::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BackplaneCorrectionAnalysisTreeMaker_cxx