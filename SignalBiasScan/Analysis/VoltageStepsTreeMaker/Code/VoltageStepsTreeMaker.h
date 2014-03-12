//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 20 17:14:55 2010 by ROOT version 5.22/00d
// from TTree ttree/ttree
// found on file: clustersTree.root
//////////////////////////////////////////////////////////

#ifndef VoltageStepsTreeMaker_h
#define VoltageStepsTreeMaker_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <map>
#include <vector>
#include <set>

#include <TH1F.h>
#include <TProfile.h>

#include "../../../interface/TreeEvent.h"

class VoltageStepsTreeMaker {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   TreeEvent		*event;
   TBranch			*branch;

   VoltageStepsTreeMaker(TTree *tree=0, std::string file="Run.txt", int subdet=1, bool timestamp=false,
    int angle = 0, int width = 0, int onstrip = 0 );
   virtual ~VoltageStepsTreeMaker();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
      
   std::set< int > Voltage;
   std::multimap< int, std::pair< int, std::pair< int,int > > > VoltageSteps_evtnumber;
   void setVoltageSteps_evtnumber(const char*);
   int getVoltageSteps_evtnumber(int run, int evt, int timestamp=0);

   std::multimap< int, std::pair< int,int > > VoltageSteps_timestamp;
   void setVoltageSteps_timestamp(const char*);
   int getVoltageSteps_timestamp(int timestamp);
   std::map< int, int > VoltageSteps_timestart;
   std::map< int, int > VoltageSteps_timeend;

   void setVoltageSteps(const char* inputfile, bool usetimestamp=false);

   std::map<ULong64_t , std::vector<TH1F*> > HistSoN_TIB;
   std::vector< std::map<ULong64_t , std::vector<TH1F*> > > HistSoN_TIB_onstrip;
   std::map<ULong64_t , std::vector<TH1F*> > HistSoN_TOB;
   std::map<ULong64_t , std::vector<TH1F*> > HistSoN_TID;
   std::vector< std::map<ULong64_t , std::vector<TH1F*> > > HistSoN_TID_onstrip;
   std::map<ULong64_t , std::vector<TH1F*> > HistSoN_TEC;
   
   unsigned int t_monitor_start;
   unsigned int t_monitor_end;
   std::map<ULong64_t, TProfile* > Monitors_TIB;
   std::map<ULong64_t, TProfile* > Monitors_TOB;
   std::map<ULong64_t, TProfile* > Monitors_TID;
   std::map<ULong64_t, TProfile* > Monitors_TEC;
   
   void FillHitInfo(std::map<ULong64_t , std::vector<TH1F*> > &, TreeHit *, int , bool , std::vector< TH1F* > );
   void FillHitInfo(std::map<ULong64_t , std::vector<TH1F*> > &, TreeFullHit *, int , bool , std::vector< TH1F* > );
   void FillHitInfo(std::map<ULong64_t , std::vector<TH1F*> > &, TreeHit *, int , bool , std::vector< TH1F* > ,
     float, float , float );

   void FillHistSoN(std::map< ULong64_t , std::vector<TH1F*> > &, std::vector<TreeHit> &, int , bool , std::vector< TH1F* >, std::map<ULong64_t, TProfile* > );
   void FillHistSoN(std::map< ULong64_t , std::vector<TH1F*> > &, std::vector<TreeFullHit> &, int , bool , std::vector< TH1F* >, std::map<ULong64_t, TProfile* > );
   
   void FitHistos(std::map< ULong64_t , std::vector<TH1F*> > &, std::string , std::vector< TH1F* >, std::map<ULong64_t, TProfile* >);

enum SubDet { All, TIB, TOB, TID, TEC};

   const int part;
   const std::string stepsfile;
   const bool usetimestamp;
   const int use_angle;
   const int use_width;
   int use_onstrip;

};

#endif

#ifdef VoltageStepsTreeMaker_cxx
VoltageStepsTreeMaker::VoltageStepsTreeMaker(TTree *tree, std::string file,
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

VoltageStepsTreeMaker::~VoltageStepsTreeMaker()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t VoltageStepsTreeMaker::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t VoltageStepsTreeMaker::LoadTree(Long64_t entry)
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

void VoltageStepsTreeMaker::Init(TTree *tree)
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

Bool_t VoltageStepsTreeMaker::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void VoltageStepsTreeMaker::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t VoltageStepsTreeMaker::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef VoltageStepsTreeMaker_cxx
