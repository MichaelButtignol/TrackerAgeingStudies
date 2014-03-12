// Code to read voltage root file get from ConditionBrowser web page
// Only for crosscheck of steps definition

#include "TString.h"
#include "TTimeStamp.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1F.h"

TString SecUTC(Int_t sec)
{
  TTimeStamp ts(sec, 0);
  TString s = ts.AsString("c");
  return s(0, 4) + "." +
  s(5, 2) + "." +
  s(8, 2) + " " +
  s(11, 2) + ":" +
  s(14, 2) + ":" +
  s(17, 2);
}

TGraph* ReadVoltage(char* filename="Data/ConditionBrowser_1343662723590.root")
{
  char NAME[200];
  Int_t NROW;
  Double_t WEIGHT[20000];
  Int_t TIME[20000];
  Double_t VALUE[20000];
  
  TFile f(filename);
  TTree *t = (TTree*)f.Get("tree");
  if(!t) return 0;
  t->SetBranchAddress("NAME", &NAME);
  t->SetBranchAddress("NROW", &NROW);
  t->SetBranchAddress("WEIGHT", WEIGHT);
  t->SetBranchAddress("TIME", TIME);
  t->SetBranchAddress("VALUE", VALUE);
  
  TGraph *g = new TGraph();
  g->SetName("Voltage");
  
  Int_t nentries = (Int_t)t->GetEntries();
  for (Int_t i=0; i<nentries; i++)
  {
    t->GetEntry(i);
    //cout << NAME << endl;
    //cout << NROW << endl;
    //cout << "ROW \tWEIGHT \tTIME \tVALUE" << endl;
    for (Int_t j=0; j<NROW; j++)
    {
      /*cout << (j+1) << " \t"
       << WEIGHT[j] << " \t"
       << SecUTC(TIME[j]) << " \t"
       << VALUE[j] << endl;*/
      g->SetPoint(j, TIME[j]-10800, VALUE[j]);
    }
  }
  f.Close();
  
  TH1F* h = g->GetHistogram();
  h->GetXaxis()->SetTimeDisplay(1);
  h->GetXaxis()->SetTimeFormat("%H:%M");
  g->SetMarkerStyle(21);
  
  return g;
}
