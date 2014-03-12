#include <string>
#include <TTree.h>
#include <TBranch.h>
#include <TGraphErrors.h>
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "iostream"

//#include CurvesFunctions.C

using namespace std;

TGraphErrors* GetGraph(string subdet, string run, ULong64_t modid, int err_type)
{

  TGraphErrors* g;

  //Int_t detid_2;
  ULong64_t detid_2;
  double volt_2;
  int id_2;
  double evolt_2;
  double mpv_2;
  double empv_2;
  int tempdetid_2;
  double chi2overndf_2;
  double slope_2;
  double eslope_2;
  double origin_2;
  double eorigin_2;  

  const int step = 200;
  double avolt[step];
  int aid[step];
  double aevolt[step];
  double ampv[step];
  double aempv[step];
  double aslope[step];
  double aeslope[step];
  double aorigin[step];
  double aeorigin[step];  

//  char* file = Form("%s_output_DecoMode%s.root", subdet.c_str(), run.c_str());
  char* file = Form("%s_output_DecoMode_prof%s.root", subdet.c_str(), run.c_str());
  TFile* f = TFile::Open(file);
  if(!f) {cout<<"Error : no file '"<<file<<"'"<<endl; return g;}

  TTree* tr = (TTree*) f->FindObjectAny("T");
  tr->SetBranchAddress("DetID",&detid_2);
  tr->SetBranchAddress("Voltage",&volt_2);
  tr->SetBranchAddress("Index",&id_2);
  tr->SetBranchAddress("errVoltage",&evolt_2);
  tr->SetBranchAddress("Slope",&slope_2);
  tr->SetBranchAddress("errSlope",&eslope_2);
  tr->SetBranchAddress("Origin",&origin_2);
  tr->SetBranchAddress("errOrigin",&eorigin_2);    
//  tr->SetBranchAddress("MPV",&mpv_2);
//  tr->SetBranchAddress("errMPV",&empv_2);
  tr->SetBranchAddress("Chi2OverNdf",&chi2overndf_2);

  UInt_t nentries = tr->GetEntries();
  ULong64_t modid_init = modid;
  int k=0;
  for(unsigned int j = 0; j <nentries; j++)
  {
	//cout << "j= " << j << " over " << nentries << " entries!" << endl;
  	tr->GetEntry(j);
  	//cout << "detId_2= " << detid_2 << " and modid= " << modid << endl;
  	
/*		modid = modid_init;   //if we want to compute both of the two sensors of the TOB module in an only graph 
  		if(subdet == "TOB")
  		{
  			modid = 10*modid + 1;	//To switch for example from 436281508 to 4362815081 (for TOB)
  			if(detid_2!=modid) modid += 1;  //To switch for example from 4362815081 to 4362815082 (the other part of the same module)
  		}*/
  		
  	
  	if(k>0 && detid_2!=modid) break;
  		
	if(detid_2==modid)
	{
		if(empv_2<5)
		{
			avolt[k]=volt_2;
			aid[k]=id_2;
	    		
	    		aevolt[k]=evolt_2;
	    		//ampv[k]=mpv_2;
	    		//aempv[k]=empv_2;
	    		aslope[k]=slope_2;
	    		aeslope[k]=eslope_2;
	    		aorigin[k]=origin_2;
	    		aeorigin[k]=eorigin_2;
		
			if(err_type==1) aempv[k]=2.5;
			if(err_type==2) aempv[k]=5;
			if(err_type==3) aempv[k]=empv_2*5;
			if(err_type==4) aempv[k]=empv_2*10;
			if(err_type==5) aempv[k]=10;
	    		k++;
		}
		else cout << "WARNING: empv > 5!! Risk of leak... " << endl;
	}
  }

  f->Close();

  //cout << " | avolt[1]= " << avolt[1] << " | aslope[1]= " << aslope[1] << endl;
  if(!k) return 0;
  
  // tell it here what variable you want to see, in which module, from what run and with what kind of error 
  g = new TGraphErrors(k, avolt, aslope, aevolt, aeslope);
  //g = new TGraphErrors(k, avolt, ampv, aevolt, aempv);  
  g->SetLineColor(2);
  g->SetMarkerStyle(20);

  return g;
}

void ShowOneCurve(string subdet, string run, ULong64_t modid, bool showfit, bool print)
{

 int err_type=0;

 // Get graph
 TGraphErrors* g = GetGraph(subdet, run, modid, err_type);
 
 // Draw graph and function
 TCanvas c1("c1","",0,0,700,500);
 TH1F* h = g->GetHistogram();
 h->SetMinimum(-0.005);
// if(subdet=="TIB") h->SetMaximum(100);//100//140//65
 if(subdet=="TIB") h->SetMaximum(0.0005);
 if(subdet=="TOB") h->SetMaximum(0.005);
 g->SetTitle(Form("Detid %llu",modid));
 g->Draw("ALP");
 
 c1.Modified();
 c1.Update();

 if(print) c1.Print(Form("detid_%llu%s.eps", modid, run.c_str()));
 getchar();
 c1.Close();
}

void ShowTIBCurves_SmallScan(string run, bool showfit=true, bool print=false)
{
  //TIBminus_1_2_2_1
  ShowOneCurve("TIB", run, 369121381, showfit, print);
  ShowOneCurve("TIB", run, 369121382, showfit, print);
  ShowOneCurve("TIB", run, 369121385, showfit, print);
  ShowOneCurve("TIB", run, 369121386, showfit, print);
  ShowOneCurve("TIB", run, 369121389, showfit, print);
  ShowOneCurve("TIB", run, 369121390, showfit, print);
 //TIBminus_1_4_2_5
/*  ShowOneCurve("TIB", run, 369121605, showfit, print);
  ShowOneCurve("TIB", run, 369121606, showfit, print);
  ShowOneCurve("TIB", run, 369121609, showfit, print);
  ShowOneCurve("TIB", run, 369121610, showfit, print);
  ShowOneCurve("TIB", run, 369121613, showfit, print);
  ShowOneCurve("TIB", run, 369121614, showfit, print);
*/ //TIBplus_1_6_2_5
//  ShowOneCurve("TIB", run, 369125861, showfit, print);
//  ShowOneCurve("TIB", run, 369125862, showfit, print);
//  //ShowOneCurve("TIB", run, 369125865, showfit, print);
//  ShowOneCurve("TIB", run, 369125866, showfit, print);
//  ShowOneCurve("TIB", run, 369125869, showfit, print);
//  ShowOneCurve("TIB", run, 369125870, showfit, print);

}

void ShowTOBCurves_SmallScan(string run, bool showfit=true, bool print=false)
{
  //TOB
  ShowOneCurve("TOB", run, 4362815081, showfit, print);
  ShowOneCurve("TOB", run, 4362815082, showfit, print);
  ShowOneCurve("TOB", run, 4362815121, showfit, print);
  ShowOneCurve("TOB", run, 4362815122, showfit, print);
  ShowOneCurve("TOB", run, 4362815161, showfit, print);
  ShowOneCurve("TOB", run, 4362815162, showfit, print);
  ShowOneCurve("TOB", run, 4362815201, showfit, print);
  ShowOneCurve("TOB", run, 4362815202, showfit, print);
  ShowOneCurve("TOB", run, 4362815241, showfit, print);
  ShowOneCurve("TOB", run, 4362815242, showfit, print);
  ShowOneCurve("TOB", run, 4362815281, showfit, print);
  ShowOneCurve("TOB", run, 4362815282, showfit, print);    
}

void ShowCurves()
{
// gROOT->ProcessLine(".L signalFunction.cpp");
// gROOT->ProcessLine(".L fitSignal.cpp");
// gSystem->Load("signalFunction.so");

 //gROOT->ProcessLine(".L CurvesFunctions.C");

 string subdet = "TOB";

// string run = "_20120928_run203832"; 
 string run = "_20110315_run160497";
 bool showfit=true;
 bool print=false;
 
 cout<<endl<<"Run "<<run<<endl<<endl;

 if(subdet=="TIB") ShowTIBCurves_SmallScan(run, showfit, print); // TIB
 if(subdet=="TOB") ShowTOBCurves_SmallScan(run, showfit, print); // TOB

}
