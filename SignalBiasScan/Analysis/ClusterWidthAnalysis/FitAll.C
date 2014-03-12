#include <iostream>
#include <time.h>
#include <sstream>
#include <map>
#include <sstream>
#include <utility>

#include <TString.h>
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TRint.h>
#include <exception>
#include <TTree.h>

#include "Fit.C"

void FitAll(string SubDet="TIB", string Run="_190459")
{
  int err_type=0;
  const int step = 55;

  string FileToOpen = SubDet+"_output_DecoMode"+Run+".root";
  TString Output = "DECO_ClusterWidth_"+SubDet+".root";

  std::cout << "Opening input file " << FileToOpen << std::endl;

  TFile* myFile = TFile::Open(FileToOpen.c_str());
  TTree* tr = (TTree*)myFile->FindObjectAny("T");
  std::cout << "tree opened" << std::endl;

  Int_t detid; // Int_t for old files
  //ULong64_t detid; // Int_t for old files
  Double_t volt;
  Int_t id;
  Double_t evolt;
  Double_t mpv;
  Double_t empv;
  ULong64_t tempdetid;
  Double_t chi2overndf;
    
  tr->SetBranchAddress("DetID",&detid);
  tr->SetBranchAddress("Voltage",&volt);
  tr->SetBranchAddress("Index",&id);
  tr->SetBranchAddress("errVoltage",&evolt);
  tr->SetBranchAddress("Mean",&mpv);
  tr->SetBranchAddress("errMean",&empv);
 

  std::cout << "Opening output file " << Output << std::endl;

  TFile* output = new TFile(Output,"recreate");
  TTree* tout = new TTree("tout","SignalSummary");

  ULong64_t odetid;
  int olayer;
  double odepvolt;
  double oerrdepvolt;
  double oplateau;
  double ofitchisquare;
  int ofitstatus;
  double olastpty;
  double ochi2;

  tout->Branch("DETID",&odetid,"DETID/l");
  tout->Branch("LAYER",&olayer,"LAYER/I");
  tout->Branch("ERRDEPVOLT",&oerrdepvolt,"ERRDEPVOLT/D");
  tout->Branch("DEPVOLT",&odepvolt,"DEPVOLT/D");
  tout->Branch("PLATEAU",&oplateau,"PLATEAU/D");
  tout->Branch("FITCHI2",&ofitchisquare,"FITCHI2/D");
  tout->Branch("FITSTATUS",&ofitstatus,"FITSTATUS/I");
  tout->Branch("LASTPOINTS",&olastpty,"LASTPOINTS/D");
  tout->Branch("CHI2",&ochi2,"CHI2/D");

  UInt_t nentries = tr->GetEntries();
  Int_t nfit = 0;
  Int_t nbadfit = 0;
  Int_t status = 0;
  double lastpty = 0;
  TCanvas *c1;
  TCanvas *c2;

  double Vdep;
  double errVdep;

  typedef std::map< ULong64_t, std::vector<double> > mapping;
  std::vector< double > tempdepvolt;
  mapping DetID_Vdep;

  double avolt[step];
  int aid[step];
  double aevolt[step];
  double ampv[step];
  double aempv[step];
  UInt_t k=0;
  
  Int_t layer = 0;
  float mean, rms;

  std::cout << nentries << std::endl;

  // Loop over entries to create curves
  //  for(Int_t i = 0; i < nentries; i++){
  for(UInt_t i = 0; i <nentries; i++){
    tr->GetEntry(i);
    tempdetid= (ULong64_t) detid;
	//if(detid!=369121606) continue;
	
    mapping::iterator iter = DetID_Vdep.find(tempdetid);
    if(iter == DetID_Vdep.end() ){
      k=0;
      for(UInt_t j = i; j< nentries; j++){
	    tr->GetEntry(j);

		if(tempdetid==(ULong64_t) detid){
		  if(empv<5){
	    	avolt[k]=volt;
	    	aid[k]=id;
	    	aevolt[k]=evolt;
	    	ampv[k]=mpv;
	    	aempv[k]=empv;
			if(err_type==1) aempv[k]=2.5;
	
	    	k++;

		  }
		}
	  
	  }//for(UInt_t j = i; j< nentries; j++)

	  if(k<3) continue; // Need enought points to compute 2nd derivative.
	  
	  if(k) lastpty = ampv[k-1];
	  if(k>2) lastpty = (ampv[k-1]+ampv[k-2]+ampv[k-3])/3.;

      TString canvasname;
      canvasname.Form("CW_%llu",tempdetid);
      c1 = new TCanvas(canvasname);
      c1->cd();
      TGraphErrors * thegraph = new TGraphErrors(k, avolt, ampv, aevolt, aempv);
      string corr_name="_"+SubDet+Run;
      int corrected = CorrectGraphForLeakageCurrent(thegraph, tempdetid, corr_name);

      GetYMeanRMS(thegraph, mean, rms);
	  //cout<<"rms "<<rms<<endl;
	  if(rms<0.1) continue;

      thegraph->SetLineColor(2);
      thegraph->SetMarkerColor(1);
      thegraph->SetMarkerStyle(20);
	  
	  layer = GetLayer(tempdetid);
	  
      std::cout << tempdetid <<" layer "<< layer << std::endl;
	  Vdep = FitCurve( thegraph );
	  
	  
	  tempdepvolt.clear();
      tempdepvolt.push_back(layer);
      tempdepvolt.push_back(Vdep);
      tempdepvolt.push_back(errVdep);
	  tempdepvolt.push_back(0);
      tempdepvolt.push_back(0);
	  tempdepvolt.push_back(0);
	  tempdepvolt.push_back(lastpty);
	  tempdepvolt.push_back(0);
      DetID_Vdep.insert(std::pair< ULong64_t , std::vector<double> >(tempdetid,tempdepvolt));
	  //if(status!=0)
	  //if(tempdetid==369121606 || tempdetid==369121365) 
	  //{ gDirectory->Append(c1); gDirectory->Append(c2);}

    }//if(iter == DetID_Vdep.end() )
  
  }//for(UInt_t i = 0; i <nentries; i++)

  for(mapping::iterator iter = DetID_Vdep.begin(); iter != DetID_Vdep.end(); ++iter){
    olayer = iter->second.at(0);
    odepvolt = iter->second.at(1);
    oerrdepvolt = iter->second.at(2);
	oplateau = iter->second.at(3);
    ofitchisquare = iter->second.at(4);
	ofitstatus = iter->second.at(5);
	olastpty = iter->second.at(6);
	ochi2 = iter->second.at(7);
    odetid = iter->first;
    tout->Fill();
  }
  
  std::cout<<nbadfit<<" bad fits over "<<nfit<<std::endl;
  
  output->cd();
  tout->Write();
  gDirectory->Write();
  output->Close();
  myFile->Close();

  std::cout << "Closing input file " << FileToOpen << std::endl;
  std::cout << "Closing output file " << Output << std::endl;


}


int main(int argc, char *argv[]) {
  std::cout << "Starting GetKinkForAll.C" << std::endl;
  FitAll();
  
  return 0;
}
