#define ClusterWidthAnalysisTreeMaker_cxx

#include "ClusterWidthAnalysisTreeMaker.h"

#include <iostream>
#include <sstream>

#include <TString.h>
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TMath.h>
#include <TRint.h>


// 1 = TIB // 2 = TOB // 3 = TID // 4 = TEC
// ****************************************

enum SubDet { All, TIB, TOB, TID, TEC};

using namespace std;


// Main method //
//-------------//
void ClusterWidthAnalysisTreeMaker::Loop()
{
  
  std::cout<<"Initializing"<<std::endl;
  
  std::string output_file;
  
  // Read definition of voltage steps
  std::cout << "Setting Deco Voltage vector : " << std::endl;
  if(usetimestamp) std::cout << " using timestamp step definition " << std::endl;
  else std::cout << " using event number step definition" << std::endl;
  
  //if(usetimestamp) VSmaker.readVoltageSteps_timestamp(stepsfile.c_str());
  //else VSmaker.readVoltageSteps_evtnumber(stepsfile.c_str());
  VSmaker.readVoltageSteps(stepsfile.c_str(), usetimestamp);
  std::cout << " done" << std::endl;
  
  VSmaker.Initialize();
  t_monitor_start = VSmaker.t_monitor_start;
  t_monitor_end = VSmaker.t_monitor_end;

  int nevent    = 0;
  int theVoltage=-1;
  
  if (fChain == 0) return;
  
  // Define histos
  std::vector< std::vector< TH1F* > > commonHistos;
  std::vector< std::vector< TH1F* > > emptyHistos;
  string subdetname[] = {"TIB", "TOB", "TID", "TEC"};
    
  for(unsigned int idet=0; idet<4; idet++)
  {
	std::vector< TH1F* > Histos;
	Histos.push_back(new TH1F(Form("hNtracks_%s", subdetname[idet].c_str()), "hNtracks", 1000, 0, 1000));
	Histos.push_back(new TH1F(Form("hTrackAngle_%s", subdetname[idet].c_str()), "hTrackAngle", 90, -180, 180));
	Histos.push_back(new TH1F(Form("hTsosxMinusClusx_%s", subdetname[idet].c_str()), "hTsosxMinusClusx", 100, -0.01, 0.01));
	Histos.push_back(new TH1F(Form("hClusMinusSeed_%s", subdetname[idet].c_str()), "hClusMinusSeed", 150, -1.5, 1.5));
	Histos.push_back(new TH1F(Form("hNevtPerStep_%s", subdetname[idet].c_str()), "hNevtPerStep", 80, 0, 400));
    commonHistos.push_back(Histos);
  }
  
  // Define monitors
  Monitors_TIB.insert(std::pair< ULong64_t, TProfile* > (369121605, new TProfile("monitor_369121605", "", (t_monitor_end-t_monitor_start)/30, 0, t_monitor_end-t_monitor_start)));
  Monitors_TIB.insert(std::pair< ULong64_t, TProfile* > (369121390, new TProfile("monitor_369121390", "", (t_monitor_end-t_monitor_start)/30, 0, t_monitor_end-t_monitor_start)));
  Monitors_TIB.insert(std::pair< ULong64_t, TProfile* > (369125870, new TProfile("monitor_369125870", "", (t_monitor_end-t_monitor_start)/30, 0, t_monitor_end-t_monitor_start)));
  Monitors_TOB.insert(std::pair< ULong64_t, TProfile* > (436281512, new TProfile("monitor_436281512", "", (t_monitor_end-t_monitor_start)/30, 0, t_monitor_end-t_monitor_start)));
  Monitors_TID.insert(std::pair< ULong64_t, TProfile* > (402664070, new TProfile("monitor_402664070", "", (t_monitor_end-t_monitor_start)/30, 0, t_monitor_end-t_monitor_start)));
  Monitors_TEC.insert(std::pair< ULong64_t, TProfile* > (470148196, new TProfile("monitor_470148196", "", (t_monitor_end-t_monitor_start)/30, 0, t_monitor_end-t_monitor_start)));
  Monitors_TEC.insert(std::pair< ULong64_t, TProfile* > (470148300, new TProfile("monitor_470148300", "", (t_monitor_end-t_monitor_start)/30, 0, t_monitor_end-t_monitor_start)));

  Long64_t nentries = fChain->GetEntriesFast();
  std::cout << "NUMBER OF ENTRIES = "<< nentries << std::endl; 


  //------------------//
  // Loop over events //
  //------------------//
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
	fChain->GetEntry(jentry);
    
	nevent++;
	if(nevent%1000 == 0) 
    { 
       std::cout << "number of events " << nevent << std::endl;
       std::cout<<"run "<<event->run_nr<<" event "<<event->ev_nr<<" timestamp "<<event->ev_timestamp<<" V "<<theVoltage<<std::endl;
    }   

    // Get voltage setting for this event
	theVoltage=-1;
	if(usetimestamp) theVoltage = VSmaker.getVoltage_timestamp(event->ev_timestamp);
	else theVoltage = VSmaker.getVoltage_evtnumber(event->run_nr, event->ev_nr, event->ev_timestamp);
	if(theVoltage<0) continue; // skip event if not on a voltage step

	int thePointNumber = VSmaker.getIndex(theVoltage);

	// check index
	if(thePointNumber < 0){
      std::cout<<" WARNING : point number out of range : "<<thePointNumber<<std::endl;
	  return;
	}


	for(unsigned int idet=0; idet<4; idet++)
	{
	  //commonHistos[idet][0]->Fill(event->tracks.size()); // N tracks per event
	  commonHistos[idet][0]->Fill(event->Ntracks); // adapted to v1.1 data format
	  commonHistos[idet][4]->Fill(theVoltage);
	}
 
    // Loop over tracks
    for(unsigned int itr=0; itr<event->tracks.size(); itr++)
	{
	  TreeTrack *track = &(event->tracks[itr]);
	  if(track->chi2>5) continue; // adapted to v1.2 data format
	  //if(track->pT<5) continue;

	  // Get number of hits
      int nTIBhits, nTOBhits, nTIDhits, nTEChits;
	  
	  if(track->TIB_hits.size()) nTIBhits = track->TIB_hits.size();
	  else nTIBhits = track->TIB_fullHits.size(); // same hits with more infos
	  if(track->TOB_hits.size()) nTOBhits = track->TOB_hits.size();
	  else nTOBhits = track->TOB_fullHits.size();
	  if(track->TID_hits.size()) nTIDhits = track->TID_hits.size();
	  else nTIDhits = track->TID_fullHits.size();
	  if(track->TEC_hits.size()) nTEChits = track->TEC_hits.size();
	  else nTEChits = track->TEC_fullHits.size();

      // int nHitsTotal = nTIBhits + nTOBhits + nTIDhits + nTEChits;
      int nHitsTotal = track->Nhits; // adapted to v1.1 data format
      if(nHitsTotal < 5) continue; // remove tracks with less than 5 hits


	  // Fill histos
	  
	  // TIB
	  if((part==TIB || part==All) && !use_onstrip)
	  {
	    if(track->TIB_hits.size()) FillHistos(HistSoN_TIB, ProfVsAngle_TIB, track->TIB_hits, thePointNumber, false, commonHistos[0], Monitors_TIB);
	    else if(track->TIB_fullHits.size()) FillHistos(HistSoN_TIB, ProfVsAngle_TIB, track->TIB_fullHits, thePointNumber, false, commonHistos[0], Monitors_TIB);
      }

	  if((part==TIB || part==All) && use_onstrip!=0)
	  {
		use_onstrip=1;
		FillHistos(HistSoN_TIB_onstrip[0], ProfVsAngle_TIB, track->TIB_fullHits, thePointNumber, false, commonHistos[0], Monitors_TIB);
		use_onstrip=2;
		FillHistos(HistSoN_TIB_onstrip[1], ProfVsAngle_TIB, track->TIB_fullHits, thePointNumber, false, commonHistos[0], Monitors_TIB);
		use_onstrip=3;
		FillHistos(HistSoN_TIB_onstrip[2], ProfVsAngle_TIB, track->TIB_fullHits, thePointNumber, false, commonHistos[0], Monitors_TIB);
		use_onstrip=1;
	  }

      // TOB
	  if(part==TOB || part==All)
	  {
	    if(track->TOB_hits.size()) FillHistos(HistSoN_TOB, ProfVsAngle_TOB, track->TOB_hits, thePointNumber, true, commonHistos[1], Monitors_TOB);
	    else if(track->TOB_fullHits.size()) FillHistos(HistSoN_TOB, ProfVsAngle_TOB, track->TOB_fullHits, thePointNumber, true, commonHistos[1], Monitors_TOB);
      }

      // TID
	  if((part==TID || part==All) && !use_onstrip)
	  {
	    if(track->TID_hits.size()) FillHistos(HistSoN_TID, ProfVsAngle_TID, track->TID_hits, thePointNumber, false, commonHistos[2], Monitors_TID);
	    else if(track->TID_fullHits.size()) FillHistos(HistSoN_TID, ProfVsAngle_TID, track->TID_fullHits, thePointNumber, false, commonHistos[2], Monitors_TID);	  
	  }
	  
	  if((part==TID || part==All) && use_onstrip!=0)
	  {
		use_onstrip=1;
		FillHistos(HistSoN_TID_onstrip[0], ProfVsAngle_TID, track->TID_fullHits, thePointNumber, false, commonHistos[2], Monitors_TID);
		use_onstrip=2;
		FillHistos(HistSoN_TID_onstrip[1], ProfVsAngle_TID, track->TID_fullHits, thePointNumber, false, commonHistos[2], Monitors_TID);
		use_onstrip=3;
		FillHistos(HistSoN_TID_onstrip[2], ProfVsAngle_TID, track->TID_fullHits, thePointNumber, false, commonHistos[2], Monitors_TID);
		use_onstrip=1;
	  }

      // TEC
	  if(part==TEC || part==All) 
	  {
	    if(track->TEC_hits.size()) FillHistos(HistSoN_TEC, ProfVsAngle_TEC, track->TEC_hits, thePointNumber, false, commonHistos[3], Monitors_TEC);
	    else if(track->TEC_fullHits.size()) FillHistos(HistSoN_TEC, ProfVsAngle_TEC, track->TEC_fullHits, thePointNumber, false, commonHistos[3], Monitors_TEC);
      }
	  

    } // End of loop over tracks
  } // End of loop over events
 
  

  // Fit of histos with Landau and store results
  std::cout << "Starting to perform fits " << std::endl;
  
  if((part==1 || part==0) && !use_onstrip) {
    FitHistos(HistSoN_TIB, "TIB_output_DecoMode.root", commonHistos[0], Monitors_TIB);
    FitProfiles(ProfVsAngle_TIB, "TIB_output_DecoMode_prof.root");
  }
  if((part==1 || part==0) && use_onstrip!=0) {
    FitHistos(HistSoN_TIB_onstrip[0], "TIB_output_DecoMode_onstrip_1.root", commonHistos[0], Monitors_TIB);
    FitHistos(HistSoN_TIB_onstrip[1], "TIB_output_DecoMode_onstrip_2.root", commonHistos[0], Monitors_TIB);
    FitHistos(HistSoN_TIB_onstrip[2], "TIB_output_DecoMode_onstrip_3.root", commonHistos[0], Monitors_TIB);
  }
  
  if((part==2 || part==0)) FitHistos(HistSoN_TOB, "TOB_output_DecoMode.root", commonHistos[1], Monitors_TOB);
  
  if((part==3 || part==0) && !use_onstrip) FitHistos(HistSoN_TID, "TID_output_DecoMode.root", commonHistos[2], Monitors_TID);
  if((part==3 || part==0) && use_onstrip!=0) {
    FitHistos(HistSoN_TID_onstrip[0], "TID_output_DecoMode_onstrip_1.root", commonHistos[2], Monitors_TID);
    FitHistos(HistSoN_TID_onstrip[1], "TID_output_DecoMode_onstrip_2.root", commonHistos[2], Monitors_TID);
    FitHistos(HistSoN_TID_onstrip[2], "TID_output_DecoMode_onstrip_3.root", commonHistos[2], Monitors_TID);
  }
  if((part==4 || part==0)) FitHistos(HistSoN_TEC, "TEC_output_DecoMode.root", commonHistos[3], Monitors_TEC);
  
  

  // If evtid selection used, print corresponding timestamp ranges
  if(!usetimestamp) VSmaker.printComputedSteps();

}

//---------------------------------------------------------------------------------------------------------


void ClusterWidthAnalysisTreeMaker::FillHitInfo(std::map<ULong64_t , std::vector<TH1F*> > &HistSoN, std::map<ULong64_t , std::vector<TProfile*> > &ProfVsAngle,
 TreeHit *hit, int thePointNumber, bool sensors, std::vector< TH1F* > commonHistos)
{
  FillHitInfo(HistSoN, ProfVsAngle, hit, thePointNumber, sensors, commonHistos, 0, 0, 0);
}

void ClusterWidthAnalysisTreeMaker::FillHitInfo(std::map<ULong64_t , std::vector<TH1F*> > &HistSoN, std::map<ULong64_t , std::vector<TProfile*> > &ProfVsAngle,
 TreeFullHit *hit, int thePointNumber, bool sensors, std::vector< TH1F* > commonHistos)
{
  FillHitInfo(HistSoN, ProfVsAngle, dynamic_cast< TreeHit* >(hit), thePointNumber, sensors, commonHistos, hit->barycenter, hit->seed, hit->seedChargeAngleCorr);
}
  
void ClusterWidthAnalysisTreeMaker::FillHitInfo(std::map<ULong64_t , std::vector<TH1F*> > &HistSoN, std::map<ULong64_t , std::vector<TProfile*> > &ProfVsAngle,
 TreeHit *hit, int thePointNumber, bool sensors, std::vector< TH1F* > commonHistos, float
barycenter, float seed, float seedChargeAngleCorr)
{
  if(thePointNumber<0) return;

  float angle = hit->angle;
  if(angle>TMath::Pi()/2) angle -= TMath::Pi();
  if(angle<-TMath::Pi()/2) angle += TMath::Pi();
  //if(fabs(angle)<0. || fabs(angle)>0.2) return;
  //if(angle>0 && angle<0.15) return;

  std::map<ULong64_t , std::vector<TH1F*> >::iterator iter;
  std::vector<TH1F*> detidvector;
  std::map<ULong64_t , std::vector<TProfile*> >::iterator iterProf;
  std::vector<TProfile*> profvector;

  int firstsensor = hit->tsosy<0 ? 1 : 2;
  
  ULong64_t detid = hit->detId;
  if(sensors)  detid = detid*10+firstsensor;
  //std::cout << "Detid " << detid << " "<< hit->detId <<" "<<firstsensor<< std::endl;
 
  iter = HistSoN.find(detid);
  if(iter == HistSoN.end() ){
	//std::cout << "Detid " << detIds[i] << " not found yet, creating vector " << std::endl;
	detidvector.clear();
	profvector.clear();
	TString histoname;
	histoname.Form("DetID_%llu",detid);
	TString profname;
	profname.Form("DetID_%llu_prof",detid);	
	if(use_onstrip!=0) histoname.Append(Form("_%i",use_onstrip)); 
	for(int j=0; j< VSmaker.getNVoltage(); j++){
	  std::string s;
	  std::stringstream out;
	  out << j;
	  s = out.str();
	  TString thestr = histoname+"_"+s;
	  detidvector.push_back( new TH1F(thestr.Data(), thestr.Data() , 20, 0, 20)  ); 
	  TString theprofstr = profname+"_"+s;
	  profvector.push_back( new TProfile(theprofstr.Data(), theprofstr.Data() , 60, -3, 3)  ); 
	}

	HistSoN.insert(std::pair<ULong64_t,std::vector<TH1F*> >(detid,detidvector));
	ProfVsAngle.insert(std::pair<ULong64_t,std::vector<TProfile*> >(detid,profvector));
    iter = HistSoN.find(detid);
  }

  if( iter == HistSoN.end() ){
	std::cout << "Error could not find Histogram with DetID: " << detid << std::endl;
	return;
  }

  if(commonHistos.size()>2) commonHistos[2]->Fill(hit->tsosx-hit->clusx);
  // before using angle info, check that traj measurement close to cluster position
  if(use_angle && abs(hit->tsosx-hit->clusx)>0.001) return; // a verifier
  // cut choisi pour TIB


  bool usehit=false;
  float angle_deg = hit->angle/3.141592*180;
  if(angle_deg>90) angle_deg = angle_deg-180;
  if(angle_deg<-90) angle_deg = angle_deg+180;
  if(commonHistos.size()>1) commonHistos[1]->Fill(angle_deg);
  angle_deg = TMath::Abs( angle_deg );

  if(use_angle==0) usehit=true; // default
  if(use_angle==1 && angle_deg<10) usehit=true;
  if(use_angle==2 && angle_deg>20 && angle_deg<40) usehit=true;
  if(use_angle==3 && angle_deg>40 && angle_deg<60) usehit=true;
  if(use_angle==4 && angle_deg>60 && angle_deg<80) usehit=true;
  if(commonHistos.size()>1) commonHistos[1]->Fill(angle_deg);

  /*if(use_width==1 && hit->width!=1 && usehit) usehit=false;
  if(use_width==2 && hit->width!=2 && usehit) usehit=false;
  if(use_width==3 && hit->width!=3 && usehit) usehit=false;
  if(use_width==4 && hit->width!=4 && usehit) usehit=false;
  if(use_width==5 && hit->width!=5 && usehit) usehit=false;*/


  // for TIB
  //float pitch = 1.;
  /*if(nstrips) {
    if(nstrips[i]==512) pitch = 0.012;
    if(nstrips[i]==768) pitch = 0.008;
  }*/
  //float offset = -3.072;
  //float cluspos = (clusx[i] - offset)/pitch -0.5; // prendre clusx ou tsosx ? 
  float cluspos = 0;
  if(barycenter) cluspos = barycenter -0.5; // ou barycenter
  cluspos -= floor(cluspos); if(cluspos>0.5) cluspos-=1;
  //if(seed) cluspos-=seed;

  if(commonHistos.size()>3 && seed && usehit) commonHistos[3]->Fill(cluspos);

  int clus_onStrip = 0;
  if(abs(cluspos)<0.125) clus_onStrip = 1;
  if(abs(cluspos)>0.125 && abs(cluspos)<0.375) clus_onStrip = (cluspos>0) ? 2 : -2;
  if(abs(cluspos)>0.375 && abs(cluspos)<0.625) clus_onStrip = (cluspos>0) ? 3 : -3;
  if(abs(cluspos)>0.625 && abs(cluspos)<0.875) clus_onStrip = (cluspos>0) ? 4 : -4;


  if(use_onstrip==1 && clus_onStrip!=1 && usehit) usehit=false;
  if(use_onstrip==2 && abs(clus_onStrip)!=2 && usehit) usehit=false;
  if(use_onstrip==3 && abs(clus_onStrip)!=3 && usehit) usehit=false;
  if(use_onstrip==4 && clus_onStrip!=4 && usehit) usehit=false;
  if(use_onstrip==-4 && clus_onStrip!=-4 && usehit) usehit=false;


  if(hit->chargeAngleCorr > 0 && usehit)
    iter->second.at(thePointNumber)->Fill(hit->width);

  iterProf = ProfVsAngle.find(detid);
  if( iterProf != ProfVsAngle.end() )
    iterProf->second.at(thePointNumber)->Fill(tan(angle), hit->width);
}

void ClusterWidthAnalysisTreeMaker::FillHistos(std::map<ULong64_t , std::vector<TH1F*> > &HistSoN, std::map<ULong64_t , std::vector<TProfile*> > &ProfVsAngle,
 std::vector< TreeHit > &Hits, int thePointNumber, bool sensors, std::vector< TH1F* > commonHistos, std::map<ULong64_t, TProfile* > Monitors)
{

  std::map<ULong64_t, TProfile* >::iterator itMon;
  for(unsigned int i=0; i<Hits.size(); i++)
  {
    FillHitInfo(HistSoN, ProfVsAngle, &(Hits[i]), thePointNumber, sensors, commonHistos);
	for(itMon=Monitors.begin(); itMon!=Monitors.end(); itMon++)
	  if(itMon->first==Hits[i].detId) itMon->second->Fill(event->ev_timestamp - t_monitor_start, Hits[i].width);
  }

}

void ClusterWidthAnalysisTreeMaker::FillHistos(std::map<ULong64_t , std::vector<TH1F*> > &HistSoN, std::map<ULong64_t , std::vector<TProfile*> > &ProfVsAngle,
 std::vector< TreeFullHit > &Hits, int thePointNumber, bool sensors, std::vector< TH1F* > commonHistos, std::map<ULong64_t, TProfile* > Monitors)
{

  std::map<ULong64_t, TProfile* >::iterator itMon;
  for(unsigned int i=0; i<Hits.size(); i++)
  {
    FillHitInfo(HistSoN, ProfVsAngle, &(Hits[i]), thePointNumber, sensors, commonHistos);
	for(itMon=Monitors.begin(); itMon!=Monitors.end(); itMon++)
	  if(itMon->first==Hits[i].detId) itMon->second->Fill(event->ev_timestamp - t_monitor_start, Hits[i].width);
  }

}

//---------------------------------------------------------------------------------------------------------


void ClusterWidthAnalysisTreeMaker::FitHistos(std::map<ULong64_t , std::vector<TH1F*> > &HistSoN, string output_file, 
 std::vector< TH1F* > commonHistos, std::map<ULong64_t, TProfile* > Monitors){

  TFile * myFile = new TFile(output_file.c_str(), "recreate");

  ULong64_t detid;
  double voltage;
  double errvoltage;
  double Mean;
  double errMean;
  double RMS;
  double errRMS;
  int index;
  int nhits;
  TTree *tree = new TTree("T", "summary information");

  tree->Branch("DetID",&detid, "DetID/l");
  tree->Branch("Voltage",&voltage,"Voltage/D");
  tree->Branch("Index",&index,"Index/I");
  tree->Branch("errVoltage",&errvoltage,"errVoltage/D");
  tree->Branch("Mean",&Mean,"Mean/D");
  tree->Branch("errMean",&errMean,"errMean/D");
  tree->Branch("RMS",&RMS,"RMS/D");
  tree->Branch("errRMS",&errRMS,"errRMS/D");
  tree->Branch("Nhits",&nhits,"Nhits/I");


  //TCanvas* c1 = new TCanvas();
  TH1F* hNhits = new TH1F("hNhits", "hNhits", 1000, 0,1000); // N hits per module
  
  unsigned int nfitrm=0;

  for(std::map<ULong64_t , std::vector<TH1F*> >::iterator iter = HistSoN.begin(); iter != HistSoN.end(); ++iter){
    
	unsigned int i=0; // voltage index    
    std::set< int >::iterator itVolt;
	std::set< int > Voltage = VSmaker.getVoltageList();
    for( itVolt=Voltage.begin(); itVolt!=Voltage.end(); itVolt++){
      
      //std::cout<<"going through the measurement: " << i << std::endl;
            
      TString thestring;
      thestring.Form("DetID_%llu_%u",iter->first,i);
 	  
	  
      //std::cout << "searching for " << thestring.Data() << std::endl;
      //TH1F*  SoNHisto= (TH1F*)gROOT->FindObject( thestring.Data() );
	  
	  if(i>=iter->second.size()) 
       { std::cout<<" Wrong number of voltage steps. "<<std::endl; i++; continue;}
 	  TH1F*  Histo = iter->second[i];
	  
	  if(!Histo) 
       { std::cout<<" Histo "<<thestring.Data()<<"_"<<i<<" not found."<<std::endl; i++; continue;}
 
      if(Histo->GetEntries()) hNhits->Fill(Histo->Integral());
	  
	  if(Histo->Integral()<20) //0.1
	   { //std::cout<<" Not enought entries for histo "<<thestring.Data()<<std::endl;
	    i++; continue;}
 
	  
	  detid = iter->first;

	  bool rmfit=false;

      if( rmfit || 
      // TIB modules
          // TIB - 1.4.2.5
      detid==369121605 || detid==369121606 || detid==369121614 || 
      detid==369121613 || detid==369121610 || detid==369121609 ||
          // TIB - 1.2.2.1
      detid==369121390 || detid==369121382 || detid==369121386 || 
      detid==369121385 || detid==369121389 || detid==369121381 ||
          // others in TIB  
      detid==369121437 || detid==369142077 || detid==369121722 || 
      detid==369125534 || detid==369137018 || detid==369121689 ||
      detid==369121765 || detid==369137045 || detid==369169740 ||
      detid==369121689 ||
      // TOB modules 
	      // TOB + 4.3.3.8
      detid/10==436281512 || detid/10==436281528 || detid/10==436281508 ||
      detid/10==436281524 || detid/10==436281520 || detid/10==436281516 ||
          // others in TOB  
      detid/10==436228249 || detid/10==436232694 || detid/10==436228805 ||
      detid/10==436244722 || detid/10==436245110 || detid/10==436249546 ||
      detid/10==436310808 || detid/10==436312136 || detid/10==436315600 ||
	      // without 'sensors' option 
      detid==436281512 || detid==436281528 || detid==436281508 ||
      detid==436281524 || detid==436281520 || detid==436281516 ||
      detid==436228249 || detid==436232694 || detid==436228805 ||
      detid==436244722 || detid==436245110 || detid==436249546 ||
      detid==436310808 || detid==436312136 || detid==436315600 || 
      // TID modules
      detid==402664070 || detid==402664110 ||
	  // TEC modules in small scans
      detid==470148196 || detid==470148200 || detid==470148204 ||
      detid==470148228 || detid==470148232 || detid==470148236 ||
      detid==470148240 || detid==470148261 || detid==470148262 ||
	  detid==470148265 || detid==470148266 || detid==470148292 ||
	  detid==470148296 || detid==470148300 || detid==470148304 ||
	  detid==470148324 || detid==470148328 || detid==470148332 ||
	  detid==470148336 || detid==470148340 )  { 
	    Histo->Write();
        std::cout << " Saving histo : " << thestring.Data() << std::endl;
      }  


	  if(rmfit) {nfitrm++; i++; continue;}

          int subdet = ((detid>>25)&0x7);
          int TECgeom=0;
          if(subdet==6) TECgeom = ((detid>>5)&0x7);

      // save values
	  detid = iter->first;
	  voltage  = *itVolt;
	  index = i;
	  errvoltage = 2 ;
	  Mean = Histo->GetMean();
	  errMean = Histo->GetMeanError();
	  RMS = Histo->GetRMS();
	  errRMS = Histo->GetRMSError();
	  nhits = (int) Histo->Integral();
	  tree->Fill();
	    
	  i++;

    }  

  }
  
  tree->Write();
  hNhits->Write();
  

  for(unsigned int ih=0; ih<commonHistos.size(); ih++) commonHistos[ih]->Write();

  std::map<ULong64_t, TProfile* >::iterator itMon;
  for(itMon=Monitors.begin(); itMon!=Monitors.end(); itMon++)
  {
    itMon->second->GetXaxis()->SetTimeDisplay(1);
	itMon->second->GetXaxis()->SetTimeFormat("%H:%M");
	itMon->second->GetXaxis()->SetTimeOffset(t_monitor_start);
	itMon->second->Write();
  }
  
  //// If you want to store all the individual detId histograms uncomments this line !!!!
  //myFile->Write();
  myFile->Close();

}


void ClusterWidthAnalysisTreeMaker::FitProfiles(std::map<ULong64_t , std::vector<TProfile*> > &ProfVsAngle, string output_file){

  TFile * myFile = new TFile(output_file.c_str(), "recreate");

  ULong64_t detid;
  double voltage;
  double errvoltage;
  double Slope;
  double errSlope;
  double Origin;
  double errOrigin;
  double Chi2;
  int index;
  TTree *tree = new TTree("T", "summary information");

  tree->Branch("DetID",&detid, "DetID/l");
  tree->Branch("Voltage",&voltage,"Voltage/D");
  tree->Branch("Index",&index,"Index/I");
  tree->Branch("errVoltage",&errvoltage,"errVoltage/D");
  tree->Branch("Slope",&Slope,"Slope/D");
  tree->Branch("errSlope",&errSlope,"errSlope/D");
  tree->Branch("Origin",&Origin,"Origin/D");
  tree->Branch("errOrigin",&errOrigin,"errOrigin/D");
  tree->Branch("Chi2",&Chi2,"Chi2/D");


  //TCanvas* c1 = new TCanvas();
  TH1F* hChi2 = new TH1F("hChi2", "hChi2", 100, 0, 100);
  
  unsigned int nfitrm=0;

  for(std::map<ULong64_t , std::vector<TProfile*> >::iterator iter = ProfVsAngle.begin(); iter != ProfVsAngle.end(); ++iter){
    
	unsigned int i=0; // voltage index    
    std::set< int >::iterator itVolt;
	std::set< int > Voltage = VSmaker.getVoltageList();
    for( itVolt=Voltage.begin(); itVolt!=Voltage.end(); itVolt++){
      
      //std::cout<<"going through the measurement: " << i << std::endl;
            
      TString thestring;
      thestring.Form("DetID_%llu_prof_%u",iter->first,i);
 	  	  
	  if(i>=iter->second.size()) 
       { std::cout<<" Wrong number of voltage steps. "<<std::endl; i++; continue;}
 	  TProfile*  Prof = iter->second[i];
	  
	  if(!Prof) 
       { std::cout<<" Profile "<<thestring.Data()<<"_"<<i<<" not found."<<std::endl; i++; continue;}
 
      //if(Histo->GetEntries()) hNhits->Fill(Histo->Integral());
	  
	  /*if(Histo->Integral()<20) //0.1
	   { //std::cout<<" Not enought entries for histo "<<thestring.Data()<<std::endl;
	    i++; continue;}*/
 
	  
	  detid = iter->first;

	  bool rmfit=false;
	  
	  TF1* pol = new TF1("pol", "pol1", -3, 3);
	  pol->SetRange(-1,0);
	  Prof->Fit("pol", "qr");
	  double chi2 = pol->GetChisquare()/pol->GetNDF();
	  hChi2->Fill(chi2);
	  if(chi2>10) rmfit=true;
	  

      if( rmfit || 
      // TIB modules
          // TIB - 1.4.2.5
      detid==369121605 || detid==369121606 || detid==369121614 || 
      detid==369121613 || detid==369121610 || detid==369121609 ||
          // TIB - 1.2.2.1
      detid==369121390 || detid==369121382 || detid==369121386 || 
      detid==369121385 || detid==369121389 || detid==369121381 ||
          // others in TIB  
      detid==369121437 || detid==369142077 || detid==369121722 || 
      detid==369125534 || detid==369137018 || detid==369121689 ||
      detid==369121765 || detid==369137045 || detid==369169740 ||
      detid==369121689 ||
      // TOB modules 
	      // TOB + 4.3.3.8
      detid/10==436281512 || detid/10==436281528 || detid/10==436281508 ||
      detid/10==436281524 || detid/10==436281520 || detid/10==436281516 ||
          // others in TOB  
      detid/10==436228249 || detid/10==436232694 || detid/10==436228805 ||
      detid/10==436244722 || detid/10==436245110 || detid/10==436249546 ||
      detid/10==436310808 || detid/10==436312136 || detid/10==436315600 ||
	      // without 'sensors' option 
      detid==436281512 || detid==436281528 || detid==436281508 ||
      detid==436281524 || detid==436281520 || detid==436281516 ||
      detid==436228249 || detid==436232694 || detid==436228805 ||
      detid==436244722 || detid==436245110 || detid==436249546 ||
      detid==436310808 || detid==436312136 || detid==436315600 || 
      // TID modules
      detid==402664070 || detid==402664110 ||
	  // TEC modules in small scans
      detid==470148196 || detid==470148200 || detid==470148204 ||
      detid==470148228 || detid==470148232 || detid==470148236 ||
      detid==470148240 || detid==470148261 || detid==470148262 ||
	  detid==470148265 || detid==470148266 || detid==470148292 ||
	  detid==470148296 || detid==470148300 || detid==470148304 ||
	  detid==470148324 || detid==470148328 || detid==470148332 ||
	  detid==470148336 || detid==470148340 )  { 
	    Prof->Write();
        std::cout << " Saving histo : " << thestring.Data() << std::endl;
      }  


	  if(rmfit) {nfitrm++; i++; continue;}

          int subdet = ((detid>>25)&0x7);
          int TECgeom=0;
          if(subdet==6) TECgeom = ((detid>>5)&0x7);

      // save values
	  detid = iter->first;
	  voltage  = *itVolt;
	  index = i;
	  errvoltage = 2 ;
	  Slope = pol->GetParameter(1);
	  errSlope = pol->GetParError(1);
	  Origin = pol->GetParameter(0);
	  errOrigin = pol->GetParError(0);
	  Chi2 = chi2;
	  tree->Fill();
	    
	  i++;

    }  

  }
  
  tree->Write();
  //hNhits->Write();

  
  //// If you want to store all the individual detId histograms uncomments this line !!!!
  //myFile->Write();
  myFile->Close();

}
