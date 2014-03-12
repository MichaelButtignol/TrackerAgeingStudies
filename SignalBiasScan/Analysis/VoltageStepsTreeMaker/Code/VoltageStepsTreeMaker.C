#define VoltageStepsTreeMaker_cxx
#include "VoltageStepsTreeMaker.h"
#include <iostream>
#include <time.h>
#include <sstream>
#include <set>
#include <utility>
#include <fstream>
#include <cmath>

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
#include <exception>


// 1 = TIB // 2 = TOB // 3 = TID // 4 = TEC
// ****************************************

enum SubDet { All, TIB, TOB, TID, TEC};

using namespace std;

Double_t fitLandauGauss(Double_t *x, Double_t *par){

  Double_t value =  (1-par[3])*TMath::Landau(x[0], par[0], par[1]) + par[3]*TMath::Gaus(x[0], par[0], par[2]);
  return value;
  
}

time_t convertTimestamp( std::string str){

  std::stringstream ss;
  tm time;
  int year=-1;
  int month=-1;
  
  ss.clear(); ss << str.substr(0, 4); ss >> year;
  time.tm_year = year - 1900;
  ss.clear(); ss << str.substr(4, 2); ss >> month;
  time.tm_mon = month - 1;
  ss.clear(); ss << str.substr(6, 2); ss >> time.tm_mday;    
  ss.clear(); ss << str.substr(8, 2); ss >> time.tm_hour;    
  ss.clear(); ss << str.substr(10, 2); ss >> time.tm_min;    
  ss.clear(); ss << str.substr(12, 2); ss >> time.tm_sec;
  //time.tm_isdst=1;
  
  if(year < 2000 || year > 2020) std::cout<<" Wrong year format : "<<year<<std::endl;
  if(time.tm_sec < 0 || time.tm_sec > 61) std::cout<<" Wrong sec format : "<<year<<std::endl;
  //cout<<" timestamp "<<year<<" "<<time.tm_mon<<" "<<time.tm_mday<<" "<<time.tm_hour<<" "<<time.tm_min<<" "<<time.tm_sec<<std::endl;

  time_t out_time = mktime( &time );
  ctime(&out_time); // Have to let this line, otherwise 1 hour shift with gcc compared to AClic !!
  //2 hours shift, due to daylight saving time + shift of local time with UTC ?
  // only 1 hour shift starting from nov 2012
  if(year>2012 || (year==2012 && month>=10)) out_time+=3600;
  else out_time+=7200;
  
  return out_time;
  
}

// Main method //
//-------------//
void VoltageStepsTreeMaker::Loop()
{
  
  std::cout<<"Initializing"<<std::endl;
  
  std::string output_file;
  
  // Read definition of voltage steps
  std::cout << "Setting Deco Voltage vector : " << std::endl;
  if(usetimestamp) std::cout << " using timestamp step definition " << std::endl;
  else std::cout << " using event number step definition" << std::endl;
  
  //if(usetimestamp) setVoltageSteps_timestamp(stepsfile.c_str());
  //else setVoltageSteps_evtnumber(stepsfile.c_str());
  setVoltageSteps(stepsfile.c_str(), usetimestamp);
  std::cout << " done" << std::endl;
  
  std::set< int >::iterator itVolt;
  for( itVolt=Voltage.begin(); itVolt!=Voltage.end(); itVolt++){
    VoltageSteps_timestart[*itVolt] = 2147483647;//int max
    VoltageSteps_timeend[*itVolt] = 0;
  }

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
	if(usetimestamp) theVoltage=getVoltageSteps_timestamp(event->ev_timestamp);
	else theVoltage=getVoltageSteps_evtnumber(event->run_nr, event->ev_nr, event->ev_timestamp);
	if(theVoltage<0) continue; // skip event if not on a voltage step

	// Define index of Voltage for storage
	int ipt=0;
	int thePointNumber=-1;
	for( itVolt=Voltage.begin(); itVolt!=Voltage.end(); itVolt++){
      if(theVoltage==(*itVolt)) thePointNumber=ipt;
      ipt++;
	}

	// check index
	if(thePointNumber < 0 && thePointNumber >= (int)Voltage.size()){
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
	    if(track->TIB_hits.size()) FillHistSoN(HistSoN_TIB, track->TIB_hits, thePointNumber, false, commonHistos[0], Monitors_TIB);
	    else if(track->TIB_fullHits.size()) FillHistSoN(HistSoN_TIB, track->TIB_fullHits, thePointNumber, false, commonHistos[0], Monitors_TIB);
      }

	  if((part==TIB || part==All) && use_onstrip!=0)
	  {
		use_onstrip=1;
		FillHistSoN(HistSoN_TIB_onstrip[0], track->TIB_fullHits, thePointNumber, false, commonHistos[0], Monitors_TIB);
		use_onstrip=2;
		FillHistSoN(HistSoN_TIB_onstrip[1], track->TIB_fullHits, thePointNumber, false, commonHistos[0], Monitors_TIB);
		use_onstrip=3;
		FillHistSoN(HistSoN_TIB_onstrip[2], track->TIB_fullHits, thePointNumber, false, commonHistos[0], Monitors_TIB);
		use_onstrip=1;
	  }

      // TOB
	  if(part==TOB || part==All)
	  {
	    if(track->TOB_hits.size()) FillHistSoN(HistSoN_TOB, track->TOB_hits, thePointNumber, true, commonHistos[1], Monitors_TOB);
	    else if(track->TOB_fullHits.size()) FillHistSoN(HistSoN_TOB, track->TOB_fullHits, thePointNumber, true, commonHistos[1], Monitors_TOB);
      }

      // TID
	  if((part==TID || part==All) && !use_onstrip)
	  {
	    if(track->TID_hits.size()) FillHistSoN(HistSoN_TID, track->TID_hits, thePointNumber, false, commonHistos[2], Monitors_TID);
	    else if(track->TID_fullHits.size()) FillHistSoN(HistSoN_TID, track->TID_fullHits, thePointNumber, false, commonHistos[2], Monitors_TID);	  
	  }
	  
	  if((part==TID || part==All) && use_onstrip!=0)
	  {
		use_onstrip=1;
		FillHistSoN(HistSoN_TID_onstrip[0], track->TID_fullHits, thePointNumber, false, commonHistos[2], Monitors_TID);
		use_onstrip=2;
		FillHistSoN(HistSoN_TID_onstrip[1], track->TID_fullHits, thePointNumber, false, commonHistos[2], Monitors_TID);
		use_onstrip=3;
		FillHistSoN(HistSoN_TID_onstrip[2], track->TID_fullHits, thePointNumber, false, commonHistos[2], Monitors_TID);
		use_onstrip=1;
	  }

      // TEC
	  if(part==TEC || part==All) 
	  {
	    if(track->TEC_hits.size()) FillHistSoN(HistSoN_TEC, track->TEC_hits, thePointNumber, false, commonHistos[3], Monitors_TEC);
	    else if(track->TEC_fullHits.size()) FillHistSoN(HistSoN_TEC, track->TEC_fullHits, thePointNumber, false, commonHistos[3], Monitors_TEC);
      }
	  

    } // End of loop over tracks
  } // End of loop over events
 
  

  // Fit of histos with Landau and store results
  std::cout << "Starting to perform fits " << std::endl;
  
  if((part==1 || part==0) && !use_onstrip) FitHistos(HistSoN_TIB, "TIB_output_DecoMode.root", commonHistos[0], Monitors_TIB);
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
  if(!usetimestamp)
  std::cout << "Time definition of steps : " << std::endl;
  if(!usetimestamp)
  for( itVolt=Voltage.begin(); itVolt!=Voltage.end(); itVolt++){
    std::cout<<(*itVolt)<<" "<<VoltageSteps_timestart[*itVolt]<<" "<<VoltageSteps_timeend[*itVolt]<<std::endl;
    time_t tt = VoltageSteps_timestart[*itVolt];
	std::cout<<"   start : "<<ctime(&tt);
    tt = VoltageSteps_timeend[*itVolt];
	std::cout<<"     end : "<<ctime(&tt);
  }

}

//---------------------------------------------------------------------------------------------------------


void VoltageStepsTreeMaker::FillHitInfo(std::map<ULong64_t , std::vector<TH1F*> > &HistSoN, TreeHit *hit, int thePointNumber, bool sensors, std::vector< TH1F* > commonHistos)
{
  FillHitInfo(HistSoN, hit, thePointNumber, sensors, commonHistos, 0, 0, 0);
}

void VoltageStepsTreeMaker::FillHitInfo(std::map<ULong64_t , std::vector<TH1F*> > &HistSoN, TreeFullHit *hit, int thePointNumber, bool sensors, std::vector< TH1F* > commonHistos)
{
  FillHitInfo(HistSoN, dynamic_cast< TreeHit* >(hit), thePointNumber, sensors, commonHistos, hit->barycenter, hit->seed, hit->seedChargeAngleCorr);
}
  
void VoltageStepsTreeMaker::FillHitInfo(std::map<ULong64_t , std::vector<TH1F*> > &HistSoN, TreeHit *hit, int thePointNumber, bool sensors, std::vector< TH1F* > commonHistos, float
barycenter, float seed, float seedChargeAngleCorr)
{
  if(thePointNumber<0) return;

  std::map<ULong64_t , std::vector<TH1F*> >::iterator iter;
  std::vector<TH1F*> detidvector;

  int firstsensor = hit->tsosy<0 ? 1 : 2;
  
  ULong64_t detid = hit->detId;
  if(sensors)  detid = detid*10+firstsensor;
  //std::cout << "Detid " << detid << " "<< hit->detId <<" "<<firstsensor<< std::endl;
 
  iter = HistSoN.find(detid);
  if(iter == HistSoN.end() ){
	//std::cout << "Detid " << detIds[i] << " not found yet, creating vector " << std::endl;
	detidvector.clear();
	TString histoname;
	histoname.Form("DetID_%llu",detid);
	if(use_onstrip!=0) histoname.Append(Form("_%i",use_onstrip)); 
	for(int j=0; j<(int)Voltage.size(); j++){
	  std::string s;
	  std::stringstream out;
	  out << j;
	  s = out.str();
	  TString thestr = histoname+"_"+s;
	  detidvector.push_back( new TH1F(thestr.Data(), thestr.Data() , 90, 0, 300)  ); // 30, 0, 100 
	  // 60, 0, 200 // not far enought for TOB -> 90 0 300
//std::cout << thestr << std::endl;
	}

	HistSoN.insert(std::pair<ULong64_t,std::vector<TH1F*> >(detid,detidvector));
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

  if(use_width==1 && hit->width!=1 && usehit) usehit=false;
  if(use_width==2 && hit->width!=2 && usehit) usehit=false;
  if(use_width==3 && hit->width!=3 && usehit) usehit=false;
  if(use_width==4 && hit->width!=4 && usehit) usehit=false;
  if(use_width==5 && hit->width!=5 && usehit) usehit=false;


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
    iter->second.at(thePointNumber)->Fill(hit->chargeAngleCorr);


}

void VoltageStepsTreeMaker::FillHistSoN(std::map<ULong64_t , std::vector<TH1F*> > &HistSoN, std::vector< TreeHit > &Hits, 
 int thePointNumber, bool sensors, std::vector< TH1F* > commonHistos, std::map<ULong64_t, TProfile* > Monitors)
{

  std::map<ULong64_t, TProfile* >::iterator itMon;
  for(unsigned int i=0; i<Hits.size(); i++)
  {
    FillHitInfo(HistSoN, &(Hits[i]), thePointNumber, sensors, commonHistos);
	for(itMon=Monitors.begin(); itMon!=Monitors.end(); itMon++)
	  if(itMon->first==Hits[i].detId) itMon->second->Fill(event->ev_timestamp - t_monitor_start, Hits[i].chargeAngleCorr);
  }

}

void VoltageStepsTreeMaker::FillHistSoN(std::map<ULong64_t , std::vector<TH1F*> > &HistSoN, std::vector< TreeFullHit > &Hits, 
 int thePointNumber, bool sensors, std::vector< TH1F* > commonHistos, std::map<ULong64_t, TProfile* > Monitors)
{

  std::map<ULong64_t, TProfile* >::iterator itMon;
  for(unsigned int i=0; i<Hits.size(); i++)
  {
    FillHitInfo(HistSoN, &(Hits[i]), thePointNumber, sensors, commonHistos);
	for(itMon=Monitors.begin(); itMon!=Monitors.end(); itMon++)
	  if(itMon->first==Hits[i].detId) itMon->second->Fill(event->ev_timestamp - t_monitor_start, Hits[i].chargeAngleCorr);
  }

}

//---------------------------------------------------------------------------------------------------------


void VoltageStepsTreeMaker::FitHistos(std::map<ULong64_t , std::vector<TH1F*> > &HistSoN, string output_file, 
 std::vector< TH1F* > commonHistos, std::map<ULong64_t, TProfile* > Monitors){

  TFile * myFile = new TFile(output_file.c_str(), "recreate");

  ULong64_t detid;
  double voltage;
  double errvoltage;
  double MPV;
  double errMPV;
  double Width;
  double errWidth;
  int index;
  double chi2overndf;
  int nhits;
  double pedestal_fraction;
  TTree *tree = new TTree("T", "summary information");

  tree->Branch("DetID",&detid, "DetID/l");
  tree->Branch("Voltage",&voltage,"Voltage/D");
  tree->Branch("Index",&index,"Index/I");
  tree->Branch("errVoltage",&errvoltage,"errVoltage/D");
  tree->Branch("MPV",&MPV,"MPV/D");
  tree->Branch("errMPV",&errMPV,"errMPV/D");
  tree->Branch("Width",&Width,"Width/D");
  tree->Branch("errWidth",&errWidth,"errWidth/D");
  tree->Branch("Chi2OverNdf",&chi2overndf,"Chi2OverNdf/D");
  tree->Branch("Nhits",&nhits,"Nhits/I");
  tree->Branch("PedFraction",&pedestal_fraction,"PedFraction/D");


  //TCanvas* c1 = new TCanvas();
  TH1F* hNhits = new TH1F("hNhits", "hNhits", 1000, 0,1000); // N hits per module
  TH1F* hChi2OverNDF = new TH1F("hChi2OverNDF", "hChi2OverNDF", 100, 0,100);
  TH2F* hChi2OverNDFvsstep = new TH2F("hChi2OverNDFvsstep", "hChi2OverNDFvsstep", 20, 0, 20, 100, 0,100);
  
  unsigned int nbadfits=0;
  unsigned int nnegpar0=0;
  unsigned int nnegpar1=0;
  unsigned int nfitrm=0;

  for(std::map<ULong64_t , std::vector<TH1F*> >::iterator iter = HistSoN.begin(); iter != HistSoN.end(); ++iter){
    
	unsigned int i=0; // voltage index    
    std::set< int >::iterator itVolt;
    for( itVolt=Voltage.begin(); itVolt!=Voltage.end(); itVolt++){
      
      //std::cout<<"going through the measurement: " << i << std::endl;
      
      TF1 *fitFunc = new TF1("fitFunc","landau", 0, 200);
      fitFunc->SetLineColor(1);
      
      TString thestring;
      thestring.Form("DetID_%llu_%u",iter->first,i);
 	  
	  
      //std::cout << "searching for " << thestring.Data() << std::endl;
      //TH1F*  SoNHisto= (TH1F*)gROOT->FindObject( thestring.Data() );
	  
	  if(i>=iter->second.size()) 
       { std::cout<<" Wrong number of voltage steps. "<<std::endl; i++; continue;}
 	  TH1F*  SoNHisto = iter->second[i];
	  
	  if(!SoNHisto) 
       { std::cout<<" Histo "<<thestring.Data()<<"_"<<i<<" not found."<<std::endl; i++; continue;}
 
      if(SoNHisto->GetEntries()) hNhits->Fill(SoNHisto->Integral());
	  
	  if(SoNHisto->Integral()<20) //0.1
	   { //std::cout<<" Not enought entries for histo "<<thestring.Data()<<std::endl;
	    i++; continue;}
 
	  //fitFunc->SetParameter(0, SoNHisto->Integral("w")); //JLA
	  //fitFunc->SetParameter(1, 10); //JLA
	  //fitFunc->SetParLimits(0, 0, 1e+10);
	  //fitFunc->SetParLimits(1, 0, 100);
	  Int_t fitStatus = SoNHisto->Fit("fitFunc","q");
	  Float_t chi2overndf_temp = 999;
	  if(fitFunc->GetNDF()>0) chi2overndf_temp = fitFunc->GetChisquare()/fitFunc->GetNDF();
	  
	  if( fitStatus!=0 || fitFunc->GetParameter(0)<0 || fitFunc->GetParameter(1)<0 ||
	   fitFunc->GetParameter(1)>200 || chi2overndf_temp > 100 )
	  {
		//fitFunc->SetParameter(0, SoNHisto->Integral());
		//fitFunc->SetParLimits(0, 0, 1e+10);
    	//fitFunc->SetParLimits(1, 0, 100);
		//fitFunc->SetParameter(1, SoNHisto->GetMean()); 
		fitFunc->SetRange(0, SoNHisto->GetMean()+3*SoNHisto->GetRMS());
		fitStatus = SoNHisto->Fit("fitFunc","rq");
	  }

	  chi2overndf = 999;
	  if(fitFunc->GetNDF()>0) chi2overndf = fitFunc->GetChisquare()/fitFunc->GetNDF();
	  hChi2OverNDF->Fill(chi2overndf);
	  hChi2OverNDFvsstep->Fill( i, chi2overndf);
	  detid = iter->first;

	  bool rmfit=false;
	  if(fitStatus!=0 || fitFunc->GetParameter(0)<0 || fitFunc->GetParameter(1)<0 ||
	   fitFunc->GetParameter(1)>200 || chi2overndf > 100) rmfit=true;

      //if(MPV<0 || (MPV>25 && k<15) || errMPV>1 || detid==369169740)
      //if(errMPV>5)
      /*if(detid==369121605 || detid==369121606 || detid==369121614 || detid==369121613 ||
      detid==369121610 || detid==369121609 || detid==369121437 || detid==369142077 ||
      detid==369121722 || detid==369125534 || detid==369137018 || detid==369121689 ||
      detid==369121765 || detid==369137045 || detid==369169740)*/
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
	    SoNHisto->Write();
        std::cout << " Saving histo : " << thestring.Data() << std::endl;
	    if(fitStatus!=0) std::cout << " fit status : " << fitStatus << std::endl;
      }  


	  if(fitStatus!=0) nbadfits++;
	  if(fitFunc->GetParameter(0)<0) nnegpar0++;
	  if(fitFunc->GetParameter(1)<0) nnegpar1++;

	  if(rmfit) {nfitrm++; i++; continue;}

          int subdet = ((detid>>25)&0x7);
          int TECgeom=0;
          if(subdet==6) TECgeom = ((detid>>5)&0x7);

      // save values
	  detid = iter->first;
	  voltage  = *itVolt;
	  index = i;
	  errvoltage = 2 ;
	  MPV = fitFunc->GetParameter(1);
	  errMPV = fitFunc->GetParError(1);
	  Width = fitFunc->GetParameter(2);
	  errWidth = fitFunc->GetParError(2);
	  // chi2overndf already set
	  nhits = (int) SoNHisto->Integral();
      if(subdet==3 || subdet==4 || (subdet==6 && TECgeom<5))
         pedestal_fraction = SoNHisto->Integral(0, 16)/SoNHisto->Integral(17, 91); // fraction of noise below 53.3 ADC counts for thin sensors
      else 
         pedestal_fraction = SoNHisto->Integral(0, 27)/SoNHisto->Integral(28, 91); // fraction of noise below 90 ADC counts for thick sensors
	  tree->Fill();
	    
	  i++;

    }  

  }
  
  tree->Write();
  hNhits->Write();
  hChi2OverNDF->Write();
  hChi2OverNDFvsstep->Write();
  
  std::cout<<" N fits wo convergence : "<<nbadfits<<std::endl;
  std::cout<<" N fits with par0<0 : "<<nnegpar0<<std::endl;
  std::cout<<" N fits with par1<0 : "<<nnegpar1<<std::endl;
  std::cout<<" N fit results removed : "<<nfitrm<<std::endl;


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

//---------------------------------------------------------------------------------------------------------

void VoltageStepsTreeMaker::setVoltageSteps(const char* inputfile, bool usetimestamp){
  
  std::string line;
  ifstream fin(inputfile);
  
  int step=0;
  int step_previous=-1;
  time_t time_start=0;
  time_t time_end=0;
  std::string str_time;
  int evt=-1;
  int evt_start=-1;
  int evt_end=-1;
  std::string str_orbit;
  float voltage=0;
  float voltage_previous=0;
  int run=-1;
  int run_previous=-1;
  
  t_monitor_start=4294967295; // val max pour unsigned int sur 32 bits
  t_monitor_end=0;
  
  
  if(!fin.is_open())
    std::cerr<<" WARNING : file '"<<inputfile<<"' not found. Steps will not be set."<<std::endl;
  else
  {
    while( getline ( fin, line) )
	{
	  if(fin.eof()) continue;
	  std::stringstream ss(line);
	  ss >> step >> str_time >> evt >> str_orbit >> voltage >> run;
	  	  
	  if(step_previous!=step) 
	  {
	    evt_start = evt;
	    time_start = convertTimestamp( str_time );
	    if(time_start<0) std::cout<<" >>> WARNING : wrong timestamp conversion."<<std::endl;
		if(t_monitor_start>time_start) t_monitor_start=time_start;
	  std::cout<<step<<" t_start "<<time_start<<" "<<t_monitor_start<<std::endl;
	  }
	  else
	  {
	    evt_end = evt;
	    time_end = convertTimestamp( str_time );
	    if(time_end<0) std::cout<<" >>> WARNING : wrong timestamp conversion."<<std::endl;
		if(t_monitor_end<time_end) t_monitor_end=time_end;
	    
		if(fabs(voltage_previous-voltage)>0.1) 
		  std::cout<<" >>> ERROR : skip step "<<step<<". voltage_end != voltage_start !"<<std::endl;
		if(run_previous!=run) 
		  std::cout<<" >>> ERROR : skip step "<<step<<". run_end != run_start !"<<std::endl;
		if(evt_end < evt_start) 
		  std::cout<<" >>> ERROR : skip step "<<step<<". evt_end < evt_start !"<<std::endl;
		if(time_end < time_start) 
		  std::cout<<" >>> ERROR : skip step "<<step<<". time_end < time_start !"<<std::endl;
		
		// Store steps definitions
		if(run_previous==run && evt_end > evt_start && time_end > time_start)
		{
          time_start+=1; // Prefer to loose 1 sec than to use data with wrong setting
          time_end-=1; // Prefer to loose 1 sec than to use data with wrong setting 
          if(usetimestamp)
		  {
		    std::cout<<" Step : "<<voltage<<" "<<time_start<<"-"<<time_end<<"  start : "<<ctime(&time_start);

   	        std::pair< int, int > time_pair( time_start, time_end );
	        VoltageSteps_timestamp.insert(make_pair( voltage, time_pair));
          }
		  else
		  {
		    std::cout<<" Step : "<<voltage<<" "<<run<<" "<<evt_start<<"-"<<evt_end<<std::endl;

	        std::pair< int, int > evt_pair( evt_start, evt_end );
	        std::pair< int, std::pair< int,int > > run_evt_pair( run, evt_pair);
	        VoltageSteps_evtnumber.insert(make_pair( voltage, run_evt_pair));
		  }

		} // end Store step
	  }
	  
	  step_previous = step;
	  voltage_previous = voltage;
	  run_previous = run;
	  
	}

	fin.close();

	// start monitoring 5 min before scan start
	t_monitor_start-=60*5;
	t_monitor_start=t_monitor_start/60*60; // truncated to the lower minute
	// stop monitoring 5 min after scan stop
	t_monitor_end+=60*5;
	t_monitor_end=t_monitor_end/60*60 + 60; // rounded to the upper minute

	std::cout<<" Monitor : t_start "<<t_monitor_start<<"  t_end "<<t_monitor_end<<std::endl;
  }
  
  
 // Store voltage values removing double values
  std::multimap< int, std::pair< int,int > >::iterator itVStep_time;
  std::multimap< int, std::pair< int, std::pair< int,int > > >::iterator itVStep_evt;
  if(usetimestamp)
	for(itVStep_time=VoltageSteps_timestamp.begin(); itVStep_time!=VoltageSteps_timestamp.end(); itVStep_time++ )
      Voltage.insert(itVStep_time->first);
  else
	for(itVStep_evt=VoltageSteps_evtnumber.begin(); itVStep_evt!=VoltageSteps_evtnumber.end(); itVStep_evt++ )
      Voltage.insert(itVStep_evt->first);
	
}

void VoltageStepsTreeMaker::setVoltageSteps_evtnumber(const char* inputfile){
  
  std::string line;
  ifstream fin(inputfile);
  float voltage=0;
  int run=-1;
  int evt_start=-1;
  int evt_end=-1;
  
  if(fin.is_open())
  {
    while( getline ( fin, line) )
	{
	  if(fin.eof()) continue;
	  std::stringstream ss(line);
	  ss >> voltage >> run >> evt_start >> evt_end;
	  if(evt_end < evt_start) std::cout<<" >>> WARNING : evt_end < evt_start !"<<std::endl;
      std::cout<<" Step : "<<voltage<<" "<<run<<" "<<evt_start<<"-"<<evt_end<<std::endl;

	  std::pair< int, int > evt_pair( evt_start, evt_end );
	  std::pair< int, std::pair< int,int > > run_evt_pair( run, evt_pair);
	  VoltageSteps_evtnumber.insert(make_pair( voltage, run_evt_pair));
	}
	fin.close();
  }
  else
  {
    std::cerr<<" WARNING : file '"<<inputfile<<"' not found. Steps will not be set."<<std::endl;
  }
  
  // Store voltage values removing double values
  std::multimap< int, std::pair< int, std::pair< int,int > > >::iterator itVStep;
  for(itVStep=VoltageSteps_evtnumber.begin(); itVStep!=VoltageSteps_evtnumber.end(); itVStep++ )
    Voltage.insert(itVStep->first);
	
}

int VoltageStepsTreeMaker::getVoltageSteps_evtnumber(int run, int evt, int timestamp){

  std::multimap< int, std::pair< int, std::pair< int,int > > >::iterator itVStep;
  int theVoltage=-1;
  for(itVStep=VoltageSteps_evtnumber.begin(); itVStep!=VoltageSteps_evtnumber.end(); itVStep++ )
  { 
	std::pair< int, std::pair< int,int > > run_evt_pair = itVStep->second;
	std::pair< int,int > evt_pair = run_evt_pair.second;
	int stepRun = run_evt_pair.first;
	int evt_start = evt_pair.first;
	int evt_end = evt_pair.second;
	if ( run==stepRun && evt>=evt_start && evt<=evt_end ) theVoltage=itVStep->first;

    // When evtid selection used, compute corresponding timestamp ranges from selected events
	if(timestamp)
	if( run==stepRun && evt>=evt_start && evt<=evt_end )
	{
	  if(timestamp<VoltageSteps_timestart[itVStep->first]) VoltageSteps_timestart[itVStep->first]=timestamp;
	  if(timestamp>VoltageSteps_timeend[itVStep->first]) VoltageSteps_timeend[itVStep->first]=timestamp;
	}
  }
  return theVoltage;
}

void VoltageStepsTreeMaker::setVoltageSteps_timestamp(const char* inputfile){
  
  std::string line;
  ifstream fin(inputfile);
  float voltage=0;
  int time_start=-1;
  int time_end=-1;
  std::string str_start;
  std::string str_end;
  
  if(fin.is_open())
  {
    while( getline ( fin, line) )
	{
	  if(fin.eof()) continue;
	  std::stringstream ss(line);
	  ss >> voltage >> str_start >> str_end;
	  time_start = convertTimestamp( str_start );
	  time_end = convertTimestamp( str_end );
	  if(time_start<0 || time_end<0) std::cout<<" >>> WARNING : wrong timestamp conversion."<<std::endl;
	  
	  if(time_end < time_start) std::cout<<" >>> WARNING : time_end < time_start !"<<std::endl;
          time_start+=1; // Prefer to loose 1 sec than to use data with wrong setting
          time_end-=1; // Prefer to loose 1 sec than to use data with wrong setting 
      std::cout<<" Step : "<<voltage<<" "<<time_start<<"-"<<time_end<<std::endl;
      time_t tt = time_start;
	  std::cout<<"  start : "<<ctime(&tt);

	  std::pair< int, int > time_pair( time_start, time_end );
	  VoltageSteps_timestamp.insert(make_pair( voltage, time_pair));
	}
	fin.close();
  }
  else
  {
    std::cerr<<" WARNING : file '"<<inputfile<<"' not found. Steps will not be set."<<std::endl;
  }
  
  // Store voltage values removing double values
  std::multimap< int, std::pair< int,int > >::iterator itVStep;
  for(itVStep=VoltageSteps_timestamp.begin(); itVStep!=VoltageSteps_timestamp.end(); itVStep++ )
    Voltage.insert(itVStep->first);

}

int VoltageStepsTreeMaker::getVoltageSteps_timestamp(int timestamp){

  std::multimap< int, std::pair< int,int > >::iterator itVStep_t;
  int theVoltage=-1;
  for(itVStep_t=VoltageSteps_timestamp.begin(); itVStep_t!=VoltageSteps_timestamp.end(); itVStep_t++ )
  { 
	std::pair< int,int > time_pair = itVStep_t->second;
	int time_start = time_pair.first;
	int time_end = time_pair.second;
	if ( timestamp>=time_start && timestamp<=time_end ) theVoltage=itVStep_t->first;
  }
  return theVoltage;
}
