#ifndef DCUCurrent_h
#define DCUCurrent_h

#include <iostream>
#include <sstream>
#include <fstream>
#include <time.h>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1F.h"
// for timestamp conversion
#include "./Steps.h"

Long64_t convertDate( std::string str_date, std::string str_time){

  std::stringstream ss;
  tm time;
  int year=-1;
  std::string str_month;
  int month=-1;
  //int day=-1;
  int hour=-1;

  ss.clear(); ss << str_date.substr(0, 2); ss >> time.tm_mday;

  int first_ = str_date.find("-");
  int second_ = str_date.find("-", first_+1);
  //cout<<first_<<"-"<<second_<<endl;
  str_month = str_date.substr(first_+1, second_-first_-1);
  //cout << str_month <<endl;
  if(str_month=="JAN") month = 1;
  if(str_month=="FEB") month = 2;
  if(str_month=="MAR") month = 3;
  if(str_month=="APR") month = 4;
  if(str_month=="MAY") month = 5;
  if(str_month=="JUN") month = 6;
  if(str_month=="JUL") month = 7;
  if(str_month=="AUG") month = 8;
  if(str_month=="SEP") month = 9;
  if(str_month=="OCT") month = 10;
  if(str_month=="NOV") month = 11;
  if(str_month=="DEC") month = 12;
  time.tm_mon = month - 1;
  if(month<0) { std::cout<<"Error : month '"<<str_month<<"' not known."<<std::endl; return -1;}

  ss.clear(); ss << str_date.substr(second_+1, 2); ss >> year;
  time.tm_year = year + 2000 - 1900;

  ss.clear(); ss << str_time.substr(0, 2); ss >> hour; time.tm_hour=hour;
  ss.clear(); ss << str_time.substr(3, 2); ss >> time.tm_min;
  ss.clear(); ss << str_time.substr(6, 2); ss >> time.tm_sec;
  //time.tm_isdst=1;

  if(year+2000 < 2000 || year+2000 > 2020) std::cout<<" Wrong year format : "<<year+2000<<std::endl;
  if(time.tm_sec < 0 || time.tm_sec > 61) std::cout<<" Wrong sec format : "<<time.tm_sec<<std::endl; 
  //cout<<" timestamp "<<year<<" "<<time.tm_mon<<" "<<time.tm_mday<<" "<<time.tm_hour<<" "<<time.tm_min<<" "<<time.tm_sec<<std::endl;
 
 
  time_t out_time = mktime( &time );
  // Next lines to avoid random shifts of 1 hour that happen sometimes when doing conversion !!
  struct tm * timeinfo;
  timeinfo = localtime(&out_time); 
  if(hour!=timeinfo->tm_hour) 
  {
    //cout<<"try to correct hour"<<endl;
    time.tm_hour=hour;
    out_time = mktime( &time );
    timeinfo = localtime(&out_time);
  }
  if(hour!=timeinfo->tm_hour) {std::cout<<"ERROR in hour : "<<hour<<" "<<timeinfo->tm_hour<<std::endl; return 0;}
  return (Long64_t) out_time-3600; //->UTC

} 

void ReadBadPeriods(char* filename, vector< int > &bad_periods_start, vector< int > &bad_periods_end)
{
  
  //Read file with bad periods definition
  
  std::string line;
  ifstream fin(filename);
  if(!fin.is_open()) { std::cout<<"Error : file "<<std::string(filename)<<" not found."<<std::endl; return;}
    
  int time_start=-1;
  std::string str_time_start;
  int time_end=-1;
  std::string str_time_end;

  if(fin.is_open())  {
    while( getline ( fin, line) )
        {
          if(fin.eof()) continue;
          std::stringstream ss(line);
          ss >> str_time_start >>str_time_end;
          //time_start = convertTimestamp( str_time_start );
          //time_end = convertTimestamp( str_time_end );
          //ss >> time_start >> time_end;
		  bad_periods_start.push_back(convertTimestamp(str_time_start));
		  bad_periods_end.push_back(convertTimestamp(str_time_end));
        }
        fin.close();
  }

	//for(unsigned int ip=0; ip<bad_periods_start.size(); ip++)
	  //cout<<"bad period "<<ip<<" :  start "<<bad_periods_start[ip]<<"  end "<<bad_periods_end[ip]<<endl;

  return;
}


TGraph* ReadDCUCurrentTxt(char* filename="Data/DCU_I_TIB_L1_20120405_run190459.txt", int detid=369121605, char* bad_periods="")
{
  
  // Read bad periods
  vector< int > bad_periods_start;
  vector< int > bad_periods_end;
  if(strcmp(bad_periods, ""))
  {
    ReadBadPeriods(bad_periods, bad_periods_start, bad_periods_end);
	if(bad_periods_start.size() != bad_periods_end.size())
	{
	  std::cerr<<"Wrong definition of bad periods : size of starts != size of ends. Will not use them."<<endl;
	  bad_periods="";
	}
  }

  std::string line;
  ifstream fin(filename);
  std::string ps;
  int modid=-1;
  std::string str_date;
  std::string str_time;
  Int_t time=-1;
  float current=-1;

  TGraph *g = new TGraph();
  g->SetName("Current");

  int i=0;
  if(fin.is_open())  {
    while( getline ( fin, line) )
        {
          if(fin.eof()) continue;
          std::stringstream ss(line);
          ss >> ps >> modid >> str_date >> str_time >> current;
          time = convertDate( str_date , str_time);

          // Remove points during bad periods
          if(strcmp(bad_periods, ""))
		  {
		     for(unsigned int ip=0; ip<bad_periods_start.size(); ip++)
			   if(time>=bad_periods_start[ip] && time<=bad_periods_end[ip]) continue;
		  }

          if(modid==detid) 
          {
            std::cout<<str_date<<" "<<str_time<<" "<<current<<std::endl;
            std::cout<<time<<std::endl;
            g->SetPoint(i, time, current);
          }
          i++;
        }
        fin.close();
  }

  TH1F* h = g->GetHistogram();
  h->GetXaxis()->SetTimeDisplay(1);
  h->GetXaxis()->SetTimeFormat("%H:%M");
  
  return g;

}

void ConvertDCUCurrentTxtToRoot(char* filename="Data/DCU_I_TIB_L1_20120405_run190459.txt")
{
  
  cout<<"Converting file "<<filename<<" to root format."<<endl;
  
  // Create output file
  TString file(filename);
  int dot=file.Index(".txt");
  file.Replace(dot,4,".root");
  TFile fout(file.Data(), "recreate");

  // Create output tree
  char ps[50];
  int detid=-1;
  Long64_t time=0;
  float current=-1;
  TTree* tree = new TTree("dcu", "dcu");
  tree->Branch("PS", &ps, "PS/C");
  tree->Branch("DETID", &detid, "DETID/I");
  tree->Branch("TIME", &time, "TIME/L");
  tree->Branch("CURRENT", &current, "CURRENT/F");

  // Read input file
  std::string line;
  ifstream fin(filename);
  std::string str_ps;
  std::string str_date;
  std::string str_time;
  int prev_detid=0;

  if(fin.is_open())  {
    while( getline ( fin, line) )
        {
          if(fin.eof()) continue;
          std::stringstream ss(line);
          ss >> str_ps >> detid >> str_date >> str_time >> current;
          sprintf(ps, "%s", str_ps.c_str());
          time = convertDate( str_date , str_time);
           //  std::cout<<str_date<<" "<<str_time<<" "<<current<<std::endl;
          //  std::cout<<time<<endl;
          if(detid!=prev_detid) std::cout<<"detid "<<detid<<std::endl;
          tree->Fill();
          prev_detid=detid;
        }
        fin.close();
  }

  // Save tree
  tree->Write();
  fout.Close();
  
}

TGraph* ReadCurrentRoot(char* filename="Data/DCU_I_TIB_L1_20120405_run190459.root", int modid=369121605, char* treename="dcu",
 char* bad_periods="", bool print=false)
{

  // Read bad periods
  vector< int > bad_periods_start;
  vector< int > bad_periods_end;
  if(strcmp(bad_periods, ""))
  {
    ReadBadPeriods(bad_periods, bad_periods_start, bad_periods_end);
	if(bad_periods_start.size() != bad_periods_end.size())
	{
	  std::cerr<<"Wrong definition of bad periods : size of starts != size of ends. Will not use them."<<endl;
	  bad_periods="";
	}
	for(unsigned int ip=0; ip<bad_periods_start.size(); ip++)
	  cout<<"bad period "<<ip<<" :  start "<<bad_periods_start[ip]<<"  end "<<bad_periods_end[ip]<<endl;
  }

  // Open file and set tree
  TFile* fin = TFile::Open(filename);
  if(!fin) {std::cout<<"Error : file '"<<std::string(filename)<<"' not found."<<std::endl; return 0;}
  
  char ps[50];
  int detid=-1;
  Long64_t time=0;
  float current=-1;

  TTree* tree = (TTree*) fin->Get(treename);
  if(!tree) {std::cout<<"Error : tree not found."<<std::endl; return 0;}
  tree->SetBranchAddress("PS", &ps);
  tree->SetBranchAddress("DETID", &detid);
  tree->SetBranchAddress("TIME", &time);
  tree->SetBranchAddress("CURRENT", &current);

  // Create graph
  TGraph *g = new TGraph();
  g->SetName("Current");

  // Read tree and fill graph
  int nentries = tree->GetEntries();
  int ipt=0;
  bool remove_point=false;
  for(int i=0; i<nentries; i++)
  {
    tree->GetEntry(i);

    // Remove points during bad periods
	if(strcmp(bad_periods, ""))
	{
       remove_point=false;
	   for(unsigned int ip=0; ip<bad_periods_start.size(); ip++)
	   {
		 if(time>=bad_periods_start[ip] && time<=bad_periods_end[ip]) remove_point=true;
		 //cout<<time<<"  start "<<bad_periods_start[ip]<<"  end "<<bad_periods_end[ip]<<" "<<remove_point<<endl;
	   }
	}
	if(remove_point) continue;

    if(detid==modid) 
	{ 
	  if(print) cout<<" I_dcu "<<current<<" t "<<time<<endl;
	  g->SetPoint(ipt, time, current);
	  ipt++;
	}
  }

  // Set time axis and markers
  TH1F* h = g->GetHistogram();
  h->GetXaxis()->SetTimeDisplay(1);
  h->GetXaxis()->SetTimeFormat("%H:%M");

  fin->Close();
  
  return g;

}

TGraph* ReadDCUCurrentRoot(char* filename="Data/DCU_I_TIB_L1_20120405_run190459.root", int modid=369121605, char* bad_periods="")
{

  TGraph* g = ReadCurrentRoot(filename, modid, "dcu", bad_periods); // Same tree format for DCU and PS currents
  if(!g) return 0;
  
  g->SetMarkerStyle(22);
  g->SetMarkerColor(4);
  
  return g;
  
}

TGraph* ReadDCUCurrentFromGB(char* filename="~/work/DCU_TIBD_TOB_from_1348837200_to_1348862400.root", int modid=369121605,
 char* bad_periods="", bool print=false, TGraph* gtemp=0)
{

  // Read bad periods
  vector< int > bad_periods_start;
  vector< int > bad_periods_end;
  if(strcmp(bad_periods, ""))
  {
    ReadBadPeriods(bad_periods, bad_periods_start, bad_periods_end);
	if(bad_periods_start.size() != bad_periods_end.size())
	{
	  std::cerr<<"Wrong definition of bad periods : size of starts != size of ends. Will not use them."<<endl;
	  bad_periods="";
	}
	for(unsigned int ip=0; ip<bad_periods_start.size(); ip++)
	  cout<<"bad period "<<ip<<" :  start "<<bad_periods_start[ip]<<"  end "<<bad_periods_end[ip]<<endl;
  }

  // Open file and set tree
  TFile* fin = TFile::Open(filename);
  if(!fin) {std::cout<<"Error : file '"<<std::string(filename)<<"' not found."<<std::endl; return 0;}
  
  double Detid;
  double Dcutimestamp;
  double TemperatureSi;
  double Ileak;

  TTree* tree = (TTree*) fin->Get("outTree");
  if(!tree) {std::cout<<"Error : tree not found."<<std::endl; return 0;}
  tree->SetBranchAddress("Detid", &Detid);
  tree->SetBranchAddress("Dcutimestamp", &Dcutimestamp);
  tree->SetBranchAddress("TemperatureSi", &TemperatureSi);
  tree->SetBranchAddress("Ileak", &Ileak);

  // Create graph
  TGraph *g = new TGraph();
  g->SetName("Current");

  // Read tree and fill graph
  int nentries = tree->GetEntries();
  int ipt=0;
  bool remove_point=false;
  double time;
  for(int i=0; i<nentries; i++)
  {
    tree->GetEntry(i);
	time = Dcutimestamp-3600*3;

    // Remove points during bad periods
	if(strcmp(bad_periods, ""))
	{
       remove_point=false;
	   for(unsigned int ip=0; ip<bad_periods_start.size(); ip++)
	   {
		 if(time>=bad_periods_start[ip] && time<=bad_periods_end[ip]) remove_point=true;
		 //cout<<time<<"  start "<<bad_periods_start[ip]<<"  end "<<bad_periods_end[ip]<<" "<<remove_point<<endl;
	   }
	}
	if(remove_point) continue;

    if(Detid==modid) 
	{ 
	  if(print) cout<<" I_dcu "<<Ileak<<" t "<<time<<endl;
	  if(Ileak<20) continue;
	  g->SetPoint(ipt, time, Ileak);
	  if(gtemp) gtemp->SetPoint(ipt, time, TemperatureSi);
	  ipt++;
	}
  }

  // Set time axis and markers
  TH1F* h = g->GetHistogram();
  h->GetXaxis()->SetTimeDisplay(1);
  h->GetXaxis()->SetTimeFormat("%H:%M");
  g->SetMarkerStyle(22);
  g->SetMarkerColor(4);

  fin->Close();
  
  return g;

}

#endif
