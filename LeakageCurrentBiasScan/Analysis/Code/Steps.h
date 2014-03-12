#include <iostream>
#include <sstream>
#include <fstream>
#include <time.h>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1F.h"

int convertTimestamp( std::string str){

  std::stringstream ss;
  tm time;
  int year=-1;
  int month=-1;
  int hour=-1;

  ss.clear(); ss << str.substr(0, 4); ss >> year;
  time.tm_year = year - 1900;
  ss.clear(); ss << str.substr(4, 2); ss >> month;
  time.tm_mon = month - 1;
  ss.clear(); ss << str.substr(6, 2); ss >> time.tm_mday;
  ss.clear(); ss << str.substr(8, 2); ss >> hour; time.tm_hour=hour;
  ss.clear(); ss << str.substr(10, 2); ss >> time.tm_min;
  ss.clear(); ss << str.substr(12, 2); ss >> time.tm_sec;
  //time.tm_isdst=1;

  if(year < 2000 || year > 2020) std::cout<<" Wrong year format : "<<year<<std::endl;
  if(time.tm_sec < 0 || time.tm_sec > 61) std::cout<<" Wrong sec format : "<<year<<std::endl; 
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
  return out_time-3600; //->UTC
  
} 


TGraph* ReadSteps(char* filename, bool print=false)
{
  
  //Read file with steps definition
  
  std::string line;
  ifstream fin(filename);
  if(!fin.is_open()) { std::cout<<"Error : file "<<std::string(filename)<<" not found."<<std::endl; return 0;}
    
  int step=-1;
  int time=-1;
  std::string str_time;
  int event=-1;
  std::string orbit;
  float voltage=0;
  int run=-1;

  TGraph *g = new TGraph();
  g->SetName("VoltageSteps");

  int i=0;
  if(fin.is_open())  {
    while( getline ( fin, line) )
        {
          if(fin.eof()) continue;
          std::stringstream ss(line);
          ss >> step >> str_time >> event >> orbit >> voltage >> run;
          time = convertTimestamp( str_time );
          if(print) std::cout<<" Step : "<<voltage<<" "<<time<<std::endl;
          //time_t tt = time;
          //std::cout<<"  time : "<<ctime(&tt);
          g->SetPoint(i, time, voltage);
          i++;
        }
        fin.close();
  }

  TH1F* h = g->GetHistogram();
  h->GetXaxis()->SetTimeDisplay(1);
  h->GetXaxis()->SetTimeFormat("%H:%M");
  g->SetLineColor(2);

  return g;
}


void ReadBadPeriodsOld(char* filename, vector< int > &bad_periods_start, vector< int > &bad_periods_end)
{
  
  //Read file with bad periods definition
  
  std::string line;
  ifstream fin(filename);
  if(!fin.is_open()) { std::cout<<"Error : file "<<std::string(filename)<<" not found."<<std::endl; return;}
    
  int time_start=-1;
  std::string str_time_start;
  int time_end=-1;
  std::string str_time_end;

  N_periods=0;
  if(fin.is_open())  {
    while( getline ( fin, line) )
        {
          if(fin.eof()) continue;
          std::stringstream ss(line);
          //ss >> str_time_start >>str_time_end;
          //time_start = convertTimestamp( str_time_start );
          //time_end = convertTimestamp( str_time_end );
          ss >> time_start >> time_end;
		  bad_periods_start.push_back(time_start);
		  bad_periods_end.push_back(time_end);
          N_periods++;
        }
        fin.close();
  }

	//for(unsigned int ip=0; ip<N_periods; ip++)
	  //cout<<"bad period "<<ip<<" :  start "<<bad_periods_start[ip]<<"  end "<<bad_periods_end[ip]<<endl;

  return;
}
