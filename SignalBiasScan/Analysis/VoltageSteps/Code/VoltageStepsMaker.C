#define VoltageStepsMaker_cxx
#include "VoltageStepsMaker.h"

#include <iostream>
#include <time.h>
#include <sstream>
#include <fstream>

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


void VoltageStepsMaker::Initialize(){

  std::set< int >::iterator itVolt;
  for( itVolt=Voltage.begin(); itVolt!=Voltage.end(); itVolt++){
    VoltageSteps_timestart[*itVolt] = 2147483647;//int max
    VoltageSteps_timeend[*itVolt] = 0;
  }
}


void VoltageStepsMaker::ComputeVoltageStepsIndex()
{
  std::set< int >::iterator itVolt;
  int istep=0;
  for( itVolt=Voltage.begin(); itVolt!=Voltage.end(); itVolt++){
    VoltageStepsIndex.insert( make_pair(*itVolt, istep) );
	istep++;
  }   
}


// Methods to read and store steps definition
//-------------------------------------------- 

void VoltageStepsMaker::readVoltageSteps(const char* inputfile, bool usetimestamp){
  
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
	    //std::cout<<step<<" t_start "<<time_start<<" "<<t_monitor_start<<std::endl;
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
		if(run_previous==run && evt_end > evt_start && time_end >= time_start)
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

    
	// Compute timestamps just before and just after the whole scan

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

  ComputeVoltageStepsIndex();
}


// Reads old steps definition file with evt numbers
void VoltageStepsMaker::readVoltageSteps_evtnumber(const char* inputfile){
  
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
	
  ComputeVoltageStepsIndex();
}


// Reads old steps definition file with timestamps
void VoltageStepsMaker::readVoltageSteps_timestamp(const char* inputfile){
  
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

  ComputeVoltageStepsIndex();
}



// Methods to know to which step an evt belong
//--------------------------------------------- 

// Returns voltage from run and evt number
int VoltageStepsMaker::getVoltage_evtnumber(int run, int evt, int timestamp){

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
	// After running on all evts, can get timestamps definition of each steps
	if(timestamp)
	if( run==stepRun && evt>=evt_start && evt<=evt_end )
	{
	  if(timestamp<VoltageSteps_timestart[itVStep->first]) VoltageSteps_timestart[itVStep->first]=timestamp;
	  if(timestamp>VoltageSteps_timeend[itVStep->first]) VoltageSteps_timeend[itVStep->first]=timestamp;
	}
  }
  return theVoltage;
}


// Returns voltage from timestamp
int VoltageStepsMaker::getVoltage_timestamp(int timestamp){

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


// Returns voltage step index
int VoltageStepsMaker::getIndex(int voltage){

  int index = -1;
  index =  VoltageStepsIndex[voltage];
  
  // throws an out-of-range exception if voltage not in VoltageStepsIndex;
  
  return index;

}
   

// Print steps timestamp definition computed from evt read and definition from evt number  
void VoltageStepsMaker::printComputedSteps(){

  std::cout << "Time definition of steps : " << std::endl;

  std::set< int >::iterator itVolt;
  for( itVolt=Voltage.begin(); itVolt!=Voltage.end(); itVolt++){
    std::cout<<(*itVolt)<<" "<<VoltageSteps_timestart[*itVolt]<<" "<<VoltageSteps_timeend[*itVolt]<<std::endl;
    time_t tt = VoltageSteps_timestart[*itVolt];
	std::cout<<"   start : "<<ctime(&tt);
    tt = VoltageSteps_timeend[*itVolt];
	std::cout<<"     end : "<<ctime(&tt);
  }

}
