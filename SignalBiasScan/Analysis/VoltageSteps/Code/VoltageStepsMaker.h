#ifndef VoltageStepsMaker_h
#define VoltageStepsMaker_h

#include <map>
#include <vector>
#include <set>

class VoltageStepsMaker {

 public :

   VoltageStepsMaker(){ t_monitor_start=0; t_monitor_end=0; }

   // initialization of VoltageSteps_timestart, VoltageSteps_timeend
   void Initialize();

   // Reads and stores steps definition
   void readVoltageSteps_evtnumber(const char*);
   void readVoltageSteps_timestamp(const char*);
   // new steps definition format
   void readVoltageSteps(const char* inputfile, bool usetimestamp=false);
   
   // Returns voltage from evt number or timestamp
   int getVoltage_evtnumber(int run, int evt, int timestamp=0);
   int getVoltage_timestamp(int timestamp);

   // Returns voltage step index
   int getIndex(int voltage);
   
   // Print steps timestamp definition computed from evt read and definition from evt number  
   void printComputedSteps();
   
   // Returns voltage list
   std::set< int > getVoltageList() { return Voltage; }
   
   // Returns number of steps
   int getNVoltage() { return Voltage.size(); }

   unsigned int t_monitor_start;
   unsigned int t_monitor_end;
   
 private :
 
   // Compute index from voltage steps list
   void ComputeVoltageStepsIndex();
   
   // Voltage steps store in 'set' to merge those done in several parts
   std::set< int > Voltage;
   // Map to get voltage step index from voltage
   std::map< int, int > VoltageStepsIndex;
   
   std::multimap< int, std::pair< int, std::pair< int,int > > > VoltageSteps_evtnumber;

   std::multimap< int, std::pair< int,int > > VoltageSteps_timestamp;
   std::map< int, int > VoltageSteps_timestart;
   std::map< int, int > VoltageSteps_timeend;
 
};

#endif 
