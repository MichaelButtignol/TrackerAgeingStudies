#include "Code/DCUCurrent.h"

TGraph* ReadPSCurrentTxt(char* filename="Data/PS_I_TIB_L1_20120405_run190459.txt", int detid=369121605, char* bad_periods="")
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
  std::string PS;
  std::string clob;
  std::string modids;
  int modid1=-1;
  int modid2=-1;
  int modid3=-1;
  std::string str_date;
  std::string str_time;
  Int_t time=-1;
  float current=-1;

  TGraph *g = new TGraph();
  g->SetName("PS_Current");

  int i=0;
  if(fin.is_open())  {
    while( getline ( fin, line) )
        {
          if(fin.eof()) continue;
          current=-1;
          std::stringstream ss(line);
          ss >> PS >> clob >> modids >> str_date >> str_time >> current;
          // some format without (CLOB)
          if(current==-1)
          { 
            std::stringstream ss(line);
            ss >> PS >> modids >> str_date >> str_time >> current;
          }
          //cout<<PS<<"; "<<modids<<"; "<<str_date<<"; "<<str_time<<"; "<<current<<endl;
          time = convertDate( str_date , str_time);

          // Remove points during bad periods
          if(strcmp(bad_periods, ""))
		  {
		     for(unsigned int ip=0; ip<bad_periods_start.size(); ip++)
			   if(time>=bad_periods_start[ip] && time<=bad_periods_end[ip]) continue;
		  }

          //cout<<modids<<endl;
          std::stringstream ss2;
          ss2.clear(); ss2 << modids.substr(0, 9); ss2 >> modid1;
          ss2.clear(); ss2 << modids.substr(10, 9); ss2 >> modid2;
          ss2.clear(); ss2 << modids.substr(20, 9); ss2 >> modid3;
          //cout<<modid1<<" "<<modid2<<" "<<modid3<<endl;

          if(detid==modid1 || detid==modid2 || detid==modid3) 
          {
            std::cout<<str_date<<" "<<str_time<<" "<<current<<std::endl;
            std::cout<<time<<std::endl;
            g->SetPoint(i, time, current);
            i++;
          }
        }
    fin.close();
  }

  TH1F* h = g->GetHistogram();
  h->GetXaxis()->SetTimeDisplay(1);
  h->GetXaxis()->SetTimeFormat("%H:%M");
  
  return g;

}

void ConvertPSCurrentTxtToRoot(char* filename="Data/PS_I_TIB_L1_20120405_run190459.txt")
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
  TTree* tree = new TTree("ps", "ps");
  tree->Branch("PS", &ps, "PS/C");
  tree->Branch("DETID", &detid, "DETID/I");
  tree->Branch("TIME", &time, "TIME/L");
  tree->Branch("CURRENT", &current, "CURRENT/F");

  // Read input file
  std::string line;
  ifstream fin(filename);
  std::string str_ps;
  std::string str_clob;
  std::string str_modids;
  int modid1=-1;
  int modid2=-1;
  int modid3=-1;
  std::string str_date;
  std::string str_time;
  std::string prev_str_ps;

  if(fin.is_open())  {
    while( getline ( fin, line) )
        {
          if(fin.eof()) continue;
          std::stringstream ss(line);
          current=-1;
          ss >> str_ps >> str_clob >> str_modids >> str_date >> str_time >> current;
          if(current==-1)
          {
            std::stringstream ss(line);
            ss >> str_ps >> str_modids >> str_date >> str_time >> current;
          }
          sprintf(ps, "%s", str_ps.c_str());
          time = convertDate( str_date , str_time);

          //cout<<modids<<endl;
          std::stringstream ss2;
          ss2.clear(); ss2 << str_modids.substr(0, 9); ss2 >> modid1;
          ss2.clear(); ss2 << str_modids.substr(10, 9); ss2 >> modid2;
          ss2.clear(); ss2 << str_modids.substr(20, 9); ss2 >> modid3;
          //cout<<modid1<<" "<<modid2<<" "<<modid3<<endl;

          if(str_ps!=prev_str_ps) std::cout<<str_ps<<"   "<<str_modids<<std::endl;
          
          detid=modid1;
          tree->Fill();
          detid=modid2;
          tree->Fill();
          detid=modid3;
          tree->Fill();
          
          prev_str_ps=str_ps;
        }
    fin.close();
  }

  // Save tree
  tree->Write();
  fout.Close();
  
}

TGraph* ReadPSCurrentRoot(char* filename="Data/PS_I_TIB_L1_20120405_run190459.root", int modid=369121605, char* bad_periods="",
bool print=false)
{
  
  TGraph* g = ReadCurrentRoot(filename, modid, "ps", bad_periods, print); // Same tree format for DCU and PS currents
  
  g->SetMarkerStyle(23);
  g->SetMarkerColor(8);
  g->SetLineColor(8);
  
  return g;
  
}
