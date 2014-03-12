#include <fstream>
#include <sstream>
#include "CommonTools/TrackerMap/interface/TrackerMap.h"

int main(int argc, char **argv)
{

  std::cout<<"Syntax : TrackerMapMaker [input vmin vmax output]"<<std::endl;

  // initialization
  std::string inputfile="tkMap.txt";
  if(argc>1) inputfile=argv[1];

  float vmin=0;
  float vmax=0; 
  if(argc>2) vmin=atof(argv[2]);
  if(argc>3) vmax=atof(argv[3]);

  std::string outputfile="tkMap.pdf";
  if(argc>4) outputfile=argv[4];

  std::string line;
  ifstream fin(inputfile);
  int detid;
  float value;

  TrackerMap tkMap;

  if(fin.is_open())
  {
    while( getline ( fin, line) )
	{
	  if(fin.eof()) continue;
	  std::stringstream ss(line);
	  ss >> detid >> value;
          tkMap.fill(detid, value);
	}
	fin.close();
    tkMap.save(false, vmin, vmax, outputfile.c_str() );
  }
  else
  {
    std::cerr<<" WARNING : file '"<<inputfile<<"' not found. Exit."<<std::endl;
  }
 


}
