#include "TChain.h"
#include "VoltageStepsTreeMaker.h"

int main()
{
  
  TChain c("ttree");
  c.Add("file:../../../clustersTree*.root/demo/ttree");

  int subdet=1; // 0 all subdet, 1 TIB, 2 TOB, 3 TID, 4 TEC
  bool usetimestamp=true; // true if input file uses timestamps
  int angle=0; // 0 all angles, 1 <20deg, ...
  int width=0; // 0 all width, 1 nstrips=1, ...
  int onstrip=0; // 0 not used, 1 on strip, 2 or -2 next to strip, 3 or -3 inter-strip
  VoltageStepsTreeMaker *t = new VoltageStepsTreeMaker(&c, "../Steps/Run_test.txt", subdet, usetimestamp, angle, width, onstrip );
  t->Loop();

}

