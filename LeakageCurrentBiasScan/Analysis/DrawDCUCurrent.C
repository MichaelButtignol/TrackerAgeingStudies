#include "Code/DCUCurrent.h"
#include "TCanvas.h"

void ConvertAllDCUToRoot()
{
  ConvertDCUCurrentTxtToRoot("Data/DCU_I_TIB_L1_20111012_run178367.txt");
  ConvertDCUCurrentTxtToRoot("Data/DCU_I_TIB_L1_20120405_run190459.txt");
  ConvertDCUCurrentTxtToRoot("Data/DCU_I_TIB_L1_20120506_run193541.txt");

  ConvertDCUCurrentTxtToRoot("Data/DCU_I_TIB_L2_20120405_run190459.txt");
  ConvertDCUCurrentTxtToRoot("Data/DCU_I_TIB_L2_20120506_run193541.txt");
}

void DrawDCUCurrent()
{
  //ConvertAllDCUToRoot();
  //TGraph* g = ReadDCUCurrentRoot("Data/DCU_I_TIB_L1_20120405_run190459.root", 369121605);
  TGraph* g = ReadCurrentFromGB();
  
  TCanvas c1;
  g->Draw("APL");
  c1.Modified();
  c1.Update();
  getchar();
}
