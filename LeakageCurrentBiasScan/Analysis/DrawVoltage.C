#include "Code/Voltage.h"
#include "TCanvas.h"

void DrawVoltage()
{

  TGraph* g = ReadVoltage();

  TCanvas c1;
  g->Draw("AP");
  c1.Modified();
  c1.Update();
  getchar();

}