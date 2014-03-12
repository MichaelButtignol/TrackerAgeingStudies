#include "Code/Steps.h"
#include "TCanvas.h"

void DrawSteps()
{

  TGraph* g = ReadSteps("Steps/Steps_20120928_run203832.txt", true);

  TCanvas c1;
  g->Draw("APL");
  c1.Modified();
  c1.Update();
  getchar();
}
