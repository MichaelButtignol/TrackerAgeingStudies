//#include "../CommonTools/CurvesFunctions.C"
#include "../../../LeakageCurrentBiasScan/Analysis/GetKink.C"
#include "TLine.h"


Double_t fitsigmo(Double_t *x, Double_t *par){
  Double_t value, value_th, value_0;
  if (x[0]<par[2]) {
    value = par[0]/(1+exp(-1*par[1]*(x[0]+par[3])));//+ par[4]*x[0]; // sigmoid
  }
  else
  {
    value_th = par[0]/(1+exp(-1*par[1]*(par[2]+par[3])));
    value_0 = value_th - par[4]*par[2];
    value = value_0 + par[4]*x[0]; // pol1
  }
  //else { value=0; TF1::RejectPoint();}
  
  return value;
}

Double_t fitpol(Double_t *x, Double_t *par){
  Double_t value, value_th, value_0;
  if (x[0]<par[4]) {
	value = par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
  }
  else
  {
    value_th = par[0] + par[1]*par[4] + par[2]*par[4]*par[4] + par[3]*par[4]*par[4]*par[4];
    value_0 = value_th - par[5]*par[4];
    value = value_0 + par[5]*x[0]; // pol1
  }
  //else { value=0; TF1::RejectPoint();}
  
  return value;
}

double FitCurve(TGraphErrors* g, int debug=0)
{
  
  double vdep=0;
  if(!g) return vdep;
  
  
  TF1* f3 = new TF1("fp1", "pol1", 20, 360);
  f3->SetLineColor(4);
  
  double y;
  double xlow=350;
  double ymax=0;
  double vdep1=0;
  // prendre le point le plus haut est robuste et marche bien pour les hauts vdep
  TGraphErrors *gsmooth = MedianFilter(g);
  gsmooth = MedianFilter(gsmooth);
  for(int ipt=0; ipt<gsmooth->GetN(); ipt++)
  {
    gsmooth->GetPoint(ipt, xlow, y);
	if(y>ymax) {ymax=y; vdep1=xlow;}
  }
  
  //g=gsmooth;
  int npt = g->GetN()-5;
  xlow=350;
  double chi2=0;
  int status = 0;
  while(chi2<5. && xlow>100 && npt>=0)
  {
    g->GetPoint(npt, xlow, y);
    f3->SetRange(xlow, 350);
    status = g->Fit("fp1", "rqn");
    chi2 = f3->GetChisquare()/f3->GetNDF();
    if(debug>=1) cout<<"xlow "<<xlow<<" chi2 "<<chi2<<endl;
	npt--;
  }
  g->GetPoint(npt+2, xlow, y);
  f3->SetRange(xlow, 350);
  g->Fit("fp1", "rqn");
  vdep = xlow;
  
  //if(vdep>230) vdep = vdep1;
  //else if(fabs(vdep-vdep1)>40) vdep = vdep1;
  //vdep=vdep1;

  cout<<" Vdepl = "<<vdep<<endl;
  
  
  
  TF1* f1 = new TF1("fitsigmo", fitsigmo, 50, 350, 5);
  f1->SetParameter(0, 3.);
  f1->SetParameter(1, 0.02);
  f1->SetParameter(2, vdep);
  f1->SetParameter(3, -10.);
  f1->SetParameter(4, -0.00001);
  f1->SetParLimits(2, 100., 300.);
  f1->SetParLimits(4, -0.1, 0);
  
  
  
  
  TF1* f2 = new TF1("fitpol", fitpol, 20, 360, 6);
  f2->SetParLimits(4, 100., 300.);
  f2->SetParLimits(5, -0.1, 0);
  
 /* xlow = 30;
  chi2 = 100;
  double xmin = xlow;
  double chi2min = chi2;
  string opt="rq";
  if(debug>=2) opt="r";
  while(chi2>0.1 && xlow<vdep-80)
  {
    f2->SetParameter(0, 1.);
    f2->SetParameter(1, 0.005);
    f2->SetParameter(2, 0.0001);
    f2->SetParameter(3, 0.);
    f2->SetParameter(4, vdep);
    f2->SetParameter(5, -0.00001);
    f2->SetRange(xlow, 350);
    status = g->Fit("fitpol", opt.c_str());
    if(status!=0) status = g->Fit("fitpol", opt.c_str());
    chi2 = f2->GetChisquare()/f2->GetNDF();
    if(debug>=1) cout<<"xlow "<<xlow<<" chi2 "<<chi2<<endl;
    if(debug>=2) cout<<" chi2 "<<f2->GetChisquare()<<" ndf "<<f2->GetNDF()<<endl;
	if(chi2<chi2min) {chi2min=chi2; xmin=xlow;}
	xlow+=10;
  }
  f2->SetParameter(0, 1.);
  f2->SetParameter(1, 0.005);
  f2->SetParameter(2, 0.0001);
  f2->SetParameter(3, 0.);
  f2->SetParameter(4, vdep);
  f2->SetParameter(5, -0.00001);
  f2->SetRange(xmin, 350);
  status = g->Fit("fitpol", opt.c_str());
  if(status!=0) status = g->Fit("fitpol", opt.c_str());
  vdep = f2->GetParameter(4);
  cout<<" Vdepl2 = "<<vdep<<" xmin "<<xmin<<endl;
  if(fabs(vdep-vdep1)>30) vdep = vdep1;
*/




      TGraphErrors* gmedian;
      gmedian = MedianFilter( g );
      int nfilt=1;
      /*while (!IsMonoton(gmedian) && nfilt<4) {
        gmedian = MedianFilter( gmedian );
        nfilt++;
      }*/
      gmedian = HanningFilter(gsmooth);
      cout<<nfilt<<" median filter applied"<<endl;

      TGraphErrors* gscurv = GetCurvatureGraph( gmedian );
      gscurv->SetMarkerStyle(20);
      TGraph* g3pts = new TGraph();
      float xopt = GetOptimalMinNPts(gscurv, g3pts);

      vdep = xopt;

      cout<<" Kink = "<<xopt<<endl;


  TCanvas *c1;
  TCanvas *c2;
  if(debug>=1)
  {
    c1 = new TCanvas();
    g->Draw("AP");
    f2->Draw("same");
    f3->Draw("same");
	//gsmooth->SetMarkerColor(17);
	//gsmooth->Draw("P");
//gmedian->SetMarkerColor(4);
//gmedian->Draw("P");	
	TLine *l = new TLine(vdep,2, vdep, ymax+0.1);
	l->SetLineStyle(3);
	l->Draw();
 
    c1->Modified();
    c1->Update();

    string detid = g->GetTitle();
    detid.erase(0,6);
    //c1->Print(Form("ClusterWidthVsVbias_detid_%s_run203832.eps", detid.c_str()));

    c2 = new TCanvas;
    gscurv->Draw("AP");
    g3pts->Draw("P");
    g3pts->Fit("pol2", "q");
    c2->Modified();
    c2->Update();

    //c2->Print(Form("ClusterWidthVsVbias_detid_%s_curv_run203832.eps", detid.c_str()));

    getchar();
    c1->Close();
    c2->Close();
  }

  return vdep;
}

void FitOneCurve(string subdet, string run, ULong64_t modid, int debug=0)
{
  
  TGraphErrors* g = GetClusterWidthGraph(subdet, run, modid, 0);
  if(!g) return;

  // Correct for voltage drop due to leakage current
  string corr_name="_"+subdet+run;
  int corrected = CorrectGraphForLeakageCurrent(g, modid, corr_name);
  if(corrected) g->SetMarkerColor(13);

  g->SetTitle(Form("DetID %i", (int) modid));
  g->SetName("g");

  FitCurve(g, debug);
}


void FitTIBCurves(string run)
{
  // L1
  FitOneCurve("TIB", run, 369121606, 1);
  FitOneCurve("TIB", run, 369124429, 1);
  FitOneCurve("TIB", run, 369124550, 1);
  FitOneCurve("TIB", run, 369121365, 1);
  FitOneCurve("TIB", run, 369121302, 1);
  FitOneCurve("TIB", run, 369170776, 1);

  // L4
  FitOneCurve("TIB", run, 369169980, 1);
  FitOneCurve("TIB", run, 369169864, 1);
  FitOneCurve("TIB", run, 369170636, 1);
  FitOneCurve("TIB", run, 369173944, 1);

}

void FitTIBSmallScan(string run)
{
 //TIBminus_1_2_2_1
  FitOneCurve("TIB", run, 369121381, 1);
  FitOneCurve("TIB", run, 369121382, 1);
  FitOneCurve("TIB", run, 369121385, 1);
  FitOneCurve("TIB", run, 369121386, 1);
  FitOneCurve("TIB", run, 369121389, 1);
  FitOneCurve("TIB", run, 369121390, 1);
 //TIBminus_1_4_2_5
  FitOneCurve("TIB", run, 369121605, 1);
  FitOneCurve("TIB", run, 369121606, 1);
  FitOneCurve("TIB", run, 369121609, 1);
  FitOneCurve("TIB", run, 369121610, 1);
  FitOneCurve("TIB", run, 369121613, 1);
  FitOneCurve("TIB", run, 369121614, 1);
 //TIBplus_1_6_2_5
  FitOneCurve("TIB", run, 369125861, 1);
  FitOneCurve("TIB", run, 369125862, 1);
  FitOneCurve("TIB", run, 369125865, 1);
  FitOneCurve("TIB", run, 369125866, 1);
  FitOneCurve("TIB", run, 369125869, 1);
  FitOneCurve("TIB", run, 369125870, 1);
}

void Fit()
{
  //FitOneCurve("TIB", "_208339", 369125869, 1);
  //FitOneCurve("TIB", "_200786", 369136953, 1);
  //FitTIBCurves("_190459");
  FitTIBSmallScan("_203832");
}
