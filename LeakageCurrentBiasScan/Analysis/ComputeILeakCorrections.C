#include "Code/Steps.h"
#include "Code/Voltage.h"
#include "Code/DCUCurrent.h"
#include "Code/PSCurrent.h"
#include "Code/FitFunctions.h"

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"

double GetMaximum(TGraph* g)
{
  double max=-999;
  double x, y;
  for(unsigned int ipt=0; ipt<g->GetN(); ipt++)
  {
    g->GetPoint(ipt, x, y);
	if(y>max) max=y;
  }
  
  return max;
}

void Scale(TGraph &* g, double scale)
{
  double x, y;
  for(unsigned int ipt=0; ipt<g->GetN(); ipt++)
  {
    g->GetPoint(ipt, x, y);
	g->SetPoint(ipt, x, y*scale);
  }

  return;  
}

void GetConditions(TGraph *&gsteps, TGraph *&gcur_DCU, TGraph *&gcur_PS, TGraph *&gvolt, 
                   char* subdet="TIB_L1", char* run="20120506_run193541", int detid=369121605, char* bad_periods="")
{

  // Read files with voltage infos
  gsteps = ReadSteps(Form("Steps/Steps_%s.txt", run));
  if(!gsteps) {std::cout<<" No voltage steps info. Exit."<<std::endl; return;}
  gvolt = 0;//ReadVoltage(Form("Data/ConditionBrowser_%s.root", run));
  //if(!gvolt) std::cout<<" ConditionBrowser file does not exist, but it is not a problem."<<std::endl;
  
  // Read files with current infos
  //gcur_DCU = ReadDCUCurrentRoot(Form("Data/DCU_I_%s_%s.root", subdet, run), detid, bad_periods);
  gcur_DCU = ReadDCUCurrentFromGB("~/work/DCU_TIBD_TOB_from_1348837200_to_1348862400.root", detid, bad_periods);
  gcur_PS = ReadPSCurrentRoot(Form("Data/PS_I_%s_%s.root", subdet, run), detid, bad_periods);
  
  if(!gcur_PS) {std::cout<<" No PS current info. Exit."<<std::endl; return;}
  //gcur_DCU = 0; // force to not use DCU
  
  // If no DCU info, will user only PS info
  float scale=1.;
  double time=0;
  double current=0;
  if(!gcur_DCU)
  { 
    std::cout<<" No DCU current info, will try to use PS current instead."<<std::endl;
    gcur_DCU = (TGraph*) gcur_PS->Clone();
    if(!gcur_DCU) {std::cout<<" No PS info neither. Exit."<<std::endl; return;}
    
    scale = 1./3; // By default consider 3 modules per PS channel (not always true !)

    // fractions computed on 20120506_run193541
    const int N_2012=17;
    int detids_2012[N_2012]={369121381,369121382,369121385,369121386,369121389,369121390,
	                         369121605,369121606,369121609,369121610,369121613,369121614,
							 369125861,369125862,369125866,369125869,369125870};  
    float dcu_fractions_2012[N_2012]={0.346, 0.329, 0.324, 0.339, 0.345, 0.329,
                                      0.345, 0.339, 0.339, 0.312, 0.320, 0.345,
								      0.400, 0.307, 0.217, 0.627, 0.476};
									  
    for(int idet=0; idet<N_2012; idet++)
	  if(detid==detids_2012[idet]) scale = dcu_fractions_2012[idet];

	cout<<"dcu/ps fraction = "<<scale<<endl;

    for(int ic=0; ic<gcur_DCU->GetN(); ic++)
    {
      gcur_DCU->GetPoint(ic, time, current);
      gcur_DCU->SetPoint(ic, time, current*scale);
    }
  }
  
}

TGraph* GetVoltageDrop(TGraph* gsteps, TGraph* gcur_DCU, TGraph* gcur_PS, bool current_ratio=false)
{
  
  // if current_ratio is true, return ratio of leakage current DCU / PS
  
  TGraphErrors* gshift = new TGraphErrors();
  TGraphErrors* gratio = new TGraphErrors();
  
  double time=0;
  double current=0;
  double voltage=0;
  double previous_voltage=300;
  int time_cur=0;
  int time_volt=0;
  int prev_time_volt=0;
  
  double current_PS=0;
  double previous_current_PS=0;
  int time_cur_PS=0;
  int prev_time_cur_PS=0;
  
  double cur_PS=0;
  double cur_PS_err=0;
  double volt=0;
  double volt_err=0;
  double shift=0;
  double shift_err=0;
  
  
  // Loop on DCU measurements
  int ipt=0;
  for(int ic=0; ic<gcur_DCU->GetN(); ic++)
  {
    gcur_DCU->GetPoint(ic, time, current);
    time_cur=int(time);
    if(time_cur==0) continue;
    // cout<<time_cur<<" "<<current<<"uA"<<endl;
    
    // Look for current in PS at time_cur
    for(int icps=0; icps<gcur_PS->GetN(); icps++)
    {
      gcur_PS->GetPoint(icps, time, current_PS);
      time_cur_PS=int(time);
      
      // Get PS current
      if(time_cur>=prev_time_cur_PS && time_cur<time_cur_PS)
      {
        // cout<<prev_time_cur_PS<<" "<<previous_current_PS<<"uA - "<<time_cur_PS<<" "<<current_PS<<"uA"<<endl;
        cur_PS=(previous_current_PS+current_PS)/2.;
        if(!prev_time_cur_PS) cur_PS=0; // time_cur before first PS meas. Set to 0 to not use it.
        cur_PS_err=fabs(current_PS-cur_PS);
        if(cur_PS_err<2) cur_PS_err=2;
      }
      
      previous_current_PS=current_PS;
      prev_time_cur_PS=time_cur_PS;
    }
    
    
    // Look for voltage step corresponding at time_cur
    previous_voltage=300;
    prev_time_volt=0;
    for(int iv=0; iv<gsteps->GetN(); iv++)
    {
      gsteps->GetPoint(iv, time, voltage);
      time_volt=int(time);
      // cout<<"  "<<time_volt<<" "<<voltage<<"V"<<endl;
      
      // Compute voltage shift
      if(time_cur>=prev_time_volt && time_cur<time_volt)
      {
        //cout<<prev_time_volt<<" "<<previous_voltage<<"V - "<<time_volt<<" "<<voltage<<"V"<<endl;
        volt=(previous_voltage+voltage)/2.;
        volt_err=fabs(voltage-volt);
        if(volt_err<1) volt_err=1;

        /*if(cur_PS_err<2.5 && fabs(current/cur_PS-0.333)>0.05) 
        {
          std::cout<<" Currents : DCU "<<current<<" PS "<<cur_PS<<"  detid "<<detid<<" "<<voltage<<"V"<<std::endl;  
        }*/
		
        if(cur_PS && prev_time_volt!=0)
        // remove points not in a voltage step
		if(volt_err<=10) //1
		{
          shift= current*13.8e-3+cur_PS*1.e-3;
          shift_err= 2*13.8e-3+cur_PS_err*1.e-3;
          gshift->SetPoint(ipt, volt, shift);
          gshift->SetPointError(ipt, volt_err, shift_err);
		  gratio->SetPoint(ipt, volt, current/cur_PS);
		  gratio->SetPointError(ipt, volt_err, 
		    current/cur_PS*sqrt(pow(cur_PS_err/cur_PS,2) + pow(2/current,2)));
          ipt++;
		}
      }
      
      previous_voltage=voltage;
      prev_time_volt=time_volt;
    }
  } // End of loop on DCU measurements
  
  TH1F* h = gshift->GetHistogram();
  h->GetXaxis()->SetTitle("V_{bias} [V]");
  h->GetYaxis()->SetTitle("V_{drop} [V]");
  
  if(current_ratio) return gratio;
  else return gshift;
  
}


//----------------------------------------------------------------------------

void DrawConditions(char* subdet="TIB", char* run="20120928_run203832", int detid=369121381, char* bad_periods="bad_periods_20120928_203832.txt")
{
    cout<<" DetID "<<detid<<endl;

    // Get currents and voltage
    TGraph* gsteps;
    TGraph* gcur_DCU;
    TGraph* gcur_PS;
    TGraph* gvolt;
    GetConditions(gsteps, gcur_DCU, gcur_PS, gvolt, subdet, run, detid, bad_periods);
	
	double Steps_max = GetMaximum(gsteps);
	double PS_max = GetMaximum(gcur_PS);
	double scale_PS = 0.9*Steps_max/PS_max;
	Scale( gcur_PS, scale_PS);
	double DCU_max = GetMaximum(gcur_DCU);
	double scale_DCU = 0.45*Steps_max/DCU_max;
	Scale( gcur_DCU, scale_DCU);
    
    // Draw conditions for monitoring
    TCanvas* c1 = new TCanvas();
    gsteps->Draw("APL");
    if(gvolt) gvolt->Draw("P");
    gcur_DCU->Draw("P");
    gcur_PS->Draw("PL");
    TH1F* h = gsteps->GetHistogram();
    h->SetTitle(Form("DetID %i", detid));
    h->GetYaxis()->SetTitle("[V or #muA]");
    h->GetXaxis()->SetTitle("time");
	
	TLegend *leg = new TLegend(0.70, 0.80, 0.98, 0.98);
	leg->AddEntry(gsteps, "Voltage", "l");
	leg->AddEntry(gcur_PS, Form("Leak. cur. from PS (x%.2f)", scale_PS), "pl");
	leg->AddEntry(gcur_DCU, Form("Leak. cur. from DCU (x%.2f)", scale_DCU), "p");
	leg->Draw();
   
   double x, y;
   for(int ipt=0; ipt<gsteps->GetN(); ipt++)
   {
     gsteps->GetPoint(ipt, x, y);
     TLine *l = new TLine(x, 0, x, 360);
     l->SetLineStyle(3);
     //l->Draw();
   } 
}

double DrawDCUOverPSRatio(char* subdet="TIB_L1", char* run="20120506_run193541", int detid=369121381, char* bad_periods="", bool print=false)
{
    cout<<" DetID "<<detid<<endl;

    // Get currents and voltage
    TGraph* gsteps;
    TGraph* gcur_DCU;
    TGraph* gcur_PS;
    TGraph* gvolt;
    GetConditions(gsteps, gcur_DCU, gcur_PS, gvolt, subdet, run, detid, bad_periods);
	
    TGraph* gratio = GetVoltageDrop(gsteps, gcur_DCU, gcur_PS, true);
	
	double mean_ratio = 0.;
	gratio->Fit("pol0");
	mean_ratio = gratio->GetFunction("pol0")->GetParameter(0);

    // Draw conditions for monitoring
    TCanvas* c1 = new TCanvas();
	TH1F* h = gratio->GetHistogram();
    h->SetTitle(Form("DetID %i", detid));
	h->GetYaxis()->SetTitle("I_{leak} ratio DCU/PS");
	h->GetYaxis()->SetTitleOffset(1.2);
	h->GetXaxis()->SetTitle("Voltage [V]");
	gratio->Draw("AP");
	if(print) c1->Print(Form("ILeakDCUOverPS_detid_%i.eps", detid));
	c1->Modified();
	c1->Update();

	
	return mean_ratio;
}

//----------------------------------------------------------------------------

Double_t fitvdrop(Double_t *x, Double_t *par)
{
  Double_t value, value_th, value_0;
  if (x[0]<par[2])
  {
  	value = par[0] + par[1]*sqrt(x[0]) + par[3]*x[0]*sqrt(x[0]); // x^1/2 + x^3/2
  }
  else
  {
  	value_th = par[0] + par[1]*sqrt(par[2]) + par[3]*par[2]*sqrt(par[2]);
  	value_0 = value_th - par[4]*par[2];
  	value = value_0 + par[4]*x[0]; // pol1
  }
  //else { value=0; TF1::RejectPoint();}
  
  return value;
}

void ComputeCorrections(char* subdet, char* run, int* detids, const int N, char* bad_periods="")
{
  
  // Histos and output file
  TFile* fout = new TFile(Form("LeakCurCorr_%s_%s.root", subdet, run),"recreate");
  TH1F* hchi2 = new TH1F("hchi2", "Chi2/NDF", 100, 0, 50);
  TGraph* g2param = new TGraph();
  TH1F* hparam0 = new TH1F("hparam0", "param0", 100, -1, 1);
  TH1F* hparam1 = new TH1F("hparam1", "param1", 500, -5, 5);
  TH1F* hparam2 = new TH1F("hparam2", "param2", 100, -10, 10);

  // Loop over modules
  int detid=0;
  int ifit=0;
  for(int idet=0; idet<N; idet++)
  {
    detid = detids[idet];
    cout<<" DetID "<<detid<<endl;

    // Get currents and voltage
    TGraph* gsteps;
    TGraph* gcur_DCU;
    TGraph* gcur_PS;
    TGraph* gvolt;
    GetConditions(gsteps, gcur_DCU, gcur_PS, gvolt, subdet, run, detid, bad_periods);
    
	// scale graphs for drawing
	double Steps_max = GetMaximum(gsteps);
	double PS_max = GetMaximum(gcur_PS);
	double scale_PS = 0.9*Steps_max/PS_max;
	TGraph* gcur_PS_clone = (TGraph*) gcur_PS->Clone();
	Scale( gcur_PS_clone, scale_PS);
	double DCU_max = GetMaximum(gcur_DCU);
	double scale_DCU = 0.45*Steps_max/DCU_max;
	TGraph* gcur_DCU_clone = (TGraph*) gcur_DCU->Clone();
	Scale( gcur_DCU_clone, scale_DCU);
    
    // Draw conditions for monitoring
    TCanvas* c1 = new TCanvas();
    gsteps->Draw("APL");
    if(gvolt) gvolt->Draw("P");
    gcur_DCU_clone->Draw("P");
    gcur_PS_clone->Draw("PL");
	TH1F* h = gsteps->GetHistogram();
	h->SetTitle(Form("DetID %i", detid));
	h->GetYaxis()->SetTitle("[V or #muA]");
	h->GetXaxis()->SetTitle("time");

	TLegend *leg = new TLegend(0.70, 0.80, 0.98, 0.98);
	leg->AddEntry(gsteps, "Voltage", "l");
	leg->AddEntry(gcur_PS_clone, Form("Leak. cur. from PS (x%.2f)", scale_PS), "pl");
	leg->AddEntry(gcur_DCU_clone, Form("Leak. cur. from DCU (x%.2f)", scale_DCU), "p");
	leg->Draw();
    
    // Compute voltage drop induced by leakage current
    TGraph* gvdrop = GetVoltageDrop(gsteps, gcur_DCU, gcur_PS);
    TCanvas* c2 = new TCanvas("c2", "V drop", 200, 0, 700, 500);
    gvdrop->SetTitle(Form("DetID %i", detid));
    gvdrop->SetMarkerStyle(20);
    gvdrop->Draw("AP");
    
    // Fit voltage drop
    // sqrt(x) function works for spring 2012 runs for modules with high Vdepl (almost no 'plateau')

    TF1* fvdrop = new TF1("fvdrop", fitvdrop , 30, 360, 5);
    fvdrop->SetParameter(0, -0.2);
    fvdrop->SetParLimits(0, -10, 0);
    fvdrop->SetParameter(1,-1);
    fvdrop->SetParLimits(1, 0, 100);
    fvdrop->SetParameter(2, 200);
    //fvdrop->SetParLimits(2, 175, 225);
    //fvdrop->SetParLimits(3, -1, 10);  
    gvdrop->Fit("fvdrop", "R");
    cout << "Chi_Square_fit = " << fvdrop->GetChisquare()/fvdrop->GetNDF() << endl;
    
        
    TF1* fit = (TF1*) gvdrop->GetListOfFunctions()->First();
    
    // Store fit result
    if(fit)
    if(fit->GetNDF()>1)
    if(fit->GetChisquare()/fit->GetNDF() < 50.)
    {
      fout->cd();
      cout<<"Storing fit for detid "<<detid<<endl;
      fit->SetName(Form("fit_%i", detid));
      fit->Write();
      gvdrop->SetName(Form("vdrop_%i", detid));
      gvdrop->Write();
      if(fit->GetNDF()) hchi2->Fill(fit->GetChisquare()/fit->GetNDF());
      g2param->SetPoint(ifit, fit->GetParameter(0), fit->GetParameter(1));
      hparam0->Fill(fit->GetParameter(0));
      hparam1->Fill(fit->GetParameter(1));
      hparam2->Fill(fit->GetParameter(2));
      ifit++;
    }
    
    c1->Modified();
    c1->Update();
    //c1->Print(Form("Ileak-Vbias_%s_%i.pdf", run, detid));
    c2->Modified();
    c2->Update();
    //c2->Print(Form("IleakEffect_%s_%i.pdf", run, detid));
    getchar();
    
    delete c1;
    delete c2;
  
  } // End of loop over modules
  

  // Write output root file
  cout<<"Chi2 mean : "<<hchi2->GetMean()<<endl;
  cout<<"Correl param : "<<g2param->GetCorrelationFactor()<<endl;
  cout<<"param0 : "<<hparam0->GetMean()<<" +/- "<<hparam0->GetRMS()<<endl;
  cout<<"param1 : "<<hparam1->GetMean()<<" +/- "<<hparam1->GetRMS()<<endl;
  cout<<"param2 : "<<hparam2->GetMean()<<" +/- "<<hparam2->GetRMS()<<endl;
  fout->cd();
  hchi2->Write();
  g2param->SetName("g2param");
  g2param->Write();
  hparam0->Write();
  hparam1->Write();
  hparam2->Write();
  fout->Close();
  
}

//------------------------------------------------------------------------

void ComputeILeakDCUFractions(char* subdet="TIB_L1", char* run="20120506_run193541")
{

  // 2012
  const int N_2012=17;
  int detids_2012[N_2012]={369121381,369121382,369121385,369121386,369121389,369121390,369121605,369121606,369121609,369121610,
                           369121613,369121614,369125861,369125862,369125866,369125869,369125870};  

  
  // Loop over modules
  int detid=0;
  int ifit=0;
  for(int idet=0; idet<N_2012; idet++)
  {
    detid = detids_2012[idet];
    cout<<" DetID "<<detid<<endl;
	
    DrawDCUOverPSRatio(subdet, run, detid);
	getchar();
  }
  
  // 20120506_run193541
  float dcu_fractions_2012[N_2012]={0.346, 0.329, 0.324, 0.339, 0.345, 0.329,
                                    0.345, 0.339, 0.339, 0.312, 0.320, 0.345,
								    0.400, 0.307, 0.217, 0.627, 0.476};
}


//void ComputeILeakCorrections(char* subdet="TIB_L1", char* run="20120506_run193541")
void ComputeILeakCorrections()
{

  // 2011
  const int N_2011=11;
  int detids_2011[N_2011]={369121605,369121606,369121609,369121610,369121613,369121614,369125861,369125862,369125866,369125869,369125870};
  // Not enough data //ComputeCorrections("TIB_L1", "20111012_run178367", detids_2011, N_2011);

  // 2012
  const int N_2012=17;
  int detids_2012[N_2012]={369121381,369121382,369121385,369121386,369121389,369121390,369121605,369121606,369121609,369121610,
                           369121613,369121614,369125861,369125862,369125866,369125869,369125870};  
  //ComputeCorrections("TIB_L1", "20120405_run190459", detids_2012, N_2012, "Steps/bad_periods_20120405_run190459.txt");
  //ComputeCorrections("TIB_L1", "20120506_run193541", detids_2012, N_2012);
  // only PS currents
  //ComputeCorrections("TIB_L1", "20120510_run193928", detids_2012, N_2012, "Steps/bad_periods_20120510_run193928.txt");

  // L2 2012
  const int N_L2_2012=6;
  int detids_L2_2012[N_L2_2012]={369136677, 369136678, 369136681, 369136682, 369136685, 369136686};
  //ComputeCorrections("TIB_L2", "20120405_run190459", detids_L2_2012, N_L2_2012);
  // Not enought data //ComputeCorrections("TIB_L2", "20120506_run193541", detids_L2_2012, N_L2_2012);
  
  // August 2012 , only PS currents (L1+L2)
  const int N=4;
  int detids[N]={369125861, 369125862, 369136681, 369136682};
  //ComputeCorrections("TIB_L1", "20120812_run200786", detids, N);

  // From Summer 2012
  const int N_2012_bis=11;
  int detids_2012_bis[N_2012_bis]={369121381,369121382,369121385,369121386,369121389,369121390,369125861,369125862,369125866,369125869,369125870};
  int detids_2012_4PS[4]={369121381,369121382,369125861,369125862};

  //ComputeCorrections("TIB_L1", "20120728_run199832", detids_2012, N_2012, "Steps/bad_periods_20120728_run199832.txt");
  //ComputeCorrections("TIB_L1", "20120812_run200786", detids_2012, N_2012);

  //ComputeCorrections("TIB_L1", "20120928_run203832", detids_2012_4PS, 4, "Steps/bad_periods_20120928_run203832.txt");   
  //ComputeCorrections("TIB_L1", "20120928_run203832", detids_2012, N_2012, "Steps/bad_periods_20120928_run203832.txt");

  //ComputeCorrections("TIB_L1", "20121130_run208339", detids_2012_bis, N_2012_bis, "Steps/bad_periods_20121130_run208339.txt");
  //ComputeCorrections("TIB_L1", "20130213_run211797", detids_2012_bis, N_2012_bis, "");
  
  
  // TOB
  const int N_TOB=6;
  int detids_TOB[N_TOB]={436281508, 436281512, 436281516, 436281520, 436281524, 436281528};
  //DrawConditions("TOB", "20121130_run208339",436281508 , "Steps/bad_periods_20121130_run208339.txt");
//  ComputeCorrections("TOB", "20110405_run190459", detids_TOB, N_TOB, "Steps/bad_periods_20120405_run190459.txt"); //seg fault because of the inconsistancy between  20110405 and 20120405

  ComputeCorrections("TOB", "20120510_run193928", detids_TOB, N_TOB, "Steps/bad_periods_20120510_run193928.txt");  
  ComputeCorrections("TOB", "20120728_run199832", detids_TOB, N_TOB, "Steps/bad_periods_20120728_run199832.txt");
  ComputeCorrections("TOB", "20120928_run203832", detids_TOB, N_TOB, "Steps/bad_periods_20120928_run203832.txt");    
  ComputeCorrections("TOB", "20121130_run208339", detids_TOB, N_TOB, "Steps/bad_periods_20121130_run208339.txt");

}
