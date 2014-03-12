#include "ComputeILeakCorrections.C"
#include "GetKink.C"

TGraph* GetIleakVsVbias(TGraph* gsteps, TGraph* gcur_DCU, TGraph* gcur_PS)
{
  
  TGraphErrors* gout = new TGraphErrors();
  
  double time=0;
  double current=0;
  double voltage=0;
  double previous_voltage=300;
  int time_cur=0;
  int time_volt=0;
  int prev_time_volt=0;
  
  double volt=0;
  double volt_err=0;
  
  // Loop on DCU measurements
  int ipt=0;
  for(int ic=0; ic<gcur_DCU->GetN(); ic++)
  {
    gcur_DCU->GetPoint(ic, time, current);
    //current*=current;
    time_cur=int(time);
    if(time_cur==0) continue;
     //cout<<time_cur<<" "<<current<<"uA"<<endl;
    //if(time_cur>1333629950 && time_cur<1333647468) continue; // between 2 parts of 2012 april scan
       
    // Look for voltage step corresponding at time_cur
    previous_voltage=300;
    prev_time_volt=0;
    for(int iv=0; iv<gsteps->GetN(); iv++)
    {
      gsteps->GetPoint(iv, time, voltage);
      time_volt=int(time);
       //cout<<"  "<<time_volt<<" "<<voltage<<"V"<<endl;
      
      // Fill graph Current vs Vbias
      if(time_cur>=prev_time_volt && time_cur<time_volt)
      {
         //cout<<prev_time_volt<<" "<<previous_voltage<<"V - "<<time_volt<<" "<<voltage<<"V"<<endl;
		volt=(previous_voltage+voltage)/2.;
        volt_err=fabs(voltage-volt);
        if(volt_err<1) volt_err=1;
      }
      
      previous_voltage=voltage;
      prev_time_volt=time_volt;
    }
	
	
    if(volt_err>1) continue; // remove points not in a voltage step
    gout->SetPoint(ipt, volt, current);
    gout->SetPointError(ipt, volt_err, 2);
    ipt++;
	
  } // End of loop on DCU measurements
  
  TH1F* h = gout->GetHistogram();
  h->GetXaxis()->SetTitle("V_{bias} [V]");
  h->GetYaxis()->SetTitle("I_{leak} [#mu\A]");
  return gout;
  
}

TGraph* AverageIleakVsVbias(TGraph* g)
{
  TGraph* gout = new TGraph();
  
  double volt, current;
  double previous_volt=0;
  double mean_current=0;
  int nmean=0;
  int ip=0;
  
  for(int ipt=0; ipt<g->GetN(); ipt++)
  {
     g->GetPoint(ipt, volt, current);
     if(previous_volt==volt)
	 {
	   mean_current+=current;
	   nmean++;
	 }
	 else
	 {
	   if(previous_volt!=0)
	   {
	     mean_current/=nmean;
	     gout->SetPoint(ip, previous_volt, mean_current);
	     ip++;
	   }
	   mean_current=current;
	   nmean=1;
	 }
	 previous_volt=volt;
  }

  return gout;
}

TGraphErrors* GetDerivative(TGraph* g)
{
  
  // Initialization
  TGraphErrors* gderivative = new TGraphErrors();
  double deriv=0;
  double volt=0;
  double voltage=0;
  double previous_voltage=0;
  double shift=0;
  double previous_shift=0;
  int ipt=0;
  
  // Loop over points
  for(int ip=0; ip<g->GetN(); ip++)
  {
    previous_voltage = voltage;
    previous_shift = shift;
    g->GetPoint(ip, voltage, shift);
    if(!ip) continue;
    if(ip && voltage==previous_voltage) continue;
    volt=(previous_voltage+voltage)/2.;
    deriv = (shift-previous_shift)/(voltage-previous_voltage);
    //cout<<ipt<<" "<<volt<<" "<<deriv<<endl;
    gderivative->SetPoint(ipt, volt, deriv);
    ipt++;
  }
  
  return gderivative;
  
}

//------------------------------------------------------------------------------

void Fit(char* subdet, char* run, int* detids, const int N, char* bad_periods="")
{
  
  // Histos and output file
  TFile* fout = new TFile(Form("LeakageCurrent_%s_%s.root", subdet, run),"recreate");
  TH1F* hchi2 = new TH1F("hchi2", "Chi2/NDF", 100, 0, 5);
  TH1F* hparam0 = new TH1F("hparam0", "param0", 100, 0, 50);
  TH1F* hparam1 = new TH1F("hparam1", "param1", 500, -5, 5);
  TH1F* hparam2 = new TH1F("hparam2", "param2", 100, 0, 400);
  TH1F* hparam3 = new TH1F("hparam2", "param2", 500, -5, 5);
  
  TFile* output = new TFile(Form("DECO_Ileak_%s_%s.root", subdet, run),"recreate");
  TTree* tout = new TTree("tout","SignalSummary");

  ULong64_t odetid;
  int olayer;
  double odepvolt;
  double oerrdepvolt;
  double oplateau;
  double ofitchisquare;
  int ofitstatus;
  double olastpty;
  double ochi2;

  tout->Branch("DETID",&odetid,"DETID/l");
  tout->Branch("LAYER",&olayer,"LAYER/I");
  tout->Branch("ERRDEPVOLT",&oerrdepvolt,"ERRDEPVOLT/D");
  tout->Branch("DEPVOLT",&odepvolt,"DEPVOLT/D");
  tout->Branch("PLATEAU",&oplateau,"PLATEAU/D");
  tout->Branch("FITCHI2",&ofitchisquare,"FITCHI2/D");
  tout->Branch("FITSTATUS",&ofitstatus,"FITSTATUS/I");
  tout->Branch("LASTPOINTS",&olastpty,"LASTPOINTS/D");
  tout->Branch("CHI2",&ochi2,"CHI2/D");

  Int_t kSubdetOffset = 25;
  Int_t subdet_i = 0;
  unsigned layerStartBit_ = 14;
  unsigned layerMask_ = 0x7;
  
  double lastpty = 0;


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
    // Only gsteps and gcur_DCU really needed
    
    // Draw conditions for monitoring
    TCanvas* c1 = new TCanvas();
    gsteps->Draw("APL");
    if(gvolt) gvolt->Draw("P");
    gcur_DCU->Draw("P");
    gcur_PS->Draw("PL");
    
    
    // Get leakage current vs Vbias
    TGraph* gIleak = GetIleakVsVbias(gsteps, gcur_DCU, gcur_PS);
    TCanvas* c2 = new TCanvas("c2", "Ileak", 200, 0, 700, 500);
	TGraph* gIleak = AverageIleakVsVbias(gIleak);
    gIleak->SetMarkerStyle(20);
    gIleak->Draw("AP");
	
	int npt = gIleak->GetN();
	if(npt<3) continue;
	double x,y;
	lastpty = 0;
	gIleak->GetPoint(npt-1,x,y);
	lastpty += y;
	gIleak->GetPoint(npt-2,x,y);
	lastpty += y;
	gIleak->GetPoint(npt-3,x,y);
	lastpty += y;
	lastpty /= 3.;
	
    
    // Fit leakage current
    // function with curve in 2 parts : x^1/2 + x^3/2 and pol1 for 'plateau'
    // 5 parameters, so need more points
    /*TF1* fvdrop = new TF1("fvdrop", fitfunction, 20, 360, 5);
    fvdrop->SetParameter(2, 150);
    fvdrop->SetParLimits(1, 0.001, 10);
    fvdrop->SetParameter(3, 0.01);
    gIleak->Fit("fvdrop");*/
    
    // function with curve in 2 parts : sigmoid and pol1
    // fit needs a lot of points to work
    /*TF1* fvdrop = new TF1("fvdrop", fitfunction2, 20, 360, 5);
     fvdrop->SetParameter(0, 200.);
     fvdrop->SetParameter(1, 0.02);
     fvdrop->SetParameter(2, 150);
     fvdrop->SetParameter(4, 0.1);
     gIleak->Fit("fvdrop");*/
    
    
    // Compute derivative of current vs voltage curve
    TGraphErrors* gderivative = GetDerivative( gIleak );
    gderivative->SetMarkerStyle(20);
    gderivative->Draw("AP");
	 
    // function with curve in 2 parts: exp(-x) and pol1
    /*TF1* fderiv = new TF1("fderiv", fitfunctionderiv, 30, 360, 5);
    fderiv->SetParameter(0, 0.1);
    fderiv->SetParameter(1, 1.);
    fderiv->SetParameter(2, 250);
    if(idet<2) fderiv->SetRange(200, 360);
    gderivative->Fit("fderiv", "R");*/
    
    // function with curve in 2 parts: 2 pol1
    /*TF1* fderiv = new TF1("fderiv", fitfunctionderiv2, 30, 360, 4);
    fderiv->SetParameter(0, 3);
    fderiv->SetParameter(1, -0.01);
    fderiv->SetParameter(2, 250);
    fderiv->SetParameter(3, -0.001);
    if(idet<2) fderiv->SetRange(200, 360);
    gderivative->Fit("fderiv", "R");
    
    TF1* fit = (TF1*) gderivative->GetListOfFunctions()->First();*/
    TF1* fit = 0;

    // Find Kink
    TCanvas* c3 = new TCanvas("c3", "", 400, 0, 700, 500);
    TGraphErrors* gmedian = MedianFilter( gIleak );
    int nfilt=1;
    while (!IsMonoton(gmedian) && nfilt<4) {
      gmedian = MedianFilter( gmedian );
      nfilt++;
    }
	//cout<<nfilt<<" median filter applied"<<endl;
    gmedian = HanningFilter(gmedian);
    gIleak->SetTitle(Form("Detid %i", detid));
    gIleak->Draw("AP");

    TGraphErrors* gscurv = GetCurvatureGraph( gmedian );//GetDerivative(gmedian);
	//gscurv = GetDerivative(gscurv);
    gscurv = HanningFilter(gscurv);
	gscurv->SetMarkerStyle(20);
	TGraph* g3pts = new TGraph();
	double *chi2 = new double;
	double *ndf = new double;
	double *n_rms = new double;
	double *err = new double;
    float xopt = GetOptimalMinNPts(gscurv, g3pts, 0, chi2, ndf, n_rms, err);

    TCanvas* c4 = new TCanvas("c4", "", 500, 0, 700, 500);
	gscurv->Draw("AP");
	g3pts->Draw("P");
	g3pts->Fit("pol2", "q");
	
    float ymin = gscurv->GetYaxis()->GetXmin();
    float ymax = gscurv->GetYaxis()->GetXmax();

    TLine* lopt = new TLine(xopt, ymin, xopt, ymax);
    lopt->SetLineStyle(1);
    lopt->SetLineColor(1);
    lopt->Draw();
    cout<<endl<<"Kink : "<<xopt<<" V"<<endl;
	
	c3->cd();
	c3->Modified();
	TLine* lopt2 = new TLine(xopt, gmedian->GetYaxis()->GetXmin(),
	                          xopt, gmedian->GetYaxis()->GetXmax());
    lopt2->Draw();

    //c3->Print(Form("IleakVsVbias_detid_%i.eps", detid));
    //c4->Print(Form("IleakVsVbias_curv_detid_%i.eps", detid));
	
    // Store fit result
    if(fit)
    {
	  if(fit->GetNDF()>1)
        if(fit->GetChisquare()/fit->GetNDF() < 5.)
        {
          fout->cd();
          cout<<"Vdepl: "<<fit->GetParameter(2)<<" V"<<endl;
          cout<<"Storing fit for detid "<<detid<<endl;
          fit->SetName(Form("fit_%i", detid));
          fit->Write();
          gIleak->SetName(Form("Ileak%i", detid));
          gIleak->Write();
          if(fit->GetNDF()) hchi2->Fill(fit->GetChisquare()/fit->GetNDF());
          hparam0->Fill(fit->GetParameter(0));
          hparam1->Fill(fit->GetParameter(1));
          hparam2->Fill(fit->GetParameter(2));
          hparam3->Fill(fit->GetParameter(3));
          ifit++;
        }
    }
	else
	{
	  odetid = detid;
      
	  Int_t detid_forlayer = detid;
      if(subdet=="TOB") detid_forlayer = detid/10;
	  subdet_i = (detid_forlayer>>kSubdetOffset)&0x7;
	  // For TIB and TOB
      if(subdet_i==3 || subdet_i==5) olayer = ((detid_forlayer>>layerStartBit_) & layerMask_);

      odepvolt = xopt;
	  if(err) oerrdepvolt = *err;
	  else oerrdepvolt = 0;
	  
      oplateau = 0;
      ofitchisquare = *chi2;
      ofitstatus = *ndf;
      olastpty = lastpty;
      ochi2 = n_rms;
	  
	  tout->Fill();
	
	}
	
	
    c1->Modified();
    c1->Update();
    //c1->Print(Form("Ileak-Vbias_%s_%i.pdf", run, detid));
    c2->Modified();
    c2->Update();
    //c2->Print(Form("IleakEffect_%s_%i.pdf", run, detid));
    c3->Modified();
    c3->Update();
    c4->Modified();
    c4->Update();
    getchar();
    
    delete c1;
    delete c2;
    delete c3;
    delete c4;
    
  } // End of loop over modules
  
  
  // Write output root file
  cout<<"Chi2/ndf mean : "<<hchi2->GetMean()<<endl;
  cout<<"param0 : "<<hparam0->GetMean()<<" +/- "<<hparam0->GetRMS()<<endl;
  cout<<"param1 : "<<hparam1->GetMean()<<" +/- "<<hparam1->GetRMS()<<endl;
  cout<<"param2 : "<<hparam2->GetMean()<<" +/- "<<hparam2->GetRMS()<<endl;
  cout<<"param3 : "<<hparam3->GetMean()<<" +/- "<<hparam3->GetRMS()<<endl;
  fout->cd();
  hchi2->Write();
  hparam0->Write();
  hparam1->Write();
  hparam2->Write();
  hparam3->Write();
  fout->Close();
  
  output->cd();
  tout->Write();
  //gDirectory->Write();
  output->Close();
}

//------------------------------------------------------------------------------

void FitLeakageCurrent(char* subdet="TIB_L1", char* run="20120506_run193541")
{
   
  // Many points needed, works well only for august 2012 scan
  
  // 2012
  const int N_2012=17;
  int detids_2012[N_2012]={369121381,369121382,369121385,369121386,369121389,369121390,369121605,369121606,369121609,369121610,
  369121613,369121614,369125861,369125862,369125866,369125869,369125870};  
  //Fit("TIB_L1", "20120405_run190459", detids_2012, N_2012);
  //Fit("TIB_L1", "20120506_run193541", detids_2012, N_2012);

  // L2 2012
  const int N_L2_2012=6;
  int detids_L2_2012[N_L2_2012]={369136677, 369136678, 369136681, 369136682, 369136685, 369136686};
  //Fit("TIB_L2", "20120405_run190459", detids_L2_2012, N_L2_2012);
  
  // August 2012 , only PS currents (L1+L2)
  const int N=4;
  int detids[N]={369125861, 369125862, 369136681, 369136682};
  //Fit("TIB", "L1_20120812_run200786", detids, N);  
  
  const int N_1PS=6;
  int detids_1PS[N_1PS]={369121605,369121606,369121609,369121610,369121613,369121614};  
  Fit("TIB_L1", "20120928_run203832", detids_1PS, N_1PS, "Steps/bad_periods_20120928_run203832.txt");
  //Fit("TIB_L1", "20120928_run203832", detids_2012, N_2012, "Steps/bad_periods_20120928_run203832.txt");
}
