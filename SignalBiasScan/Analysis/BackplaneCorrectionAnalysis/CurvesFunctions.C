#include <string>

int GetSubdet(ULong64_t modid)
{
  if(modid>999999999) modid/=10; // sensor option for TOB
  // 3->6 : TIB, TID, TOB, TEC
  static const Int_t kSubdetOffset = 25;
  int subdet = (modid>>kSubdetOffset)&0x7;

  return subdet;
}

char* GetSubdetString(ULong64_t modid)
{
  char* subdet;
  sprintf(subdet, "");
  int isubdet = GetSubdet(modid);
  if(isubdet==3) sprintf(subdet, "TIB");
  if(isubdet==4) sprintf(subdet, "TID");
  if(isubdet==5) sprintf(subdet, "TOB");
  if(isubdet==6) sprintf(subdet, "TEC");

  return subdet;
}
int GetLayer(ULong64_t modid)
{
  if(modid>999999999) modid/=10; // sensor option for TOB
  Int_t subdet = GetSubdet(modid);
  
  static const unsigned layerStartBit = 14;
  static const unsigned layerMask = 0x7;
  Int_t layer = 0;
  
  // For TIB and TOB
  if(subdet==3 || subdet==5) layer = ((modid>>layerStartBit) & layerMask);

  return layer;

}

void GetFitParams(string subdet, string run, ULong64_t modid, double & param0, double & param1, double & param2)
{ 

  ULong64_t detid_1;
  int layer_1;
  double depvolt_1=0;
  double errdepvolt_1;
  double plateau_1=0;
  double fitchi2_1;
  int fitstatus_1;
  double chi2_1;

  char* file = Form("DECO_Vdeplchi_%s%s.root", subdet.c_str(), run.c_str());
  TFile* f = TFile::Open(file);
  if(!f) {cout<<"Error : no file '"<<file<<"'"<<endl; param0=0; param1=0; param2=0; return;}

  TTree *T1  = (TTree*)f->Get("tout");

  T1->SetBranchAddress("DETID",&detid_1);
  //T1->SetBranchAddress("LAYER",&layer_1);
  T1->SetBranchAddress("DEPVOLT",&depvolt_1);
  T1->SetBranchAddress("ERRDEPVOLT",&errdepvolt_1);
  T1->SetBranchAddress("PLATEAU",&plateau_1);
  T1->SetBranchAddress("FITCHI2",&fitchi2_1);
  T1->SetBranchAddress("FITSTATUS",&fitstatus_1);
  //T1->SetBranchAddress("CHI2",&chi2_1);
  unsigned int nentries1 = T1->GetEntries();

  for(unsigned int j=0; j<nentries1; j++)
  {
	T1->GetEntry(j);
	if(detid_1==modid) break;  
  }


  //double params[2];
  param0= depvolt_1;
  param1= plateau_1;
  //cout<<"Fit : plateau "<<plateau_1<<"  depvolt "<<depvolt_1<<endl; 
 
  if(detid_1!=modid) cout<<"Error : no function in '"<<modid<<"' in "<<run<<endl;
  f->Close();

  return;
}

void GetFitInfos(string subdet, string run, ULong64_t modid, double & info0, double & info1)
{ 

  ULong64_t detid_1;
  int layer_1;
  double depvolt_1=0;
  double errdepvolt_1;
  double plateau_1=0;
  double fitchi2_1;
  int fitstatus_1;
  double chi2_1;

  char* file = Form("DECO_Vdeplchi_%s%s.root", subdet.c_str(), run.c_str());
  TFile* f = TFile::Open(file);
  if(!f) {cout<<"Error : no file '"<<file<<"'"<<endl; return;}

  TTree *T1  = (TTree*)f->Get("tout");

  T1->SetBranchAddress("DETID",&detid_1);
  //T1->SetBranchAddress("LAYER",&layer_1);
  T1->SetBranchAddress("DEPVOLT",&depvolt_1);
  T1->SetBranchAddress("ERRDEPVOLT",&errdepvolt_1);
  T1->SetBranchAddress("PLATEAU",&plateau_1);
  T1->SetBranchAddress("FITCHI2",&fitchi2_1);
  T1->SetBranchAddress("FITSTATUS",&fitstatus_1);
  T1->SetBranchAddress("CHI2",&chi2_1);
  unsigned int nentries1 = T1->GetEntries();

  for(unsigned int j=0; j<nentries1; j++)
  {
	T1->GetEntry(j);
	if(detid_1==modid) break;  
  }

  info0= depvolt_1;
  info1= fitchi2_1;
 
  if(detid_1!=modid) cout<<"Error : no function in '"<<modid<<"' in "<<run<<endl;
  f->Close();
  
  return;
}

TGraphErrors* GetGraph(string subdet, string run, ULong64_t modid, int err_type)
{

  TGraphErrors* g;

  //Int_t detid_2;
  ULong64_t detid_2;
  double volt_2;
  int id_2;
  double evolt_2;
  double mpv_2;
  double empv_2;
  int tempdetid_2;
  double chi2overndf_2;

  const int step = 200;
  double avolt[step];
  int aid[step];
  double aevolt[step];
  double ampv[step];
  double aempv[step];

  char* file = Form("%s_output_DecoMode_prof_%s.root", subdet.c_str(), run.c_str());
  TFile* f = TFile::Open(file);
  if(!f) {cout<<"Error : no file '"<<file<<"'"<<endl; return g;}

  TTree* tr = (TTree*) f->FindObjectAny("T");
  tr->SetBranchAddress("DetID",&detid_2);
  tr->SetBranchAddress("Voltage",&volt_2);
  tr->SetBranchAddress("Index",&id_2);
  tr->SetBranchAddress("errVoltage",&evolt_2);
  tr->SetBranchAddress("MPV",&mpv_2);
  tr->SetBranchAddress("errMPV",&empv_2);
  tr->SetBranchAddress("Chi2OverNdf",&chi2overndf_2);

  UInt_t nentries = tr->GetEntries();

  int k=0;
  for(unsigned int j = 0; j <nentries; j++){
    tr->GetEntry(j);

    if(k>0 && detid_2!=modid) break;
	if(detid_2==modid){
	  if(empv_2<5){
	    avolt[k]=volt_2;
	    aid[k]=id_2;
	    aevolt[k]=evolt_2;
	    ampv[k]=mpv_2;
	    aempv[k]=empv_2;
		if(err_type==1) aempv[k]=2.5;
		if(err_type==2) aempv[k]=5;
		if(err_type==3) aempv[k]=empv_2*5;
		if(err_type==4) aempv[k]=empv_2*10;
		if(err_type==5) aempv[k]=10;
	    k++;

	  }
	}
  }
  
  f->Close();
  
  if(!k) return 0;
  
  g = new TGraphErrors(k, avolt, ampv, aevolt, aempv);
  g->SetLineColor(2);
  g->SetMarkerStyle(20);

  return g;
}

TGraphErrors* GetClusterWidthGraph(string subdet, string run, ULong64_t modid, int err_type)
{

  TGraphErrors* g;

  //Int_t detid_2;
  ULong64_t detid_2;
  double volt_2;
  int id_2;
  double evolt_2;
  double mean_2;
  double emean_2;
  int tempdetid_2;

  const int step = 200;
  double avolt[step];
  int aid[step];
  double aevolt[step];
  double amean[step];
  double aemean[step];

  char* file = Form("%s_output_DecoMode_%s.root", subdet.c_str(), run.c_str());
  TFile* f = TFile::Open(file);
  if(!f) {cout<<"Error : no file '"<<file<<"'"<<endl; return g;}

  TTree* tr = (TTree*) f->FindObjectAny("T");
  tr->SetBranchAddress("DetID",&detid_2);
  tr->SetBranchAddress("Voltage",&volt_2);
  tr->SetBranchAddress("Index",&id_2);
  tr->SetBranchAddress("errVoltage",&evolt_2);
  tr->SetBranchAddress("Mean",&mean_2);
  tr->SetBranchAddress("errMean",&emean_2);

  UInt_t nentries = tr->GetEntries();

  int k=0;
  for(unsigned int j = 0; j <nentries; j++){
    tr->GetEntry(j);
    if(k>0 && detid_2!=modid) break;
	if(detid_2==modid){
	  if(emean_2<5){
	    avolt[k]=volt_2;
	    aid[k]=id_2;
	    aevolt[k]=evolt_2;
	    amean[k]=mean_2;
	    aemean[k]=emean_2;
		if(err_type==1) aemean[k]=2.5;
		if(err_type==2) aemean[k]=5;
		if(err_type==3) aemean[k]=emean_2*5;
		if(err_type==4) aemean[k]=emean_2*10;
		if(err_type==5) aemean[k]=10;
	    k++;

	  }
	}
  }
  
  f->Close();
  if(!k) return 0;
  
  g = new TGraphErrors(k, avolt, amean, aevolt, aemean);
  g->SetLineColor(2);
  g->SetMarkerStyle(20);

  return g;
}

TGraphErrors* GetOneOverS2Graph(TGraphErrors* g)
{

  TGraphErrors* g2 = new TGraphErrors();
  double voltage=0;
  double signal=0;
  int ipt=0;
  for(unsigned int ip=0; ip<g->GetN(); ip++)
  {
   g->GetPoint(ip, voltage, signal);
   g2->SetPoint(ipt, voltage, 1./pow(signal,2));
   ipt++;
  }
  g2->SetMarkerStyle(20);
 
  return g2;
}

TGraphErrors* GetDerivativeGraph(TGraphErrors* g)
{

  TGraphErrors* gderivative = new TGraphErrors();
  double voltage=0;
  double previous_voltage=0;
  double signal=0;
  double previous_signal=0;
  int ipt=0;
  double volt=0;
  double deriv=0;
  for(unsigned int ip=0; ip<g->GetN(); ip++)
  {
   previous_voltage = voltage;
   previous_signal = signal;
   g->GetPoint(ip, voltage, signal);
   if(!ip) continue;
   if(ip && voltage==previous_voltage) continue;
   volt=(previous_voltage+voltage)/2.;
   deriv = (signal-previous_signal)/(voltage-previous_voltage);
   gderivative->SetPoint(ipt, volt, deriv);
   ipt++;
  }
  gderivative->SetMarkerStyle(20);
 
  return gderivative;
}

TGraphErrors* GetCurvatureGraph(TGraphErrors* g)
{
  TGraphErrors* gd = GetDerivativeGraph( g );
  TGraphErrors* gd2 = GetDerivativeGraph( gd );
  TGraphErrors* gcurv = new TGraphErrors();
  
  //gd->Print("all");
  //gd2->Print("all");
  
  double x1, y1, x2, y2, curv;
  double xmin=-1;
  double ymin=0;
  double xmax=-1;
  double ymax=0;
  for(int i=0; i<gd2->GetN(); i++)
  {
    gd2->GetPoint(i, x2, y2);
    
    // Need to interpolate to get y'
    int id=0;
    while(x1<=x2 && id<gd->GetN())
    {
      gd->GetPoint(id, x1, y1);
      if(x1<=x2) {xmin=x1; ymin=y1;}
      if(x1>x2) {xmax=x1; ymax=y1;}
      id++;
    }
    // no point with x1<=x2 found
    if(xmin==-1) continue;
    // point with x1==x2 found
    if(xmin==x2) y1=ymin;
    else
      // linear interpolatation
      if(xmin!=-1 && xmax!=-1)
      {
        float mu=(x2-xmin)/(xmax-xmin);
        if(mu<0 || mu>1) continue;
        y1=ymin*(1-mu)+ymax*mu;
      }
    
    // compute curvature
    curv=y2/TMath::Power(1+y1*y1, 3./2);
    gcurv->SetPoint(i, x2, curv);
  }
  
  /*cout<<"------"<<endl;
  g->Print("all");
  gd->Print("all");
  gd2->Print("all");
  gcurv->Print("all");
  cout<<"------"<<endl;*/

  return gcurv;
}


/////////////// LISSAGE //////////////////

TGraphErrors* MedianFilter(TGraphErrors* g)
{
  TGraphErrors* gsmooth = new TGraphErrors();
  
  const int N=3; // FUNCTION WRITTEN FOR N=3, CHECKS NEEDED IF WANT TO GENERALIZE
  double* voltage = new double[N];
  double* signal = new double[N];
  int* index = new int[N];
  int ipt=1;
  double signal_mean=0;
  double voltage_mean=0;
  if(N<2 || N>g->GetN()) return gsmooth;
  
  for(unsigned int ip=0; ip<g->GetN()+1-N; ip++)
  {
    voltage_mean = 0;
    signal_mean = 0;
    for(int i=0; i<N; i++)
    {
      g->GetPoint(ip+i, voltage[i], signal[i]);
    }
    TMath::Sort(N , signal, index, false); // increasing order
    // For N=3 median is the 2nd value
    gsmooth->SetPoint(ipt, voltage[1], signal[index[1]]);
    ipt++;
  }
  
  // For first point
  g->GetPoint(0, voltage[0], signal[0]);
  gsmooth->GetPoint(1, voltage[1], signal[1]);
  gsmooth->GetPoint(2, voltage[2], signal[2]);
  //cout<<signal[2]<<" "<<signal[1]<<" "<<signal[0]<<endl;
  signal[2] = 3*signal[1] - 2*signal[2];
  TMath::Sort(N , signal, index, false);
  gsmooth->SetPoint(0, voltage[0], signal[index[1]]);
  //cout<<signal[2]<<" "<<signal[index[1]]<<endl;
  
  // For last point
  g->GetPoint(g->GetN()+2-N, voltage[0], signal[0]);
  gsmooth->GetPoint(gsmooth->GetN()-1, voltage[1], signal[1]);
  gsmooth->GetPoint(gsmooth->GetN()-2, voltage[2], signal[2]);
  signal[2] = 3*signal[1] - 2*signal[2];
  TMath::Sort(N , signal, index, false);
  gsmooth->SetPoint(ipt, voltage[0], signal[index[1]]);
  
  gsmooth->SetMarkerStyle(20);
  //gsmooth->Print("all");
  
  return gsmooth;  
}

bool IsMonoton(TGraphErrors* g)
{
  bool ismonoton=true;
  bool isincreasing=true;
  double voltage=0;
  double previous_voltage=0;
  double signal=0;
  double previous_signal=0;
  
  int ip=0;
  while(ismonoton && ip<g->GetN())
  {
    previous_signal = signal;
    previous_voltage = voltage;
    g->GetPoint(ip, voltage, signal);
    if(ip==1 && previous_signal>signal) isincreasing=false;
    if(ip>1 && 
       ((isincreasing && previous_signal>signal) ||
        (!isincreasing && previous_signal<signal))) ismonoton=false;
    ip++;
  }
  
  return ismonoton;
}

TGraphErrors* HanningFilter(TGraphErrors* g, int npt=3)
{
  TGraphErrors* gsmooth = new TGraphErrors();
  
  if(npt!=3 && npt!=5) { cout<<"Error : Lissage de Hanning prévu seulement sur 3 ou 5 points."<<endl;; return gsmooth; }
  
  const int N=npt;// FUNCTION WRITTEN FOR N=3, CHECKS NEEDED IF WANT TO GENERALIZE
  // coef for hanning
  double* coef = new double[N];
  if(npt==3)
  {
   coef[0]=1./4;
   coef[1]=1./2;
   coef[2]=1./4;
  }
  if(npt==5)
  {
   coef[0]=1./8;
   coef[1]=1./4;
   coef[2]=1./4;
   coef[3]=1./4;
   coef[4]=1./8;
  }
  
  double voltage=0;
  double signal=0;
  //int ipt=(N-1)/2;
  int ipt=1;
  double signal_mean=0;
  double voltage_mean=0;
  if(N<2 || N>g->GetN()) return gsmooth;
  
  for(int ip=0; ip<g->GetN()-N+1; ip++)
  {
    voltage_mean = 0;
    signal_mean = 0;
    for(int i=0; i<N; i++)
    {
      g->GetPoint(ip+i, voltage, signal);
      voltage_mean+=voltage;
      signal_mean+=signal*coef[i];
    }
    voltage_mean/=N;
    gsmooth->SetPoint(ipt, voltage_mean, signal_mean);
    ipt++;
  }
  //gsmooth->Print("all");
  
  // first point : linear extrapolation of smoothed curve (no slope change)
  double* voltages = new double[N];
  double* signals = new double[N];
  
  //g->Print("all");
  //gsmooth->Print("all");
  
  g->GetPoint(0, voltages[0], signals[0]);
  gsmooth->GetPoint(1, voltages[1], signals[1]);
  gsmooth->GetPoint(2, voltages[2], signals[2]);
  voltage_mean=2*voltages[1]-voltages[2];
  signal_mean=(signals[2]-signals[1])/(voltages[2]-voltages[1])*(voltage_mean-voltages[1])+signals[1];
  gsmooth->SetPoint(0, voltage_mean, signal_mean);
  
  // last point : linear extrapolation of smoothed curve (no slope change)
  g->GetPoint(g->GetN()-N+2, voltages[0], signals[0]);
  gsmooth->GetPoint(g->GetN()-N, voltages[1], signals[1]);
  gsmooth->GetPoint(g->GetN()-N+1, voltages[2], signals[2]);
  voltage_mean=2*voltages[2]-voltages[1];
  signal_mean=(signals[2]-signals[1])/(voltages[2]-voltages[1])*(voltage_mean-voltages[2])+signals[2];
  gsmooth->SetPoint(ipt, voltage_mean, signal_mean);
  
  gsmooth->SetMarkerStyle(20);
  //gsmooth->Print("all");
  
  return gsmooth;
  
}


TGraphErrors* GetSmoothGraph(TGraphErrors* g)
{

  TGraphErrors* gsmooth = new TGraphErrors();
  double voltage=0;
  double previous_voltage=0;
  double next_voltage=0;
  double signal=0;
  double previous_signal=0;
  double next_signal=0;
  int ipt=0;
  double volt=0;
  double signal_mean=0;
  double voltage_mean=0;
  
  for(unsigned int ip=0; ip<g->GetN(); ip++)
  {
   previous_signal = signal;
   previous_voltage = voltage;
   g->GetPoint(ip, voltage, signal);
   g->GetPoint(ip+1, next_voltage, next_signal);
   signal_mean = ( previous_signal + signal + next_signal )/3.;
   voltage_mean = ( previous_voltage + voltage + next_voltage )/3.;
   if(!ip){
    signal_mean = ( signal + next_signal )/2.;
    voltage_mean = ( voltage + next_voltage )/2.;
   }
   if(ip==g->GetN()-1){
    signal_mean = ( previous_signal + signal )/2.;
    voltage_mean = ( previous_voltage + voltage )/2.;	
   }
   gsmooth->SetPoint(ipt, voltage_mean, signal_mean);
   ipt++;
  }
  gsmooth->SetMarkerStyle(20);
 
  return gsmooth;
}

TGraphErrors* GetSmoothGraph_NPtsMean(TGraphErrors* g, int N=4)
{

  TGraphErrors* gsmooth = new TGraphErrors();
  double voltage=0;
  double signal=0;
  int ipt=0;
  double signal_mean=0;
  double voltage_mean=0;
  if(N<2 || N>g->GetN()) return gsmooth;
  
  for(unsigned int ip=0; ip<g->GetN()+1-N; ip++)
  {
    voltage_mean = 0;
    signal_mean = 0;
	for(int i=0; i<N; i++)
	{
      g->GetPoint(ip+i, voltage, signal);
	  voltage_mean+=voltage;
	  signal_mean+=signal;
	}
	voltage_mean/=N;
    signal_mean/=N;
    gsmooth->SetPoint(ipt, voltage_mean, signal_mean);
    ipt++;
  }
  gsmooth->SetMarkerStyle(20);
 
  return gsmooth;
}

TGraphErrors* Interpolate(TGraphErrors* g)
{
  TGraphErrors* gsmooth = new TGraphErrors();
  
  // WRITTEN TO DOUBLE N POINTS
  
  const int N=4; // FUNCTION WRITTEN FOR CUBIC INTERPOLATION -> 4 POINTS NEEDED
  double* voltage = new double[N];
  double* signal = new double[N];
  int ipt=0;
  
  for(unsigned int ip=0; ip<g->GetN(); ip++)
  {
    // Set existing point
    g->GetPoint(ip, voltage[0], signal[0]);
    gsmooth->SetPoint(ipt, voltage[0], signal[0]);
    ipt++;
    
    // Add interpolated point
    if(ip>=(N/2-1) && ip<=(g->GetN()-1-N/2))
    {
      for(int i=0; i<N; i++)
      {
        g->GetPoint(ip+1-N/2+i, voltage[i], signal[i]);
      }
      double mu=0.5; // add 1 point in the middle of the interval
      double a0,a1,a2,a3; 
      a0 = -0.5*signal[0] + 1.5*signal[1] - 1.5*signal[2] + 0.5*signal[3];
      a1 = signal[0] - 2.5*signal[1] + 2*signal[2] - 0.5*signal[3];
      a2 = -0.5*signal[0] + 0.5*signal[2];
      a3 = signal[1];
      gsmooth->SetPoint(ipt, voltage[1]*(1-mu)+voltage[2]*mu, a0*mu*mu*mu+a1*mu*mu+a2*mu+a3);
      ipt++;
    }
  }
  
  gsmooth->SetMarkerStyle(20);
  //gsmooth->Print("all");
  
  return gsmooth;  
}

int CorrectGraphForLeakageCurrent(TGraph* g, ULong64_t detid, string run)
{
  // Get correction function
  string filename="../../../LeakageCurrentBiasScan/Analysis/LeakCurCorr_"+run+".root";
  TFile* f = new TFile(filename.c_str(), "read");
  if(!f) {cout<<"No leakage current correction file : "<<filename<<endl; return 0;}
  TF1* func = (TF1*) f->Get(Form("fit_%llu", detid));   
  if(!func) {cout<<"No leakage current correction for detid : "<<detid<<endl; return 0;}

  // Correct points
  double x,y;
  for(int i=0; i<g->GetN(); i++)
  {
    g->GetPoint(i, x, y);
    g->SetPoint(i, x-func->Eval(x), y);
  }
  
  f->Close();
  
  return 1;
}
