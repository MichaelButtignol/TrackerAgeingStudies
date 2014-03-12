{

 char* subdet="TEC";
 TFile *f = TFile::Open(Form("%s_output_DecoMode.root", subdet));
 TTree* t = (TTree*) f->Get("T");

  ULong64_t detid;
  double volt;
  int id;
  double evolt;
  double mpv;
  double empv;
  ULong64_t tempdetid=0;
  double chi2overndf;
  double pedfraction;

  t->SetBranchAddress("DetID",&detid);
  t->SetBranchAddress("Voltage",&volt);
  t->SetBranchAddress("Index",&id);
  t->SetBranchAddress("errVoltage",&evolt);
  t->SetBranchAddress("MPV",&mpv);
  t->SetBranchAddress("errMPV",&empv);
  t->SetBranchAddress("Chi2OverNdf",&chi2overndf);
  t->SetBranchAddress("PedFraction",&pedfraction);

 ofstream fout(Form("PedFraction_%s.txt", subdet)); 

 double volt1, volt2;

 TH1F* h = new TH1F("vdiff","",50,0,0.05);
 TH2F* h2 = new TH2F("h2","",50,0,0.1,50,0,0.05);
 unsigned int nent = t->GetEntriesFast();
 for(unsigned int i=0; i<nent; i++)
 {
  t->GetEntry(i);

  if(detid==tempdetid)
  {
    if(volt==5) volt1=pedfraction;
    if(volt==300) volt2=pedfraction;
  }
  else
  {
   if(volt2>0.04)h->Fill(-volt2+volt1);
   h2->Fill(volt1, -volt2+volt1); 
   tempdetid=detid;
  }
  if(volt==300) fout<<detid<<" "<<pedfraction<<endl; 
 }
 fout.close();

 gStyle->SetPalette(1);
 h2->Draw("colz");
}
