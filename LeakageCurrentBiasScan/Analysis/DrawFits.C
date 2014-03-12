{
  char* run="TIB_L1_20120506_run193541";
  
  TFile* f = new TFile(Form("Results/LeakCurCorr_%s.root", run), "read");
  if(!f) cout<<"File not found."<<endl;
  
  TH1F* h = new TH1F("h", "", 35, 0, 350);
  h->Draw();
  h->SetMaximum(6);
  gStyle->SetOptStat(0);
  
  int PS=0;
  int prev_PS=0;
  int icolor=1;
  
  gDirectory->cd();
  TIter nextkey( f->GetListOfKeys() );
  TKey *key;
  while( (key = (TKey*)nextkey()))
  {
    TObject *obj = key->ReadObj();
    if( obj->IsA()->InheritsFrom("TF1") )
    {
      TF1 *fit = (TF1*)obj;
      string name(fit->GetName());
      cout<<name<<endl;
      // get prefix indicating power group
      name.replace(0,4,"");
      name.replace(7,2,"");
      stringstream ss(name);
      ss >> PS;
      
      //if(PS!=prev_PS && prev_PS) icolor++;
      //if(icolor==3) icolor++;
      
      // From 3 different PS
      if(PS==3691213) icolor=1;
      if(PS==3691216) icolor=2;
      if(PS==3691258) icolor=4;
      
      fit->SetLineColor(icolor);
      fit->Draw("same");
      prev_PS=PS;
    } 
  }
      
  
}
