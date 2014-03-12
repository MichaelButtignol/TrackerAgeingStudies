{

  gROOT->Reset();

  gROOT->ProcessLine(".L ../../interface/TreeEvent.h+");  
  gROOT->ProcessLine(".L Code/VoltageStepsTreeMaker.C+");

  if(gROOT->GetClass("TreeEvent") == 0) return;
  if(gROOT->GetClass("VoltageStepsTreeMaker") == 0) return;

  TChain c("ttree");
  c.Add("file:/afs/cern.ch/user/j/jlagram/work/ClustersTrees/clustersTree*.root/demo/ttree");

  int subdet=1; // 0 all, 1 TIB, 2 TOB, 3 TID, 4 TEC
  bool usetimestamp=true; // true -> uses timestamps for steps definition instead of event numbers
  VoltageStepsTreeMaker *t = new VoltageStepsTreeMaker(&c, "Steps/Steps_20120812_run200786.txt", subdet, usetimestamp);
  t->Loop();

}

