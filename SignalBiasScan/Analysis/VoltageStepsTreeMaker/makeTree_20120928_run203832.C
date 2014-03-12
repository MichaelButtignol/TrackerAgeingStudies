{

  gROOT->Reset();

  gROOT->ProcessLine(".L ../../interface/TreeEvent.h+");  
  gROOT->ProcessLine(".L Code/VoltageStepsTreeMaker.C+");

  if(gROOT->GetClass("TreeEvent") == 0) return;
  if(gROOT->GetClass("VoltageStepsTreeMaker") == 0) return;

  TChain c("ttree");
  c.Add("root://eoscms//eos/cms/store/group/comm_tracker/Strip/RadMonitoring/SignalBiasScan/ClustersTrees/DecoSmallHVScan_20120928_run203832_v1.2/clustersTrees*.root/demo/ttree");

  int subdet=1; // 0 all, 1 TIB, 2 TOB, 3 TID, 4 TEC
  bool usetimestamp=true; // true -> uses timestamps for steps definition instead of event numbers
  VoltageStepsTreeMaker *t = new VoltageStepsTreeMaker(&c, "Steps/Steps_20120928_run203832_merged_steps.txt", subdet, usetimestamp);
  t->Loop();

}

