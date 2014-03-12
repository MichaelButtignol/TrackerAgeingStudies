{

  gROOT->Reset();

  gROOT->ProcessLine(".L ../../interface/TreeEvent.h+");  
  gROOT->ProcessLine(".L Code/VoltageStepsTreeMaker.C+");

  if(gROOT->GetClass("TreeEvent") == 0) return;
  if(gROOT->GetClass("VoltageStepsTreeMaker") == 0) return;

  TChain c("ttree");
  c.Add("root://eoscms//eos/cms/store/group/comm_tracker/Strip/RadMonitoring/SignalBiasScan/ClustersTrees/DecoSmallHVScan_20130213_run211797_v1_1/clustersTree*.root/demo/ttree");

  int subdet=0; // 0 all, 1 TIB, 2 TOB, 3 TID, 4 TEC
  bool usetimestamp=false; // true -> uses timestamps for steps definition instead of event numbers
  VoltageStepsTreeMaker *t = new VoltageStepsTreeMaker(&c, "Steps/Steps_20130213_run211797.txt", subdet, usetimestamp);
  t->Loop();

}

