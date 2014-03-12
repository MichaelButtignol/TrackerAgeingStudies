{

  gROOT->Reset();

  gROOT->ProcessLine(".L ../../interface/TreeEvent.h+");  
  gROOT->ProcessLine(".L Code/VoltageStepsTreeMaker.C+");

  if(gROOT->GetClass("TreeEvent") == 0) return;
  if(gROOT->GetClass("VoltageStepsTreeMaker") == 0) return;

  TChain c("ttree");
  //c.Add("root://eoscms//eos/cms/store/group/comm_tracker/Strip/RadMonitoring/SignalBiasScan/ClustersTrees/DecoSmallHVScan_20121130_run208339_alldetids_v1_1/clustersTree*.root/demo/ttree");
  c.Add("root://eoscms//eos/cms/store/group/comm_tracker/Strip/RadMonitoring/SignalBiasScan/ClustersTrees/DecoSmallHVScan_20121130_run208339_v1_1/clustersTree*.root/demo/ttree");

  int subdet=2; // 0 all, 1 TIB, 2 TOB, 3 TID, 4 TEC
  bool usetimestamp=false; // true -> uses timestamps for steps definition instead of event numbers
  VoltageStepsTreeMaker *t = new VoltageStepsTreeMaker(&c, "Steps/Steps_20121130_run208339_merged_steps.txt", subdet, usetimestamp);
  t->Loop();

}

