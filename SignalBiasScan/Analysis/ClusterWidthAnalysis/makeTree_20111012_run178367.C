{

  gROOT->Reset();

  gROOT->ProcessLine(".L ../../interface/TreeEvent.h+");  
  gROOT->ProcessLine(".L ../VoltageSteps/Code/VoltageStepsMaker.C+");
  gROOT->ProcessLine(".L Code/ClusterWidthAnalysisTreeMaker.C+");

  if(gROOT->GetClass("TreeEvent") == 0) return;
  if(gROOT->GetClass("VoltageStepsMaker") == 0) return;
  if(gROOT->GetClass("ClusterWidthAnalysisTreeMaker") == 0) return;

  TChain c("ttree");
  c.Add("root://eoscms//eos/cms/store/group/comm_tracker/Strip/RadMonitoring/SignalBiasScan/ClustersTrees/DecoHVScan_20111012_run178367_v1_2/clustersTree*.root/demo/ttree");

  int subdet=1; // 0 all, 1 TIB, 2 TOB, 3 TID, 4 TEC
  bool usetimestamp=true; // true -> uses timestamps for steps definition instead of event numbers
  ClusterWidthAnalysisTreeMaker *t = new ClusterWidthAnalysisTreeMaker(&c, "../VoltageSteps/Steps/Steps_20111012_run178367.txt", subdet, usetimestamp);
  t->Loop();

}

