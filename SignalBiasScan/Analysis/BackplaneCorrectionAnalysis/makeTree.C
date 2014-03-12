




// Still in processing----------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//---------------------------------------------------------------------------------




void main()
{
  	gROOT->Reset();

  	gROOT->ProcessLine(".L ../../interface/TreeEvent.h+");  
	gROOT->ProcessLine(".L ../VoltageSteps/Code/VoltageStepsMaker.C+");
	gROOT->ProcessLine(".L Code/BackplaneCorrectionAnalysisTreeMaker.C+");

	if(gROOT->GetClass("TreeEvent") == 0) return;
	if(gROOT->GetClass("VoltageStepsMaker") == 0) return;
	if(gROOT->GetClass("BackplaneCorrectionAnalysisTreeMaker") == 0) return;	

	TChain c("ttree");
	c.Add("root://eoscms//eos/cms/store/group/comm_tracker/Strip/RadMonitoring/SignalBiasScan/ClustersTrees/DecoHVScan_20111027_run180076_v1_2/clustersTree*.root/demo/ttree");
	int subdet=1; // 0 all, 1 TIB, 2 TOB, 3 TID, 4 TEC
  bool usetimestamp=true; // true -> uses timestamps for steps definition instead of event numbers
  BackplaneCorrectionAnalysisTreeMaker *t = new BackplaneCorrectionAnalysisTreeMaker(&c, "../VoltageSteps/Steps/Steps_20111027_run180076.txt", subdet, usetimestamp);
  t->Loop();
  
  
	TChain c("ttree");
  c.Add("root://eoscms//eos/cms/store/group/comm_tracker/Strip/RadMonitoring/SignalBiasScan/ClustersTrees/DecoHVScan_20100728_run141865_v1_2/clustersTree*.root/demo/ttree");

  int subdet=1; // 0 all, 1 TIB, 2 TOB, 3 TID, 4 TEC
  bool usetimestamp=true; // true -> uses timestamps for steps definition instead of event numbers
  BackplaneCorrectionAnalysisTreeMaker *t = new BackplaneCorrectionAnalysisTreeMaker(&c, "../VoltageSteps/Steps/Steps_20100728_run141865.txt", subdet, usetimestamp);
  t->Loop();

  TChain c("ttree");
  c.Add("root://eoscms//eos/cms/store/group/comm_tracker/Strip/RadMonitoring/SignalBiasScan/ClustersTrees/DecoHVScan_20110315_run160497_v1_2/clustersTree*.root/demo/ttree");

  int subdet=1; // 0 all, 1 TIB, 2 TOB, 3 TID, 4 TEC
  bool usetimestamp=false; // true -> uses timestamps for steps definition instead of event numbers
  BackplaneCorrectionAnalysisTreeMaker *t = new BackplaneCorrectionAnalysisTreeMaker(&c, "../VoltageSteps/Steps/Steps_20110315_run160497.txt", subdet, usetimestamp);
  t->Loop();
  
    TChain c("ttree");
  c.Add("root://eoscms//eos/cms/store/group/comm_tracker/Strip/RadMonitoring/SignalBiasScan/ClustersTrees/DecoHVScan_20111012_run178367_v1_2/clustersTree*.root/demo/ttree");

  int subdet=1; // 0 all, 1 TIB, 2 TOB, 3 TID, 4 TEC
  bool usetimestamp=true; // true -> uses timestamps for steps definition instead of event numbers
  BackplaneCorrectionAnalysisTreeMaker *t = new BackplaneCorrectionAnalysisTreeMaker(&c, "../VoltageSteps/Steps/Steps_20111012_run178367.txt", subdet, usetimestamp);
  t->Loop();
  
    TChain c("ttree");
  c.Add("root://eoscms//eos/cms/store/group/comm_tracker/Strip/RadMonitoring/SignalBiasScan/ClustersTrees/DecoHVScan_20120405_run190459_v1_2/clustersTree_*.root/demo/ttree");
  int subdet=0; // 0 all, 1 TIB, 2 TOB, 3 TID, 4 TEC
  bool usetimestamp=true; // true -> uses timestamps for steps definition instead of event numbers
  BackplaneCorrectionAnalysisTreeMaker *t = new BackplaneCorrectionAnalysisTreeMaker(&c, "../VoltageSteps/Steps/Steps_20120405_run190459.txt", subdet, usetimestamp);
  t->Loop();
  
    TChain c("ttree");
  c.Add("root://eoscms//eos/cms/store/group/comm_tracker/Strip/RadMonitoring/SignalBiasScan/ClustersTrees/DecoHVScan_20120812_200786_v1_2/clustersTree_*.root/demo/ttree");
  int subdet=0; // 0 all, 1 TIB, 2 TOB, 3 TID, 4 TEC
  bool usetimestamp=true; // true -> uses timestamps for steps definition instead of event numbers
  BackplaneCorrectionAnalysisTreeMaker *t = new BackplaneCorrectionAnalysisTreeMaker(&c, "../VoltageSteps/Steps/Steps_20120812_run200786.txt", subdet, usetimestamp);
  t->Loop();
  
    TChain c("ttree");
  c.Add("root://eoscms//eos/cms/store/group/comm_tracker/Strip/RadMonitoring/SignalBiasScan/ClustersTrees/DecoSmallHVScan_20120928_run203832_v1_2/clustersTree*.root/demo/ttree");
  int subdet=1; // 0 all, 1 TIB, 2 TOB, 3 TID, 4 TEC
  bool usetimestamp=true; // true -> uses timestamps for steps definition instead of event numbers
  BackplaneCorrectionAnalysisTreeMaker *t = new BackplaneCorrectionAnalysisTreeMaker(&c, "../VoltageSteps/Steps/Steps_20120928_run203832_merged_steps.txt", subdet, usetimestamp);
  t->Loop();
	  TChain c("ttree");
  c.Add("root://eoscms//eos/cms/store/group/comm_tracker/Strip/RadMonitoring/SignalBiasScan/ClustersTrees/DecoSmallHVScan_20120928_run203832_v1_2/clustersTree*.root/demo/ttree");
  int subdet=1; // 0 all, 1 TIB, 2 TOB, 3 TID, 4 TEC
  bool usetimestamp=true; // true -> uses timestamps for steps definition instead of event numbers
  BackplaneCorrectionAnalysisTreeMaker *t = new BackplaneCorrectionAnalysisTreeMaker(&c, "../VoltageSteps/Steps/Steps_20120928_run203832_merged_steps.txt", subdet, usetimestamp);
  t->Loop();
  
   TChain c("ttree");
  c.Add("root://eoscms//eos/cms/store/group/comm_tracker/Strip/RadMonitoring/SignalBiasScan/ClustersTrees/DecoSmallHVScan_20121130_run208339_v1_2/clustersTree_*.root/demo/ttree");
  int subdet=0; // 0 all, 1 TIB, 2 TOB, 3 TID, 4 TEC
  bool usetimestamp=true; // true -> uses timestamps for steps definition instead of event numbers
  BackplaneCorrectionAnalysisTreeMaker *t = new BackplaneCorrectionAnalysisTreeMaker(&c, "../VoltageSteps/Steps/Steps_20121130_run208339_merged_steps.txt", subdet, usetimestamp);
  t->Loop();
  
  
    TChain c("ttree");
  c.Add("root://eoscms//eos/cms/store/group/comm_tracker/Strip/RadMonitoring/SignalBiasScan/ClustersTrees/DecoSmallHVScan_20130213_run211797_v1_1/clustersTree_*.root/demo/ttree");
  int subdet=0; // 0 all, 1 TIB, 2 TOB, 3 TID, 4 TEC
  bool usetimestamp=false; // true -> uses timestamps for steps definition instead of event numbers
  BackplaneCorrectionAnalysisTreeMaker *t = new BackplaneCorrectionAnalysisTreeMaker(&c, "../VoltageSteps/Steps/Steps_20130213_run211797.txt", subdet, usetimestamp);
  t->Loop();
}
