{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Tue Feb 18 10:44:59 2014) by ROOT version5.32/00
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",130,75,1146,823);
   Canvas_1->Range(-78.42278,-3744.635,421.9058,25463.51);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetLeftMargin(0.1567426);
   Canvas_1->SetRightMargin(0.04378284);
   Canvas_1->SetTopMargin(0.07081807);
   Canvas_1->SetBottomMargin(0.1282051);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   
   TH1F *hNevtPerStep_TIB = new TH1F("hNevtPerStep_TIB","hNevtPerStep",80,0,400);
   hNevtPerStep_TIB->SetBinContent(7,10331);
   hNevtPerStep_TIB->SetBinContent(10,10570);
   hNevtPerStep_TIB->SetBinContent(13,10538);
   hNevtPerStep_TIB->SetBinContent(16,10447);
   hNevtPerStep_TIB->SetBinContent(19,10776);
   hNevtPerStep_TIB->SetBinContent(22,10970);
   hNevtPerStep_TIB->SetBinContent(25,11153);
   hNevtPerStep_TIB->SetBinContent(28,11048);
   hNevtPerStep_TIB->SetBinContent(31,11357);
   hNevtPerStep_TIB->SetBinContent(34,11543);
   hNevtPerStep_TIB->SetBinContent(37,11812);
   hNevtPerStep_TIB->SetBinContent(40,12005);
   hNevtPerStep_TIB->SetBinContent(43,22281);
   hNevtPerStep_TIB->SetBinContent(46,17263);
   hNevtPerStep_TIB->SetBinContent(49,12965);
   hNevtPerStep_TIB->SetBinContent(52,13068);
   hNevtPerStep_TIB->SetBinContent(55,12875);
   hNevtPerStep_TIB->SetBinContent(58,13017);
   hNevtPerStep_TIB->SetBinContent(61,13626);
   hNevtPerStep_TIB->SetBinContent(66,13616);
   hNevtPerStep_TIB->SetBinContent(71,16008);
   hNevtPerStep_TIB->SetEntries(267269);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.755,0.98,0.995,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   TText *text = ptstats->AddText("hNevtPerStep_TIB");
   text->SetTextSize(0.0368);
   text = ptstats->AddText("Entries = 267269 ");
   text = ptstats->AddText("Mean  =  192.8");
   text = ptstats->AddText("RMS   =  91.61");
   text = ptstats->AddText("Underflow =      0");
   text = ptstats->AddText("Overflow  =      0");
   ptstats->SetOptStat(111111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   hNevtPerStep_TIB->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(hNevtPerStep_TIB);
   hNevtPerStep_TIB->GetXaxis()->SetTitle("V_{bias} [V]");
   hNevtPerStep_TIB->GetXaxis()->SetTitleSize(0.05);
   hNevtPerStep_TIB->GetXaxis()->SetTitleOffset(1.2);
   hNevtPerStep_TIB->GetYaxis()->SetTitle("N_{events}");
   hNevtPerStep_TIB->GetYaxis()->SetTitleSize(0.05);
   hNevtPerStep_TIB->GetYaxis()->SetTitleOffset(1.5);
   hNevtPerStep_TIB->Draw("");
   
   TPaveText *pt = new TPaveText(0.07880911,0.9426129,0.2425569,0.995116,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(1);
   pt->SetFillColor(0);
   text = pt->AddText("hNevtPerStep");
   pt->Draw();
   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
}
