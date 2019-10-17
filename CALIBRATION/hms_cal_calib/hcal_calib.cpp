#include "TCanvas.h"
#include "TLine.h"
#include <TStyle.h>
#include <TH1.h>
#include <TEllipse.h>
#include <TF1.h>
#include "THcShowerCalib.h"

//
// A steering Root script for the HMS calorimeter calibration.
//

void hcal_calib(string Prefix, int nstop=-1, int nstart=0) {

  bool DRAW = 1;

  // Initialize the analysis clock
  clock_t t = clock();
  
  cout << "Calibrating file " << Prefix << ".root, events "
       << nstart << " -- " << nstop << endl;

  THcShowerCalib theShowerCalib(Prefix, nstart, nstop);

 theShowerCalib.ReadThresholds();  // Read in threshold param-s and intial gains
 theShowerCalib.Init();            // Initialize constants and variables
 theShowerCalib.CalcThresholds();  // Thresholds on the uncalibrated Edep/P
 theShowerCalib.ComposeVMs();      // Compute vectors amd matrices for calib.
 theShowerCalib.SolveAlphas();     // Solve for the calibration constants
 theShowerCalib.SaveAlphas();      // Save the constants
 // theShowerCalib.SaveRawData();  // Save raw data into file for debug purposes
 theShowerCalib.FillHEcal();       // Fill histograms
 // theShowerCalib.FillHEcalNoCor();  // Removes Y correction (pulseInt * gain * yCor) = E
 theShowerCalib.FillCutBranch();       // Fill diagnostic histos
 theShowerCalib.FillHitsGains();       // Plot hits and gains

 cout << "Test: fDeltamax:"<<theShowerCalib.GetDeltaMax()<<endl;
 // Plot histograms

 TCanvas* Canvas =
   new TCanvas("Canvas", "HMS Shower Counter calibration", 1000, 667);

 Canvas->Divide(2,2);

 Canvas->cd(1);

 // Normalized uncalibrated energy deposition.

 theShowerCalib.hEunc->DrawCopy();
  
 theShowerCalib.hEuncSel->SetFillColor(kGreen);
 theShowerCalib.hEuncSel->DrawCopy("same");

 // E(Tail) vs E(Preshower) plot.

 Canvas->cd(2);
 theShowerCalib.hETAvsEPR->Draw("colz");

 // Normalized energy deposition after calibration.

 Canvas->cd(3);
 gStyle->SetOptFit();

 //wph Fit over range centered by max bin
 Double_t maxBin= theShowerCalib.hEcal->GetMaximumBin();
 Double_t maxValue= theShowerCalib.hEcal->GetBinCenter(maxBin);
 theShowerCalib.hEcal->Fit("gaus","O","",maxValue-.1, maxValue+.2);
 // cout << "Max Value:" << maxValue <<endl;
 // cout << "At Bin #:" << maxBin <<endl;

 // theShowerCalib.hEcal->Fit("gaus","O","",0.5,1.5);

 TF1 *fit = theShowerCalib.hEcal->GetFunction("gaus");
 Double_t gmean  = fit->GetParameter(1);
 Double_t gsigma = fit->GetParameter(2);
 double gLoThr = gmean - 1.*gsigma;
 double gHiThr = gmean + 2.*gsigma;
 cout << "gLoThr=" << gLoThr << "  gHiThr=" << gHiThr << endl;
 theShowerCalib.hEcal->Fit("gaus","","",gLoThr,gHiThr);

 theShowerCalib.hEcal->GetFunction("gaus")->SetLineColor(2);
 theShowerCalib.hEcal->GetFunction("gaus")->SetLineWidth(1);
 theShowerCalib.hEcal->GetFunction("gaus")->SetLineStyle(1);

 // theShowerCalib.hEcalNoCor->SetLineColor(kGreen);
 // theShowerCalib.hEcalNoCor->Draw("same");

 // HMS delta(P) versus the calibrated energy deposition.

 Canvas->cd(4);
 theShowerCalib.hDPvsEcal->Draw("colz");

 // Save canvas in a pdf format.
 Canvas->Print(Form("%s_%d_%d.pdf",Prefix.c_str(),nstart,nstop));

 if(DRAW==1){
 TCanvas* Canvas2 =
   new TCanvas("Canvas2", "HMS Shower Counter calibration", 1000, 667);

 TCanvas* Canvas3 =
   new TCanvas("Canvas3", "HMS Shower Counter calibration", 1000, 667);

 TCanvas* Canvas4 =
   new TCanvas("Canvas4", "HMS Shower Counter calibration", 1000, 667);

 Canvas2->Divide(3,2);
 Canvas2->cd(1);
 theShowerCalib.hCaloPos->Draw("COLZ");

 Canvas2->cd(2);

 TEllipse *el1=new TEllipse( 0, 2.8, 46.507, 46.507);
 el1->SetLineColor(kRed);
 theShowerCalib.hExitPos->Draw("COLZ");
 el1->Draw("same");
 theShowerCalib.hExitPos->Draw("COLZ SAME");

 Canvas2->cd(3);
 theShowerCalib.pr1->SetStats(0);
 theShowerCalib.pr1->Draw("TEXT COLZ");
 theShowerCalib.pr1a->Draw("TEXT SAME");
 Canvas2->cd(4);
 theShowerCalib.ta2->Draw("TEXT COLZ");
 theShowerCalib.ta2->SetStats(0);
 theShowerCalib.ta3a->Draw("TEXT SAME");
 Canvas2->cd(5);
 theShowerCalib.ta3->Draw("TEXT COLZ");
 theShowerCalib.ta3->SetStats(0);
 theShowerCalib.ta3a->Draw("TEXT SAME");
 Canvas2->cd(6);
 theShowerCalib.ta4->SetStats(0);
 theShowerCalib.ta4->Draw("TEXT COLZ");
 theShowerCalib.ta4a->Draw("TEXT SAME");

 Double_t maxHit=1.3*theShowerCalib.pr1->GetMaximum();

 theShowerCalib.pr1->GetZaxis()->SetRangeUser(0,maxHit);
 theShowerCalib.ta2->GetZaxis()->SetRangeUser(0,maxHit);
 theShowerCalib.ta3->GetZaxis()->SetRangeUser(0,maxHit);
 theShowerCalib.ta4->GetZaxis()->SetRangeUser(0,maxHit);

 theShowerCalib.pr1a->SetMarkerColor(kRed);
 theShowerCalib.ta2a->SetMarkerColor(kRed);
 theShowerCalib.ta3a->SetMarkerColor(kRed);
 theShowerCalib.ta4a->SetMarkerColor(kRed);

 Canvas3->Divide(3,2);
 // draw 1D npeSum 
 Canvas3->cd(1);
 theShowerCalib.hCer->Draw();
 gPad->SetLogy();
 Double_t top=theShowerCalib.hCer->GetMaximum();
 Double_t locx=theShowerCalib.GetCerMin();
 TLine *cerl=new TLine(locx,0,locx,top);
 cerl->SetLineColor(kMagenta);
 cerl->Draw("same");

 Canvas3->cd(2);
 theShowerCalib.hP->Draw();


 // draw 1D Delta 
 Canvas3->cd(3);
 theShowerCalib.hDelta->Draw();
 top=theShowerCalib.hDelta->GetMaximum();
 locx=theShowerCalib.GetDeltaMin();
 TLine *dplmin=new TLine(locx,0,locx,top);
 dplmin->SetLineColor(kMagenta);
 dplmin->Draw("same");
 locx=theShowerCalib.GetDeltaMax();
 TLine *dplmax=new TLine(locx,0,locx,top);
 dplmax->SetLineColor(kMagenta);
 dplmax->Draw("same");

 // draw 1D Beta 
 Canvas3->cd(4);
 theShowerCalib.hBeta->Draw();
 gPad->SetLogy();
 top=theShowerCalib.hBeta->GetMaximum();
 locx=theShowerCalib.GetBetaMin();
 TLine *betalmin=new TLine(locx,0,locx,top);
 betalmin->SetLineColor(kMagenta);
 betalmin->Draw("same");
 locx=theShowerCalib.GetBetaMax();
 TLine *betalmax=new TLine(locx,0,locx,top);
 betalmax->SetLineColor(kMagenta);
 betalmax->Draw("same");
 Canvas3->cd(5);
 // theShowerCalib.hClusTrk->Draw("BOX");
 gPad->SetLogy();
 theShowerCalib.hNclust->Draw();
 Canvas3->cd(6);
 gPad->SetLogy();
 theShowerCalib.hNtrack->Draw();

 Canvas4->Divide(2,2);
 Canvas4->cd(1);
 theShowerCalib.pmtList->Draw("COLZ");
 Canvas4->cd(2);
 
Float_t adcmax=0;
Float_t temp=0;

 for(UInt_t i=1; i<78; i++){
   temp=theShowerCalib.hAdc[i]->GetMaximum();
   if (temp>adcmax)adcmax=temp;
}
 theShowerCalib.hAdc[0]->Draw();
 theShowerCalib.hAdc[0]->GetYaxis()->SetRangeUser(0,adcmax);

 for(UInt_t i=1; i<78; i++){theShowerCalib.hAdc[i]->Draw("same");}

  // Calculate the analysis rate
 //   TFile *oFile=new TFile("hcal_out.root","RECREATE");
 //  theShowerCalib.hEcal->Write();
 //  oFile->Close();
 Canvas4->cd(3);
 theShowerCalib.yCalVsEp->Draw("COLZ");
 Canvas4->cd(4);
 theShowerCalib.xCalVsEp->Draw("COLZ");

 Canvas2->Print(Form("hits_%s_%d_%d.pdf",Prefix.c_str(),nstart,nstop));
 Canvas3->Print(Form("cuts_%s_%d_%d.pdf",Prefix.c_str(),nstart,nstop));

 }


 t = clock() - t;
 printf ("The analysis took %.1f seconds \n", ((float) t) / CLOCKS_PER_SEC);


}
