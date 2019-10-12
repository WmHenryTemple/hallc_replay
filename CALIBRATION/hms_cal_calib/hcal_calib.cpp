#include "TCanvas.h"
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
 theShowerCalib.FillHEcalNoCor();       // Fill histograms
 theShowerCalib.FillCutBranch();       // Fill diagnostic histos
 theShowerCalib.FillHitsGains();       // Plot hits and gains

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

 theShowerCalib.hEcalNoCor->SetLineColor(kGreen);
 theShowerCalib.hEcalNoCor->Draw("same");

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
 theShowerCalib.pr1->Draw("TEXT COLZ");
 Canvas2->cd(4);
 theShowerCalib.ta2->Draw("TEXT COLZ");
 Canvas2->cd(5);
 theShowerCalib.ta3->Draw("TEXT COLZ");
 Canvas2->cd(6);
 theShowerCalib.ta4->Draw("TEXT COLZ");

 Canvas3->Divide(3,2);
 Canvas3->cd(1);
 theShowerCalib.hCer->Draw();
 Canvas3->cd(2);
 theShowerCalib.hP->Draw();
 Canvas3->cd(3);
 theShowerCalib.hDelta->Draw();
 Canvas3->cd(4);
 theShowerCalib.hBeta->Draw();
 Canvas3->cd(5);

Canvas4->cd();
 theShowerCalib.pmtList->Draw();

 }
  // Calculate the analysis rate

 

   TFile *oFile=new TFile("hcal_out.root","RECREATE");
   theShowerCalib.hEcal->Write();
   oFile->Close();

 t = clock() - t;
 printf ("The analysis took %.1f seconds \n", ((float) t) / CLOCKS_PER_SEC);


}
