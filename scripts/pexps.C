#include "TF1.h"
#include "TH1F.h"
#include "TFile.h"
#include "TText.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TFractionFitter.h"
#include "TFitResult.h"

#include <iostream>
using std::cout;
using std::endl;


TH1F* t1;
TH1F* t2;
TH1F* tbg;

Double_t fitf(Double_t *x,Double_t *par) {
  int ibin = t1->GetXaxis()->FindBin(x[0]);
  Double_t p1 = t1->GetBinContent(ibin);
  Double_t p2 = t2->GetBinContent(ibin);
  Double_t b  = tbg->GetBinContent(ibin);
  Double_t pred = par[0]*( par[1]*p1 + (1.-par[1])*p2 ) + b;
  return pred;
}


void pexps( int ntoyExps=100 ) {

  const int nfile=8;
  TFile* fin[nfile];

  fin[0] = new TFile("figs/ILD_l5_o1_v02_eL80pR30/fithistos.root","read");
  fin[1] = new TFile("figs/ILD_l5_o1_v02_eR80pL30/fithistos.root","read");
  fin[2] = new TFile("figs/ILD_l5_o1_v02_eL80pL30/fithistos.root","read");
  fin[3] = new TFile("figs/ILD_l5_o1_v02_eR80pR30/fithistos.root","read");

  fin[4] = new TFile("figs/ILD_s5_o1_v02_eL80pR30/fithistos.root","read");
  fin[5] = new TFile("figs/ILD_s5_o1_v02_eR80pL30/fithistos.root","read");
  fin[6] = new TFile("figs/ILD_s5_o1_v02_eL80pL30/fithistos.root","read");
  fin[7] = new TFile("figs/ILD_s5_o1_v02_eR80pR30/fithistos.root","read");

  TString lab[nfile]={"L_lr", "L_rl", "L_ll", "L_rr", "S_lr", "S_rl", "S_ll", "S_rr"};

  TText tt;

  TH1F* hFrac[nfile][2];
  TH1F* hFracErr[nfile][2];
  TH1F* hStatus[nfile][2];
  TH1F* hFitProb[nfile][2];

  TCanvas* cc = new TCanvas("cc","cc",400,400);
  cc->Print("fit.pdf[");

  TRandom3 trand3;

  float scale = 1; // SET TO 1 FOR FINAL RESULTS!!!!

  for (int iset=0; iset<nfile; iset++) {

    int icol = iset<4 ? 4 : 2;
    int ilty = iset<4 ? 1 : 2;

    for (int idec=0; idec<2; idec++) {

      TString dec = idec==0 ? "pi" : "rho";

      hFrac   [iset][idec] = new TH1F(lab[iset]+"_fitFrac_"+dec,    lab[iset]+"_fitFrac_"+dec,    100, 0, 1);
      hFracErr[iset][idec] = new TH1F(lab[iset]+"_fitFracErr_"+dec, lab[iset]+"_fitFracErr_"+dec, 5000, 0., 0.02);
      hStatus [iset][idec] = new TH1F(lab[iset]+"_fitStatus_"+dec,  lab[iset]+"_fitStatus_"+dec,  10, 0, 10);
      hFitProb[iset][idec] = new TH1F(lab[iset]+"_fitProb_"+dec,    lab[iset]+"_fitProb_"+dec,    100, 0, 1);


      TString hname="SEL_rec_" + dec + "_pol";

      TH1F* sigNeg = (TH1F*) fin[iset]->Get("group0_"+hname+"_MCneg");
      sigNeg->Add( (TH1F*) fin[iset]->Get("group1_"+hname+"_MCneg") );

      TH1F* sigPos = (TH1F*) fin[iset]->Get("group0_"+hname+"_MCpos");
      sigNeg->Add( (TH1F*) fin[iset]->Get("group1_"+hname+"_MCpos") );

      TH1F* bg = (TH1F*) fin[iset]->Get("group0_"+hname+"_MCoth");
      bg->Add( (TH1F*) fin[iset]->Get("group1_"+hname+"_MCoth") );
      bg->Add( (TH1F*) fin[iset]->Get("group2_"+hname ) );
      bg->Add( (TH1F*) fin[iset]->Get("group3_"+hname ) );
      bg->Add( (TH1F*) fin[iset]->Get("group4_"+hname ) );
      bg->Add( (TH1F*) fin[iset]->Get("group5_"+hname ) );
      bg->Add( (TH1F*) fin[iset]->Get("group6_"+hname ) );

      TH1F* total = (TH1F*) bg->Clone("sum");
      total->Add(sigNeg);
      total->Add(sigPos);

      sigNeg->SetLineColor(6);
      sigPos->SetLineColor(8);
      bg    ->SetLineColor(2);
      total ->SetLineColor(1);

      sigNeg->SetMarkerColor(6);
      sigPos->SetMarkerColor(8);
      bg    ->SetMarkerColor(2);
      total ->SetMarkerColor(1);

      int irebin=5;

      if ( scale!=1.0) {
        sigNeg->Scale(scale);
        sigPos->Scale(scale);
        bg    ->Scale(scale);
        total ->Scale(scale);
      }

      if ( irebin>1 ) {
        sigNeg->Rebin(irebin);
        sigPos->Rebin(irebin);
        bg    ->Rebin(irebin);
        total ->Rebin(irebin);
      }


      cc->Clear();

      total->SetTitle("");
      total->GetXaxis()->SetTitle("polarimeter");
      total->SetMinimum(0);
      total->Draw();
      bg->Draw("same");
      sigNeg->Draw("same");
      sigPos->Draw("same");

      tt.DrawTextNDC(0.55, 0.95, "IDR-"+lab[iset]+" ; tau -> "+dec);
      tt.DrawTextNDC(0.15, 0.95, "ILD preliminary");

      cc->Print(lab[iset]+"_"+dec+"_templates.eps");
      cc->Print("fit.pdf");


      float nevt=total->Integral();
      cout << nevt << endl;

      TH1F* htoy = (TH1F*) total->Clone("htoy");

      cc->Clear();
      cc->Divide(5,5);

      t1 = sigNeg;
      t2 = sigPos;
      tbg = bg;

      t1->Scale( 1./t1->Integral() );
      t2->Scale( 1./t2->Integral() );

      TF1 *func = new TF1("myfit",fitf,-1.,1.,2);
      func->SetParNames("total","negFrac");

      // TObjArray *mc = new TObjArray(2);
      // mc->Add( hsn );
      // mc->Add( hsp );
      // TFractionFitter* fit = new TFractionFitter(htoy, mc,"q");
      // fit->Constrain(0,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
      // fit->SetRangeX( total->GetXaxis()->FindBin(-1.), total->GetXaxis()->FindBin(1.) );

      for (int itoy=0; itoy<ntoyExps; itoy++) {

        htoy->Reset();

        // decide how many events
        int thisevt = trand3.Poisson(nevt);

        // fill histo with this # events randomly taken from total histo
        for (int i=0; i<thisevt; i++) {
          htoy->Fill( total->GetRandom() );
        }

        func->SetParameter(0, thisevt);
        //      func->FixParameter(0, thisevt);
        func->SetParameter(1, 0.5);
        func->SetParError(1, 0.5);

        TFitResultPtr r = htoy->Fit("myfit","Slq","",-1,1);

        hStatus [iset][idec]->Fill(r->Status());
        hFrac   [iset][idec]->Fill(r->Parameter(1));
        hFracErr[iset][idec]->Fill(r->ParError(1));
        hFitProb[iset][idec]->Fill(r->Prob());

      }

      cc->Clear();
      cc->Divide(2,2);
      cc->cd(1);      hStatus[iset] [idec]->Draw();
      cc->cd(2);      hFrac[iset]   [idec]->Draw();
      cc->cd(3);      

      float mm = hFracErr[iset][idec]->GetMean();
      float rr = hFracErr[iset][idec]->GetRMS();
      hFracErr[iset][idec]->GetXaxis()->SetRangeUser( 0.8*mm , 1.2*mm );
      hFracErr[iset][idec]->Draw();

      hFracErr[iset][idec]->Fit("gaus","ql");

      float mean =  hFracErr[iset][idec]->GetFunction("gaus")->GetParameter(1);

      hFracErr[iset][idec]->SetLineColor(icol);
      hFracErr[iset][idec]->SetLineStyle(ilty);
      hFracErr[iset][idec]->GetFunction("gaus")->SetLineColor(icol);

      TString ff; ff.Form("mean = %8.5f", mean );
      //      TString ff; ="mean="; ff+=hFracErr[iset][idec]->GetMean();
      tt.DrawTextNDC( 0.2, 0.8, ff );


      cc->cd(4);      hFitProb[iset][idec]->Draw();

      cout << hFracErr[iset][idec]->GetName() << " from gaus fit: " << mean << "  mean of histo: " << hFracErr[iset][idec]->GetMean() << endl;

      cc->Print("fit.pdf");

    }

  }
  for (int idec=0; idec<2; idec++) {
    cc->Clear();
    cc->Divide(2,2);
    for (int i=0; i<4; i++) {
      cc->cd(i+1);
      TH1F* h1 = hFracErr[i][idec];
      TH1F* h2 = hFracErr[i+4][idec];
      h1->Draw();
      h2->Draw("same");
    }
    cc->Print("fit.pdf");
  }


  cc->Print("fit.pdf]");


}
