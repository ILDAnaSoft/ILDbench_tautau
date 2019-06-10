#include "TF1.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TText.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TFractionFitter.h"
#include "TFitResult.h"

#include <fstream>
#include <iostream>
using std::cout;
using std::endl;


TH1F* tNeg;
TH1F* tPos;
TH1F* tbg;

Double_t fitf(Double_t *x,Double_t *par) {
  int ibin = tNeg->GetXaxis()->FindBin(x[0]);
  Double_t pNeg = tNeg->GetBinContent(ibin);
  Double_t pPos = tPos->GetBinContent(ibin);
  Double_t b  = tbg->GetBinContent(ibin);

  float pol=par[1];
  float fracNeg=(1.-pol)/2.;
  float fracPos=(1.+pol)/2.;

  Double_t pred = par[0]*( fracNeg*pNeg + fracPos*pPos ) + b;
  return pred;
}

int getNLS( TString cht, TString dec ) {

  //  TString cheatStrs[ncheat]={"exactCheatNoBG", "approxCheatNoBG", "approxCheatSELNoBG", "CheatGamSEL", "noCheatSEL"};
  //  TString decStrs[2] = { "pi", "rho"};

  int ILS = 2; //  both L and S

  if ( dec=="pi" ) {
    //    if (  cht == "CheatGamSELNoBG" || cht == "CheatGamSEL" ||  cht == "approxCheatNoBG" ) ILS=0; // dont consider
    if (  cht == "CheatGamSELNoBG" || cht == "CheatGamSEL" ) ILS=0; // dont consider
    else if (  cht == "exactCheatNoBG" || cht == "approxCheatNoBG" ) ILS=1;
  } else if ( dec=="rho" ) {
    if (  cht == "exactCheatNoBG" ||  cht == "approxCheatNoBG" ) ILS=1;
  } 

  return ILS;
}


void pexps( int ntoyExps=100 ) {

  cout << "nTOYS = " << ntoyExps << endl;

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

  const int ncheat=7;

  //  TString cheatStrs[ncheat]={"exactCheatNoBG", "approxCheatNoBG", "approxCheatSELNoBG", "CheatGamSELNoBG", "CheatGamSEL",  "noCheatSELNoBG", "noCheatSEL"};
  TString cheatStrs[ncheat]={"exactCheatNoBG", "approxCheatNoBG", "approxCheatSELNoBG", "CheatGamSELNoBG", "noCheatSELNoBG", "noCheatSEL"};
  TString decStrs[2] = { "pi", "rho"};


  float nevents[ncheat][nfile][2];

  float meanFitValue[ncheat][nfile][2];
  float meanFitValueErr[ncheat][nfile][2];
  float meanFitError[ncheat][nfile][2];

  float inputPol[ncheat][nfile][2];
  float pullMean[ncheat][nfile][2];
  float pullMeanErr[ncheat][nfile][2];
  float pullWidth[ncheat][nfile][2];
  float pullWidthErr[ncheat][nfile][2];

  for (int i=0; i<ncheat; i++) {
    for (int j=0; j<nfile; j++) {
      for (int k=0; k<2; k++) {
	nevents[i][j][k]=0;
	meanFitValue[i][j][k]=0;
	meanFitValueErr[i][j][k]=0;
	meanFitError[i][j][k]=0;
	inputPol[i][j][k]=0;
	pullMean[i][j][k]=0;
	pullMeanErr[i][j][k]=0;
	pullWidth[i][j][k]=0;
	pullWidthErr[i][j][k]=0;
      }
    }
  }


  TH1F* hPol   [ncheat][nfile][2];
  TH1F* hPolErr[ncheat][nfile][2];
  TH1F* hPolPull[ncheat][nfile][2];
  TH1F* hStatus [ncheat][nfile][2];
  TH1F* hFitProb[ncheat][nfile][2];

  TCanvas* cc = new TCanvas("cc","cc",400,400);
  cc->Print("fit.pdf[");

  TRandom3 trand3;

  float scale = 1; // SET TO 1 FOR FINAL RESULTS!!!!

  for (int icheat=0; icheat<ncheat; icheat++) {

    TString cheatSt = cheatStrs[icheat];

    cout << "ICHEAT " << icheat << " " << cheatSt << endl;

    for (int iset=0; iset<nfile; iset++) {
      
      int icol = iset<4 ? 4 : 2;
      int ilty = iset<4 ? 1 : 2;

      for (int idec=0; idec<2; idec++) {

	TString dec = decStrs[idec];

	hPol    [icheat][iset][idec] = new TH1F(cheatSt+lab[iset]+"_fitPol_"+dec,    cheatSt+lab[iset]+"_fitPol_"+dec,    1000, -1, 1);
	hPolErr [icheat][iset][idec] = new TH1F(cheatSt+lab[iset]+"_fitPolErr_"+dec, cheatSt+lab[iset]+"_fitPolErr_"+dec, 5000, 0., 0.05);
	hPolPull[icheat][iset][idec] = new TH1F(cheatSt+lab[iset]+"_fitPolPull_"+dec,cheatSt+lab[iset]+"_fitPolPull_"+dec,100, -5, 5);
	hStatus  [icheat][iset][idec] = new TH1F(cheatSt+lab[iset]+"_fitStatus_"+dec,  cheatSt+lab[iset]+"_fitStatus_"+dec,  10, 0, 10);
	hFitProb [icheat][iset][idec] = new TH1F(cheatSt+lab[iset]+"_fitProb_"+dec,    cheatSt+lab[iset]+"_fitProb_"+dec,    100, 0, 1);

	TH1F* sigNeg(0);
	TH1F* sigPos(0);
	TH1F* bg(0);
	if ( cheatStrs[icheat]== "noCheatSEL" ) {

	  TString hname="SEL_rec_" + dec + "_pol";

	  sigNeg = (TH1F*) (fin[iset]->Get("group0_"+hname+"_MCneg")->Clone("nchsn"));
	  sigPos = (TH1F*) (fin[iset]->Get("group0_"+hname+"_MCpos")->Clone("nchsp"));

	  bg = (TH1F*) (fin[iset]->Get("group0_"+hname+"_MCoth")->Clone("nchsb") );
	  bg->Add( (TH1F*) fin[iset]->Get("group1_"+hname ) );
	  bg->Add( (TH1F*) fin[iset]->Get("group2_"+hname ) );
	  bg->Add( (TH1F*) fin[iset]->Get("group3_"+hname ) );
	  bg->Add( (TH1F*) fin[iset]->Get("group4_"+hname ) );
	  bg->Add( (TH1F*) fin[iset]->Get("group5_"+hname ) );
	  bg->Add( (TH1F*) fin[iset]->Get("group6_"+hname ) );

	} else if ( cheatStrs[icheat]=="CheatGamSEL" ) {

	  if ( decStrs[idec]=="rho" ) {
	    TString hname="SEL_recCheatGam_" + dec + "_pol";

	    sigNeg = (TH1F*) ( fin[iset]->Get("group0_"+hname+"_MCneg") -> Clone("chgsn") );

	    sigPos = (TH1F*) ( fin[iset]->Get("group0_"+hname+"_MCpos")->Clone("chgsp") );

	    bg = (TH1F*) (fin[iset]->Get("group0_"+hname+"_MCoth")->Clone("chgsb") );
	    bg->Add( (TH1F*) fin[iset]->Get("group1_"+hname ) );
	    bg->Add( (TH1F*) fin[iset]->Get("group2_"+hname ) );
	    bg->Add( (TH1F*) fin[iset]->Get("group3_"+hname ) );
	    bg->Add( (TH1F*) fin[iset]->Get("group4_"+hname ) );
	    bg->Add( (TH1F*) fin[iset]->Get("group5_"+hname ) );
	    bg->Add( (TH1F*) fin[iset]->Get("group6_"+hname ) );

	  }

	} else if ( cheatStrs[icheat]== "noCheatSELNoBG" ) {

	  TString hname="SEL_rec_" + dec + "_pol";


	  sigNeg = (TH1F*) (fin[iset]->Get("group0_"+hname+"_MCneg")->Clone("nchsnbn") );
	  sigPos = (TH1F*) (fin[iset]->Get("group0_"+hname+"_MCpos")->Clone("nchsnbp") );

	  //	  sigPos = (TH1F*) fin[iset]->Get("group0_SEL_"+dec+"_mcPolarApprox_helPos");
	  bg =  (TH1F*) sigNeg->Clone("bg");
	  bg->Reset();

	} else if ( cheatStrs[icheat]== "CheatGamSELNoBG" ) {

	  if ( decStrs[idec]=="rho" ) {

	    TString hname="SEL_recCheatGam_" + dec + "_pol";

	    sigNeg = (TH1F*) (fin[iset]->Get("group0_"+hname+"_MCneg") -> Clone("chgsnbn") );
	    sigPos = (TH1F*) (fin[iset]->Get("group0_"+hname+"_MCpos") -> Clone("chgsnbp") );

	    //	  sigPos = (TH1F*) fin[iset]->Get("group0_SEL_"+dec+"_mcPolarApprox_helPos");
	    bg =  (TH1F*) sigNeg->Clone("bg");
	    bg->Reset();
	    
	  }

	} else if ( cheatStrs[icheat]== "approxCheatSELNoBG" ) {

	  sigNeg = (TH1F*) (fin[iset]->Get("group0_SEL_"+dec+"_mcPolarApprox_helNeg") -> Clone("achsnbn") );
	  sigPos = (TH1F*) (fin[iset]->Get("group0_SEL_"+dec+"_mcPolarApprox_helPos") -> Clone("achsnbp") );
	  bg =  (TH1F*) sigNeg->Clone("bg");
	  bg->Reset();

	} else if ( cheatStrs[icheat]=="exactCheatNoBG" ) {

	  sigNeg = (TH1F*) (fin[iset]->Get("group0_"+dec+"_mcPolarExact_helNeg") -> Clone("echgnbn") );
	  sigPos = (TH1F*) (fin[iset]->Get("group0_"+dec+"_mcPolarExact_helPos") -> Clone("echgnbp") );
	  bg =  (TH1F*) sigNeg->Clone("bg");
	  bg->Reset();

	} else if ( cheatStrs[icheat]=="approxCheatNoBG" ) {
	  
	  sigNeg = (TH1F*) (fin[iset]->Get("group0_"+dec+"_mcPolarApprox_helNeg") -> Clone("achgnbn") );
	  sigPos = (TH1F*) (fin[iset]->Get("group0_"+dec+"_mcPolarApprox_helPos") -> Clone("achgnbp") );
	  bg =  (TH1F*) sigNeg->Clone("bg");
	  bg->Reset();

	}

	if ( ! sigNeg ) continue;


	inputPol[icheat][iset][idec]= ( sigPos->Integral() - sigNeg->Integral() ) / ( sigPos->Integral() + sigNeg->Integral() );



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
	tt.DrawTextNDC(0.4, 0.85, cheatSt);

	cc->Print("figs/"+cheatSt+"_"+lab[iset]+"_"+dec+"_templates.eps");
	cc->Print("figs/"+cheatSt+"_"+lab[iset]+"_"+dec+"_templates.C");
	cc->Print("fit.pdf");


	float nevt=total->Integral();
	cout << "nevt " << nevt << endl;

	nevents[icheat][iset][idec] = nevt;

	TH1F* htoy = (TH1F*) total->Clone("htoy");

	cc->Clear();
	cc->Divide(5,5);

	tNeg = sigNeg;
	tPos = sigPos;
	tbg = bg;

	tNeg->Scale( 1./tNeg->Integral() );
	tPos->Scale( 1./tPos->Integral() );

	TF1 *func = new TF1("myfit",fitf,-1.,1.,2);
	func->SetParNames("total","polarisation");

	for (int itoy=0; itoy<ntoyExps; itoy++) {

	  htoy->Reset();

	  // decide how many events
	  int thisevt = trand3.Poisson(nevt);

	  // fill histo with this # events randomly taken from total histo
	  for (int i=0; i<thisevt; i++) {
	    htoy->Fill( total->GetRandom() );
	  }

	  func->SetParameter(0, thisevt);
	  func->SetParameter(1, 0.0);
	  func->SetParError(1, 0.5);

	  TFitResultPtr r = htoy->Fit("myfit","Slq","",-1,1);

	  float fitpol=r->Parameter(1); 
	  float fitpolErr = r->ParError(1);

	  //	  float fitpol=( (1.-jhjh) - jhjh );

	  hStatus [icheat][iset][idec]->Fill(r->Status());
	  hPol   [icheat][iset][idec]->Fill(fitpol);
	  hPolErr[icheat][iset][idec]->Fill(fitpolErr);
	  

	  hPolPull[icheat][iset][idec]->Fill( (fitpol -  inputPol[icheat][iset][idec])/fitpolErr );
	  hFitProb[icheat][iset][idec]->Fill(r->Prob());

	}


	cc->Clear();
	cc->Divide(3,2);
	cc->cd(1);      hStatus[icheat][iset][idec]->Draw();
	cc->cd(2);      

	float mm =  hPol  [icheat][iset][idec]->GetMean();
	float rr =  hPol  [icheat][iset][idec]->GetRMS();

	hPol  [icheat][iset][idec]->GetXaxis()->SetRangeUser(mm-5*rr, mm+5*rr );
	hPol  [icheat][iset][idec]->Draw();

	TLine* tl = new TLine(inputPol[icheat][iset][idec], 0, inputPol[icheat][iset][idec],  hPol  [icheat][iset][idec]->GetMaximum() );
	tl->SetLineColor(2);
	tl->Draw();

	cc->cd(3);      

	mm = hPolErr[icheat][iset][idec]->GetMean();
	rr = hPolErr[icheat][iset][idec]->GetRMS();

	hPolErr[icheat][iset][idec]->Fit("gaus","ql");

	float mean =  hPolErr[icheat][iset][idec]->GetFunction("gaus")->GetParameter(1);

	hPolErr[icheat][iset][idec]->SetLineColor(icol);
	hPolErr[icheat][iset][idec]->SetLineStyle(ilty);
	hPolErr[icheat][iset][idec]->GetFunction("gaus")->SetLineColor(icol);

	hPolErr[icheat][iset][idec]->GetXaxis()->SetRangeUser( mm-5*hPolErr[icheat][iset][idec]->GetRMS() , mm+5*hPolErr[icheat][iset][idec]->GetRMS() );
	hPolErr[icheat][iset][idec]->Draw();


	TString ff; ff.Form("mean = %8.5f", mean );
	//      TString ff; ="mean="; ff+=hPolErr[iset][idec]->GetMean();
	tt.DrawTextNDC( 0.2, 0.8, ff );


	cc->cd(4);      hFitProb[icheat][iset][idec]->Draw();

	cout << hPolErr[icheat][iset][idec]->GetName() << " from gaus fit: " << mean << "  mean of histo: " << hPolErr[icheat][iset][idec]->GetMean() << endl;

	meanFitValue[icheat][iset][idec] = hPol[icheat][iset][idec]->GetMean();
	meanFitValueErr[icheat][iset][idec] = hPol[icheat][iset][idec]->GetRMS()/sqrt(hPol[icheat][iset][idec]->GetEntries());
	meanFitError[icheat][iset][idec] = hPolErr[icheat][iset][idec]->GetMean();


	cc->cd(5); 

	hPolPull[icheat][iset][idec]->Fit("gaus","l;q");
	hPolPull[icheat][iset][idec]->Draw();

	pullMean    [icheat][iset][idec] = hPolPull[icheat][iset][idec]->GetFunction("gaus")->GetParameter(1);
	pullMeanErr [icheat][iset][idec] = hPolPull[icheat][iset][idec]->GetFunction("gaus")->GetParError(1);
	pullWidth   [icheat][iset][idec] = hPolPull[icheat][iset][idec]->GetFunction("gaus")->GetParameter(2);
	pullWidthErr[icheat][iset][idec] = hPolPull[icheat][iset][idec]->GetFunction("gaus")->GetParError(2);

	//	cout << "pull " << hPolPull[icheat][iset][idec]->GetFunction("gaus")->GetParameter(1) << " " <<  hPolPull[icheat][iset][idec]->GetFunction("gaus")->GetParameter(2) << endl;

	cc->Print("fit.pdf");


	//cout << "rikka1" << endl;

      }

      //cout << "rikka2" << endl;

    }

    //cout << "rikka3" << endl;

    for (int idec=0; idec<2; idec++) {
      cc->Clear();
      cc->Divide(2,2);
      for (int i=0; i<4; i++) { // 4 pols, *2 detectors
	cc->cd(i+1);
	TH1F* h1 = hPolErr[icheat][i][idec];   // large
	TH1F* h2 = hPolErr[icheat][i+4][idec]; // small

	float mmx = std::max( h1->GetMean()+5*h1->GetRMS(),  h2->GetMean()+5*h2->GetRMS() );
	float mmn = std::min( h1->GetMean()-5*h1->GetRMS(),  h2->GetMean()-5*h2->GetRMS() );

	h1->GetXaxis()->SetRangeUser( mmn, mmx );

	// cout << icheat << " " << i << " " << idec << " " << h1 << " " << h2 << endl;

	h1->Draw();
	h2->Draw("same");
      }
      cc->Print("fit.pdf");
    }

    // cout << "rikka4" << endl;

  } // cheat

  cc->Print("fit.pdf]");



  std::ofstream fitresultsfile;
  fitresultsfile.open("fitresults.txt");

  //  cout << "rikka5" << endl;
  for (int idec=0; idec<2; idec++) {
    fitresultsfile << decStrs[idec] << endl;

    for (int iverb=0; iverb<2; iverb++) {
      for (int iset=0; iset<nfile/2; iset++) {
	if ( iverb==0 ) fitresultsfile << " set " << iset << " : ";
	//      for (int idec=0; idec<2; idec++) {
	if ( iverb==0 ) fitresultsfile << decStrs[idec] << " ";
	for (int icheat=0; icheat<ncheat; icheat++) {
	  if ( iverb==0 ) fitresultsfile << cheatStrs[icheat] << " ( ";
	  
	  //	  int ILS= cheatStrs[icheat]== "noCheat" || cheatStrs[icheat]== "CheatGam" ? 2 : 1;
	  int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
	  
	  int ibest = meanFitError[icheat][iset+4][idec] < meanFitError[icheat][iset][idec] ? 1 : 0 ;
	  
	  for (int ils=0; ils<ILS; ils++) {
	    
	    TString col = "black";
	    
	    if ( cheatStrs[icheat]== "noCheatSEL" && ils==ibest) {
	      col=  ils==0 ? "blue" : "red";
	    }
	    TString hh; hh.Form(" & \\textcolor{%5s}{%3.1f}", col.Data(), 100*meanFitError[icheat][iset+ils*4][idec] );
	    fitresultsfile << hh;
	  }
	  if ( iverb==0 ) fitresultsfile << " ) ";
	  //}
	}
	
	for (int ils=0; ils<2; ils++) {
	  fitresultsfile << 100./ sqrt( pow(1./meanFitError[2][iset+ils*4][0],2) + pow(1./meanFitError[2][iset+ils*4][1],2) ) << " ";
	}
      
	fitresultsfile << " \\\\ " << endl;
      }
    }


    fitresultsfile << "\\multicolumn{9}{|c|}{INPUT POL \\%} \\\\" << endl;
    for (int iset=0; iset<nfile/2; iset++) { // pol set
      //    for (int idec=0; idec<2; idec++) { // tau decays
      for (int icheat=0; icheat<ncheat; icheat++) { // cheating
	int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
	for (int ils=0; ils<ILS; ils++) { // L/S
	  TString hh; hh.Form(" & $ %5.1f $", 100*inputPol[icheat][iset+ils*4][idec] );
	  fitresultsfile << hh;
	}
	//      }
      }
      fitresultsfile << " \\\\ " << endl;
    }

    fitresultsfile << "\\multicolumn{9}{|c|}{NEVENTS} \\\\" << endl;
    for (int iset=0; iset<nfile/2; iset++) { // pol set
      //    for (int idec=0; idec<2; idec++) { // tau decays
      for (int icheat=0; icheat<ncheat; icheat++) { // cheating
	int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
	for (int ils=0; ils<ILS; ils++) { // L/S
	  TString hh; hh.Form(" & $ %7.0f $", nevents[icheat][iset+ils*4][idec] );
	  fitresultsfile << hh;
	}
	//      }
      }
      fitresultsfile << " \\\\ " << endl;
    }
    
    fitresultsfile << "\\multicolumn{9}{|c|}{FITTED POL \\%} \\\\" << endl;
    for (int iset=0; iset<nfile/2; iset++) {
      //    for (int idec=0; idec<2; idec++) {
      for (int icheat=0; icheat<ncheat; icheat++) {
	int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
	for (int ils=0; ils<ILS; ils++) {
	  TString hh; hh.Form(" & $ %5.1f \\pm  %5.1f $", 100*meanFitValue[icheat][iset+ils*4][idec], 100*meanFitValueErr[icheat][iset+ils*4][idec] );
	  fitresultsfile << hh;
	}
	//      }
      }
      fitresultsfile << " \\\\ " << endl;
    }

    fitresultsfile << "\\multicolumn{9}{|c|}{MEAN FITTED POL Error \\%} \\\\" << endl;
    for (int iset=0; iset<nfile/2; iset++) {
      //    for (int idec=0; idec<2; idec++) {
      for (int icheat=0; icheat<ncheat; icheat++) {
	int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
	for (int ils=0; ils<ILS; ils++) {
	  TString hh; hh.Form(" & $ %4.2f $", 100*meanFitError[icheat][iset+ils*4][idec] );
	  fitresultsfile << hh;
	}
	//      }
      }
      fitresultsfile << " \\\\ " << endl;
    }

    fitresultsfile << "\\multicolumn{9}{|c|}{PULL MEAN}\\\\" << endl;
    for (int iset=0; iset<nfile/2; iset++) {
      //    for (int idec=0; idec<2; idec++) {
      for (int icheat=0; icheat<ncheat; icheat++) {
	int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
	for (int ils=0; ils<ILS; ils++) {
	  TString hh; hh.Form(" & $ %5.2f \\pm  %5.2f $", pullMean    [icheat][iset+ils*4][idec], pullMeanErr [icheat][iset+ils*4][idec] );
	  fitresultsfile << hh;
	}
	//      }
      }
      fitresultsfile << " \\\\ " << endl;
    }

    fitresultsfile << "\\multicolumn{9}{|c|}{PULL WIDTH}\\\\" << endl;
    for (int iset=0; iset<nfile/2; iset++) {
      //    for (int idec=0; idec<2; idec++) {
      for (int icheat=0; icheat<ncheat; icheat++) {
	int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
	for (int ils=0; ils<ILS; ils++) {
	  TString hh; hh.Form(" & $ %5.2f \\pm %5.2f $", pullWidth    [icheat][iset+ils*4][idec], pullWidthErr [icheat][iset+ils*4][idec] );
	  fitresultsfile << hh;
	}
	//      }
      }
      fitresultsfile << " \\\\ " << endl;
    }
  }

  fitresultsfile.close();


  // draw plots

  cout << "hi1" << endl;

  int idec=0;
  float pi_prec[4][8];
  for (int iset=0; iset<nfile/2; iset++) {
    int jj=0;
    for (int icheat=0; icheat<ncheat; icheat++) {
      int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
      for (int ils=0; ils<ILS; ils++) {
	pi_prec[iset][jj++]=100*meanFitError[icheat][iset+ils*4][idec];
      }
    }
  }

  cout << "hi2" << endl;

  idec=1;
  float rho_prec[4][10];
  for (int iset=0; iset<nfile/2; iset++) {
    int jj=0;
    for (int icheat=0; icheat<ncheat; icheat++) {
      int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
      for (int ils=0; ils<ILS; ils++) {
	rho_prec[iset][jj++]=100*meanFitError[icheat][iset+ils*4][idec];
      }
    }
  }

  cout << "hi3" << endl;

  gStyle->SetMarkerSize(1.5);
  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetPadRightMargin(0.23);

  const int npP=5;
  const int npR=6;

  float pi_prec_L[4][npP];
  float pi_prec_S[4][npP];

  float rho_prec_L[4][npR];
  float rho_prec_S[4][npR];

  TGraph* gr_pi_L[4];
  TGraph* gr_pi_S[4];

  TGraph* gr_rho_L[4];
  TGraph* gr_rho_S[4];

  float indexP[npP] = {1,2,3,5,6};
  float indexR[npR] = {1,2,3,4,5,6};

  TString axLabs[npR] = {"MC-optimal", "MC-approx (no Sel/no Bg)", "MC-approx (Sel/no Bg)", "cheat ECAL (Sel/no Bg)", "no cheat (Sel/no Bg)", "no cheat (Sel/Bg)"};

  for (int i=0; i<4; i++) {

    cout << "hi4" << endl;

    pi_prec_L[i][0] = pi_prec[i][0];
    pi_prec_L[i][1] = pi_prec[i][1];
    pi_prec_L[i][2] = pi_prec[i][2];
    pi_prec_L[i][3] = pi_prec[i][4];
    pi_prec_L[i][4] = pi_prec[i][6];

    gr_pi_L[i] = new TGraph(npP, indexP, &(pi_prec_L[i][0]) );
    gr_pi_L[i]->SetLineColor(4);
    gr_pi_L[i]->SetMarkerColor(4);
    gr_pi_L[i]->SetMarkerStyle(23);

    pi_prec_S[i][0] = pi_prec[i][0];
    pi_prec_S[i][1] = pi_prec[i][1];
    pi_prec_S[i][2] = pi_prec[i][3];
    pi_prec_S[i][3] = pi_prec[i][5];
    pi_prec_S[i][4] = pi_prec[i][7];

    gr_pi_S[i] = new TGraph(npP, indexP, &(pi_prec_S[i][0]) );
    gr_pi_S[i]->SetLineColor(2);
    gr_pi_S[i]->SetLineStyle(2);
    gr_pi_S[i]->SetMarkerColor(2);
    gr_pi_S[i]->SetMarkerStyle(23);

    rho_prec_L[i][0] = rho_prec[i][0];
    rho_prec_L[i][1] = rho_prec[i][1];
    rho_prec_L[i][2] = rho_prec[i][2];
    rho_prec_L[i][3] = rho_prec[i][4];
    rho_prec_L[i][4] = rho_prec[i][6];
    rho_prec_L[i][5] = rho_prec[i][8];

    gr_rho_L[i] = new TGraph(npR, indexR, &(rho_prec_L[i][0]) );
    gr_rho_L[i]->SetLineColor(4);
    gr_rho_L[i]->SetMarkerColor(4);
    gr_rho_L[i]->SetMarkerStyle(24);

    rho_prec_S[i][0] = rho_prec[i][0];
    rho_prec_S[i][1] = rho_prec[i][1];
    rho_prec_S[i][2] = rho_prec[i][3];
    rho_prec_S[i][3] = rho_prec[i][5];
    rho_prec_S[i][4] = rho_prec[i][7];
    rho_prec_S[i][5] = rho_prec[i][9];

    gr_rho_S[i] = new TGraph(npR, indexR, &(rho_prec_S[i][0]) );
    gr_rho_S[i]->SetLineColor(2);
    gr_rho_S[i]->SetLineStyle(2);
    gr_rho_S[i]->SetMarkerColor(2);
    gr_rho_S[i]->SetMarkerStyle(24);
  }

  cout << "hi5" << endl;

  TH2F* hh = new TH2F("prec","",npR,0.5,npR+0.5,10,0,3);
  hh->GetYaxis()->SetTitle("polarisation precision [%]");

  for (int i=1; i<=npR; i++) {
    hh->GetXaxis()->SetBinLabel(i, axLabs[i-1]);
  }


  cc = new TCanvas("c2","c2",400,400);
  cc->Print("drawpexpresults.pdf[");

  TText ttf;
  ttf.SetTextFont(62);
  ttf.SetTextSize(0.06);


  TLegend* tl=0;

  for (int i=0; i<4; i++) {
    cc->Clear();

    hh->Draw();

    TString ttx;
    if      (i==0 ) {hh->GetYaxis()->SetRangeUser(0, 1. ); ttx="eLpR";}
    else if (i==1 ) {hh->GetYaxis()->SetRangeUser(0, 1. ); ttx="eRpL";}
    else if (i==2 ) {hh->GetYaxis()->SetRangeUser(0, 3. ); ttx="eLpL";}
    else if (i==3 ) {hh->GetYaxis()->SetRangeUser(0, 3. ); ttx="eRpR";}

    gr_rho_L[i]->Draw("pl");
    gr_rho_S[i]->Draw("pl");
    gr_pi_L[i]->Draw("pl");
    gr_pi_S[i]->Draw("pl");

    tt.DrawTextNDC(0.25, 0.85, ttx);

    if (!tl) {
      tl=new TLegend(0.4, 0.25, 0.73, 0.45);
      tl->AddEntry( gr_pi_L[i] , "#tau #rightarrow #pi IDR-L", "pl" );
      tl->AddEntry( gr_pi_S[i] , "#tau #rightarrow #pi IDR-S", "pl" );
      tl->AddEntry( gr_rho_L[i], "#tau #rightarrow #rho IDR-L", "pl" );
      tl->AddEntry( gr_rho_S[i], "#tau #rightarrow #rho IDR-S", "pl" );
    }

    tl->Draw();

    ttf.DrawTextNDC(0.38, 0.95, "ILD preliminary");


    TString dd="drawpexpresults"; dd+=i;
    cc->Print("figs/"+dd+".eps");
    cc->Print("figs/"+dd+".C");

    cc->Print("drawpexpresults.pdf");

  }

  cc->Print("drawpexpresults.pdf]");

  

}
