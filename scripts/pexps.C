#include "TF1.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TText.h"
#include "TLatex.h"
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


bool useRaw;

TH1F* tNeg;
TH1F* tPos;
TH1F* tbg;

TH1F* tNegFn;
TH1F* tPosFn;
TH1F* tBgFn;

float fitRangeLow;
float fitRangeUp;

Double_t fitf(Double_t *x,Double_t *par) {

  //  if ( x[0]<fitRangeLow || x[0]>fitRangeUp ) return 0;

  TH1F* hp = useRaw ? tPos : tPosFn;
  TH1F* hn = useRaw ? tNeg : tNegFn;
  TH1F* hb = useRaw ? tbg : tBgFn;

  int ibin =  hp->GetXaxis()->FindBin(x[0]);

  Double_t pNeg = hn->GetBinContent(ibin);
  Double_t pPos = hp->GetBinContent(ibin);
  Double_t b    = hb->GetBinContent(ibin);

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

  useRaw = false;
  cout << "using raw histos ? " << useRaw << endl;

  const int irebin=1;

  TString indir = "data_ALLEVToneHad/";

  fitRangeLow=-1.; fitRangeUp=1.;

  TString label = "fit";
  if (useRaw)   label+="RAW";
  else          label+="FIT";
  label+="_RebinNN";
  label+=irebin;

  TString outdir = indir + label + "/";


  TString plotfile = outdir + label + ".pdf";


  const int nfile=8;
  // const int nfile=2;
  TFile* fin[nfile];
  
  // fin[0] = new TFile(indir + "ILD_l5_o1_v02_eL80pR30/fithistos.root","read");
  // fin[1] = new TFile(indir + "ILD_s5_o1_v02_eL80pR30/fithistos.root","read");
  // TString lab[nfile]={"L_lr", "S_lr"};

  fin[0] = new TFile(indir + "ILD_l5_o1_v02_eL80pR30/fithistos.root","read");
  fin[1] = new TFile(indir + "ILD_l5_o1_v02_eR80pL30/fithistos.root","read");
  fin[2] = new TFile(indir + "ILD_l5_o1_v02_eL80pL30/fithistos.root","read");
  fin[3] = new TFile(indir + "ILD_l5_o1_v02_eR80pR30/fithistos.root","read");
  fin[4] = new TFile(indir + "ILD_s5_o1_v02_eL80pR30/fithistos.root","read");
  fin[5] = new TFile(indir + "ILD_s5_o1_v02_eR80pL30/fithistos.root","read");
  fin[6] = new TFile(indir + "ILD_s5_o1_v02_eL80pL30/fithistos.root","read");
  fin[7] = new TFile(indir + "ILD_s5_o1_v02_eR80pR30/fithistos.root","read");
  TString lab[nfile]={"L_lr", "L_rl", "L_ll", "L_rr", "S_lr", "S_rl", "S_ll", "S_rr"};

  TString LAB[nfile]={
    "IDR-L e^{#minus}_{L80}e^{#plus}_{R30}", 
    "IDR-L e^{#minus}_{R80}e^{#plus}_{L30}", 
    "IDR-L e^{#minus}_{L80}e^{#plus}_{L30}", 
    "IDR-L e^{#minus}_{R80}e^{#plus}_{R30}", 
    "IDR-S e^{#minus}_{L80}e^{#plus}_{R30}", 
    "IDR-S e^{#minus}_{R80}e^{#plus}_{L30}", 
    "IDR-S e^{#minus}_{L80}e^{#plus}_{L30}", 
    "IDR-S e^{#minus}_{R80}e^{#plus}_{R30}" };


  const int npols = nfile/2; // assumes we have LS for each polarisation!

  TString pollab[npols]={"eLpR","eRpL","eLpL", "eRpR"};
  // TString pollab[npols]={"eLpR"};




  tNegFn=0;
  tPosFn=0;
  tBgFn=0;

  TText tt;
  TLatex tltx;

  const int ncheat=6;
  TString cheatStrs[ncheat]={"exactCheatNoBG", "approxCheatNoBG", "approxCheatSELNoBG", "CheatGamSELNoBG", "noCheatSELNoBG", "noCheatSEL"};

  TString cheatLabs[ncheat]={"MC-optimal","MC-approx (no Sel/no Bg)", "MC-approx (Sel/no Bg)",
			     "cheat ECAL (Sel/no Bg)", "no cheat (Sel/no Bg)", "no cheat (Sel/Bg)"};


  //const int ncheat=2;
  //TString cheatStrs[ncheat]={"exactCheatNoBG", "CheatGamSELNoBG"};
  //const int ncheat=2;
  //TString cheatStrs[ncheat]={"exactCheatNoBG", "noCheatSEL"};

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

  float maxmeanfiterror[npols]={0};

  TCanvas* cc = new TCanvas("cc","cc",400,400);
  cc->Print(plotfile+"[");

  TRandom3 trand3;

  float neventscale = 1.; // scale the int lumi [SET TO 1 FOR FINAL RESULTS!!!!]
  
  TF1 *func = new TF1("myfit",fitf,-2.,2.,2);
  func->SetParNames("total","polarisation");



  TF1 * ff1 = new TF1("ff1", "gaus(0)+pol3(3)", -2., 2.);
  TF1 * ff2 = new TF1("ff2", "gaus(0) + gaus(3) + [6]*sin( (x + 1)*TMath::Pi()/2 ) + [7]*sin( 2*(x + 1)*TMath::Pi()/2 ) + [8]*sin( 3*(x + 1)*TMath::Pi()/2 ) + [9]*sin( 4*(x + 1)*TMath::Pi()/2 ) + [10]*sin( 5*(x + 1)*TMath::Pi()/2 ) ", -2., 2.);
  TF1 * ff3 = new TF1("ff3", "gaus(0) + [3]*sin( (x + 1)*TMath::Pi()/2 ) + [4]*sin( 2*(x + 1)*TMath::Pi()/2 ) + [5]*sin( 3*(x + 1)*TMath::Pi()/2 ) + [6]*sin( 4*(x + 1)*TMath::Pi()/2 ) + [7]*sin( 5*(x + 1)*TMath::Pi()/2 ) ", -2., 2.);

  for (int icheat=0; icheat<ncheat; icheat++) { // cheating scenarios
    TString cheatSt = cheatStrs[icheat];
    TString cheatLb = cheatLabs[icheat];
    for (int iset=0; iset<nfile; iset++) { // polarisations
      int icol = iset<npols ? 4 : 2;
      int ilty = iset<npols ? 1 : 2;
      for (int idec=0; idec<2; idec++) { // tau decay modes


	cout << "ICHEAT " << icheat << " " << cheatSt << " set " << iset << " dec " << idec << endl;
	
	// define the fit range
	if ( idec==0 && cheatStrs[icheat].Contains("noCheat") ) { // uncheated pi -> restrict range
	  fitRangeLow=-0.9;
	  fitRangeUp=0.9;
	} else { // rho
	  fitRangeLow=-1;
	  fitRangeUp=1;
	}

	// define histograms

	TString dec = decStrs[idec];
	hPol    [icheat][iset][idec] = new TH1F(cheatSt+lab[iset]+"_fitPol_"+dec,    cheatSt+lab[iset]+"_fitPol_"+dec,    1000, -1, 1);

	hPolErr [icheat][iset][idec] = new TH1F(cheatSt+lab[iset]+"_fitPolErr_"+dec, cheatSt+lab[iset]+"_fitPolErr_"+dec, 5000, 0., 0.05/sqrt(neventscale) );
	hPolPull[icheat][iset][idec] = new TH1F(cheatSt+lab[iset]+"_fitPolPull_"+dec,cheatSt+lab[iset]+"_fitPolPull_"+dec,400, -10, 10);
	hStatus  [icheat][iset][idec] = new TH1F(cheatSt+lab[iset]+"_fitStatus_"+dec,  cheatSt+lab[iset]+"_fitStatus_"+dec,  10, 0, 10);
	hFitProb [icheat][iset][idec] = new TH1F(cheatSt+lab[iset]+"_fitProb_"+dec,    cheatSt+lab[iset]+"_fitProb_"+dec,    100, 0, 1);


	// get the raw histos from input file
	//  group them as appropriate

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

	  bg =  (TH1F*) sigNeg->Clone("bg");
	  bg->Reset();

	} else if ( cheatStrs[icheat]== "CheatGamSELNoBG" ) {

	  if ( decStrs[idec]=="rho" ) {

	    TString hname="SEL_recCheatGam_" + dec + "_pol";

	    sigNeg = (TH1F*) (fin[iset]->Get("group0_"+hname+"_MCneg") -> Clone("chgsnbn") );
	    sigPos = (TH1F*) (fin[iset]->Get("group0_"+hname+"_MCpos") -> Clone("chgsnbp") );

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

	// check all is OK
	if ( ! sigNeg ) continue;


	// now get the number of input events, and MC truth polarisation

	//	inputPol[icheat][iset][idec]= ( sigPos->Integral() - sigNeg->Integral() ) / ( sigPos->Integral() + sigNeg->Integral() );
	//float npos = sigPos->Integral(-1,1);
	//float nneg = sigNeg->Integral(-1,1);
	//float npos = sigPos->Integral(  sigPos->FindBin(fitRangeLow), sigPos->FindBin(fitRangeUp) );
	//float nneg = sigNeg->Integral(  sigNeg->FindBin(fitRangeLow), sigNeg->FindBin(fitRangeUp) );

	// here we get actual total events (even beyond fitting limits)
	float npos = sigPos->Integral( );
	float nneg = sigNeg->Integral( );
	inputPol[icheat][iset][idec]= (npos-nneg) / (npos+nneg);

	cout << "input polarisation: " << inputPol[icheat][iset][idec] << " , " << npos << " "  << nneg << endl;

	// sum the templates to get the total
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

	// rebin if requested
	if ( irebin>1 ) {
	  sigNeg->Rebin(irebin);
	  sigPos->Rebin(irebin);
	  bg    ->Rebin(irebin);
	  total ->Rebin(irebin);
	}


	// draw the input templates
	cc->Clear();
	total->SetTitle("");
	total->GetXaxis()->SetTitle("polarimeter");
	total->SetMinimum(0);
	total->SetMaximum( 1.25*total->GetMaximum() );
	total->Draw();
	bg->Draw("same");
	sigNeg->Draw("same");
	sigPos->Draw("same");
	
	// tt.DrawTextNDC(0.55, 0.95, "IDR-"+lab[iset]+" ; tau -> "+dec);
	// tt.DrawTextNDC(0.15, 0.95, "ILD preliminary");
	// tt.DrawTextNDC(0.25, 0.85, cheatLb);
	tt.DrawTextNDC(0.15, 0.95, "ILD preliminary");
	tltx.DrawLatexNDC(0.48, 0.95, LAB[iset]+"  #tau #rightarrow #"+dec);
	tt.DrawTextNDC(0.35, 0.85, cheatLb);

	cc->Print(outdir+cheatSt+"_"+lab[iset]+"_"+dec+"_rawtemplates.eps");
	cc->Print(outdir+cheatSt+"_"+lab[iset]+"_"+dec+"_rawtemplates.C");

	cc->Print(plotfile);


	// now try to fit the input templates
	TF1* fPos;
	TF1* fNeg;
	TF1* fBg;

	if ( idec==0 || cheatStrs[icheat].Contains("exact") ) { // pi OR exact -> use pol1

	  // do this shape fit in restricted range to remove edge effects...

	  sigPos->Fit("pol1","q","",-0.9, 0.9);
	  sigNeg->Fit("pol1","q","",-0.9, 0.9);
	  bg    ->Fit("pol1","q","",-0.9, 0.9);

	  // sigPos->Fit("pol1","q","",fitRangeLow, fitRangeUp);
	  // sigNeg->Fit("pol1","q","",fitRangeLow, fitRangeUp);
	  // bg->Fit("pol1","q","",fitRangeLow, fitRangeUp);

	  fPos = sigPos->GetFunction( "pol1" );
	  //if (fPos) fPos->SetLineColor( sigPos->GetLineColor() );
	  if (fPos) fPos->SetLineColor( 1 );

	  fNeg = sigNeg->GetFunction( "pol1" );
	  //if ( fNeg) fNeg->SetLineColor( sigNeg->GetLineColor() );
	  if ( fNeg) fNeg->SetLineColor( 1 );

	  fBg  = bg->GetFunction( "pol1" );
	  if ( fBg ) fBg->SetLineColor( 1 );
	  //if ( fBg ) fBg->SetLineColor(  bg->GetLineColor() );

	} else if ( idec==1 ) { // rho -> use more complex fitting function
	  
	  ff3->SetParameter(0, sigPos->GetMaximum()/2 );
	  ff3->SetParLimits(0, 0, sigPos->GetMaximum() );
	  ff3->SetParameter(1, -0.5 );
	  ff3->SetParLimits(1, -0.6, -0.4);
	  //ff3->SetParameter(2, 0.1 );
	  //ff3->SetParLimits(2, 0.02, 0.2);
	  ff3->SetParameter(2, 0.2 );
	  ff3->SetParLimits(2, 0.1, 0.5);
	  
	  sigPos->Fit("ff3","q","",fitRangeLow,fitRangeUp);

	  fPos = sigPos->GetFunction( "ff3" );
	  // if ( fPos) fPos->SetLineColor( sigPos->GetLineColor() );
	  if ( fPos) fPos->SetLineColor( 1 );

	  ff2->SetParameter(0, sigNeg->GetMaximum()/2 );
	  ff2->SetParLimits(0, 0, sigNeg->GetMaximum() );
	  ff2->SetParameter(1, -0.5 );
	  ff2->SetParLimits(1, -0.6, -0.4);
	  ff2->SetParameter(2, 0.1 );
	  ff2->SetParLimits(2, 0.02, 0.2);

	  ff2->SetParameter(3, sigPos->GetMaximum()/2 );
	  ff2->SetParLimits(3, 0, sigPos->GetMaximum() );
	  ff2->SetParameter(4, -0.5 );
	  ff2->SetParLimits(4, -0.6, -0.4);
	  ff2->SetParameter(5, 0.3 );
	  ff2->SetParLimits(5, 0.1, 0.5);
	  

	  sigNeg->Fit("ff2","q","",fitRangeLow,fitRangeUp);

	  fNeg = sigNeg->GetFunction( "ff2" );
	  //if ( fNeg) fNeg->SetLineColor( sigNeg->GetLineColor() );
	  if ( fNeg) fNeg->SetLineColor( 1 );

	  if ( bg->GetEntries()>0 ) {
	    ff3->SetParameter(0, bg->GetMaximum()/2 );
	    ff3->SetParLimits(0, 0, bg->GetMaximum() );
	    ff3->SetParameter(1, -0.5 );
	    ff3->SetParLimits(1, -0.6, -0.4);
	    ff3->SetParameter(2, 0.1 );
	    ff3->SetParLimits(2, 0.02, 0.2);
	    bg->Fit("ff3","q","",fitRangeLow,fitRangeUp);
	    fBg  = bg->GetFunction( "ff3" );
	    //if ( fBg) fBg->SetLineColor(  bg->GetLineColor() );
	    if ( fBg) fBg->SetLineColor( 1 );
	  }

	}


	// draw the fitted templates
	cc->Clear();
	total->Draw();
	bg->Draw("same");
	sigNeg->Draw("same");
	sigPos->Draw("same");

	tt.DrawTextNDC(0.15, 0.95, "ILD preliminary");
	tltx.DrawLatexNDC(0.48, 0.95, LAB[iset]+"  #tau #rightarrow #"+dec);
	tt.DrawTextNDC(0.35, 0.85, cheatLb);

	cc->Print(outdir+cheatSt+"_"+lab[iset]+"_"+dec+"_templates.eps");
	cc->Print(outdir+cheatSt+"_"+lab[iset]+"_"+dec+"_templates.C");
	cc->Print(plotfile);


	// now we get the total number of events in the fitting range

	float nevt_fitRange= neventscale * total->Integral( total->FindBin(fitRangeLow), total->FindBin(fitRangeUp) );
	nevents[icheat][iset][idec] = nevt_fitRange;


	// set the templates

	tNeg = sigNeg;
	tPos = sigPos;
	tbg = bg;


	// normalise them to have integral 1.000  -> pdf

	//tNeg->Scale( 1./tNeg->Integral( tNeg->FindBin(fitRangeLow), tNeg->FindBin(fitRangeUp) ) ); // this was nonsense
	//tPos->Scale( 1./tPos->Integral( tPos->FindBin(fitRangeLow), tPos->FindBin(fitRangeUp) ) );
	tNeg->Scale( 1./tNeg->Integral( ) );
	tPos->Scale( 1./tPos->Integral( ) );
	tbg->Scale( neventscale );  // except the bg, which we leave at full stats -> this is not floated in the fit

	// --- here we make templates from the fitted functions
	//  smoothed versions of the original
	if ( !tNegFn ) {
	  tNegFn = (TH1F*) tNeg->Clone("tNegFn");
	  tPosFn = (TH1F*) tNeg->Clone("tPosFn");
	  tBgFn = (TH1F*) tNeg->Clone("tBgFn");
	}

	for (int jb=1; jb<=tNegFn->GetNbinsX(); jb++) {
	  if (  tNegFn->GetBinCenter(jb) < fitRangeLow ||  tNegFn->GetBinCenter(jb) > fitRangeUp ) { // set the pdf to 0 outside the range
	    tNegFn->SetBinContent( jb, 0. );
	    tPosFn->SetBinContent( jb, 0. );
	    tBgFn->SetBinContent( jb, 0. );
	  } else {
	    if ( fNeg->Eval( tNegFn->GetBinCenter(jb) ) > 0 ) {
	      tNegFn->SetBinContent( jb, fNeg->Eval( tNegFn->GetBinCenter(jb) ) );
	    } else {
	      tNegFn->SetBinContent( jb, 0 );
	    }

	    if (  fPos->Eval( tPosFn->GetBinCenter(jb) ) > 0 ) {
	      tPosFn->SetBinContent( jb, fPos->Eval( tPosFn->GetBinCenter(jb) ) );
	    } else {
	      tPosFn->SetBinContent( jb, 0 );
	    }

	    if ( bg->GetEntries()>0 && fBg->Eval( tBgFn->GetBinCenter(jb) ) > 0 ) {
	      tBgFn->SetBinContent( jb, fBg->Eval( tBgFn->GetBinCenter(jb) ) );
	    } else {
	      tBgFn->SetBinContent( jb, 0. );
	    }
	  }
	}

	
	TH1F* totalFn = (TH1F*) tNegFn->Clone("totalFn");
	totalFn->Add( tPosFn );
	totalFn->Add( tBgFn );

	totalFn->SetLineColor(1);
	totalFn->SetMarkerColor(1);

	tBgFn->SetLineColor(2);
	tBgFn->SetMarkerColor(2);

	tNegFn->SetLineColor(6);
	tNegFn->SetMarkerColor(6);

	tPosFn->SetLineColor(8);
	tPosFn->SetMarkerColor(8);

	cc->Clear();
	totalFn->SetMinimum(0);
	//	totalFn->GetXaxis()->SetRangeUser(fitRangeLow,fitRangeUp);
	totalFn->Draw("hist");
	tBgFn ->Draw("hist;same");
	tNegFn->Draw("hist;same");
	tPosFn->Draw("hist;same");
	cc->Print(plotfile);


	// now we need to be careful what we scale to!
	// I think the integral of the smoother histo in the fitting range should be equal to the same integral of the unsmoothed normalised distribution

	float intNeg = tNeg->Integral( tNeg->FindBin(fitRangeLow), tNeg->FindBin(fitRangeUp) );
	float intPos = tPos->Integral( tPos->FindBin(fitRangeLow), tPos->FindBin(fitRangeUp) );


	//tNegFn->Scale( 1./tNegFn->Integral( tNegFn->FindBin(fitRangeLow), tNegFn->FindBin(fitRangeUp) ) ); // this was nonsense
	//tPosFn->Scale( 1./tPosFn->Integral( tPosFn->FindBin(fitRangeLow), tPosFn->FindBin(fitRangeUp) ) );
	//tNegFn->Scale( 1./tNegFn->Integral( ) );
	//tPosFn->Scale( 1./tPosFn->Integral( ) );
	tNegFn->Scale( intNeg/tNegFn->Integral( ) );
	tPosFn->Scale( intPos/tPosFn->Integral( ) );
	tBgFn->Scale( neventscale ); 

	cc->Clear();
	cc->Divide(2,2);
	cc->cd(1);
	tNeg->Draw("hist");
	tNegFn->Draw("hist;same");
	cc->cd(2);
	tPos->Draw("hist");
	tPosFn->Draw("hist;same");
	cc->cd(3);
	tbg->Draw("hist");
	tBgFn->Draw("hist;same");
	cc->Print(plotfile);

	// -----------------------------
	// now start to do the fits

	TH1F* hclone[25];
	TH1F* htoy = (TH1F*) total->Clone("htoy");

	for (int itoy=0; itoy<ntoyExps; itoy++) {

	  htoy->Reset();

	  // decide how many events
	  int thisevt = trand3.Poisson(nevt_fitRange);

	  // fill histo with this # events randomly taken from total histo (in the requested range)
	  int NEV(0);
	  while (NEV<thisevt) {
	    float x = useRaw ? total->GetRandom() : totalFn->GetRandom() ;
	    if ( x > fitRangeLow && x < fitRangeUp ) {
	      htoy->Fill( x );
	      NEV++;
	    }
	  }

	  // for (int i=0; i<thisevt; i++) {
	  //   if ( useRaw ) htoy->Fill( total->GetRandom() );
	  //   else          htoy->Fill( totalFn->GetRandom() );
	  // }

	  func->SetParameter(0, thisevt);
	  func->SetParLimits(1, -1.0, 1.0);
	  func->SetParameter(1, -0.50);
	  func->SetParError(1, 0.1);

	  // then fit the toy histo in the requested range
	  TString opt="SLQ";
	  TFitResultPtr r = htoy->Fit("myfit",opt,"",fitRangeLow,fitRangeUp);

	  // get the fit parameters
	  float fitpol    = r->Parameter(1); 
	  float fitpolErr = r->ParError(1);

	  // save some of the histos for plotting
	  if ( itoy<25 ) {
	    TString dd = "toy"; dd+=itoy;
	    hclone[itoy]=(TH1F*) htoy->Clone(dd);
	  }

	  hStatus [icheat][iset][idec]->Fill(r->Status());
	  hPol   [icheat][iset][idec]->Fill(fitpol);
	  hPolErr[icheat][iset][idec]->Fill(fitpolErr);
	  
	  hPolPull[icheat][iset][idec]->Fill( (fitpol -  inputPol[icheat][iset][idec])/fitpolErr );
	  hFitProb[icheat][iset][idec]->Fill(r->Prob());

	}

	cc->Clear();
	cc->Divide(5,5);
	for (int i=0; i<25; i++) {
	  cc->cd(i+1);
	  if (hclone[i]) hclone[i]->Draw("e;hist");
	}
	cc->Print(plotfile);
	for (int i=0; i<25; i++) {
	  if (hclone[i]) delete hclone[i];
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

	if ( meanFitError[icheat][iset][idec] > maxmeanfiterror[iset] ) maxmeanfiterror[iset]=meanFitError[icheat][iset][idec];

	cc->cd(5); 

	hPolPull[icheat][iset][idec]->Fit("gaus","l;q");
	hPolPull[icheat][iset][idec]->Draw();

	pullMean    [icheat][iset][idec] = hPolPull[icheat][iset][idec]->GetFunction("gaus")->GetParameter(1);
	pullMeanErr [icheat][iset][idec] = hPolPull[icheat][iset][idec]->GetFunction("gaus")->GetParError(1);
	pullWidth   [icheat][iset][idec] = hPolPull[icheat][iset][idec]->GetFunction("gaus")->GetParameter(2);
	pullWidthErr[icheat][iset][idec] = hPolPull[icheat][iset][idec]->GetFunction("gaus")->GetParError(2);

	cc->Print(plotfile);

      }
    }

    for (int idec=0; idec<2; idec++) {
      cc->Clear();
      cc->Divide(2,2);
      for (int i=0; i<npols; i++) { // 4 pols, *2 detectors
	cc->cd(i+1);
	TH1F* h1 = hPolErr[icheat][i][idec];   // large
	//	TH1F* h2 = hPolErr[icheat][i+4][idec]; // small
	TH1F* h2 = hPolErr[icheat][i+npols][idec]; // small

	float mmx = std::max( h1->GetMean()+5*h1->GetRMS(),  h2->GetMean()+5*h2->GetRMS() );
	float mmn = std::min( h1->GetMean()-5*h1->GetRMS(),  h2->GetMean()-5*h2->GetRMS() );

	h1->GetXaxis()->SetRangeUser( mmn, mmx );

	h1->Draw();
	h2->Draw("same");
      }
      cc->Print(plotfile);
    }

  } // cheat



  std::ofstream fitresultsfile;
  fitresultsfile.open(outdir+"fitresults.txt");

  for (int idec=0; idec<2; idec++) {
    fitresultsfile << decStrs[idec] << endl;

    for (int iverb=0; iverb<2; iverb++) {
      for (int iset=0; iset<npols; iset++) {
	if ( iverb==0 ) fitresultsfile << " set " << iset << " : ";
	if ( iverb==0 ) fitresultsfile << decStrs[idec] << " ";
	for (int icheat=0; icheat<ncheat; icheat++) {
	  if ( iverb==0 ) fitresultsfile << cheatStrs[icheat] << " ( ";
	  int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
	  int ibest = meanFitError[icheat][iset+npols][idec] < meanFitError[icheat][iset][idec] ? 1 : 0 ;
	  for (int ils=0; ils<ILS; ils++) {
	    TString col = "black";
	    if ( cheatStrs[icheat]== "noCheatSEL" && ils==ibest) {
	      col=  ils==0 ? "blue" : "red";
	    }
	    TString hh; hh.Form(" & \\textcolor{%5s}{%3.1f}", col.Data(), 100*meanFitError[icheat][iset+ils*npols][idec] );
	    fitresultsfile << hh;
	  }
	  if ( iverb==0 ) fitresultsfile << " ) ";
	}
	
	for (int ils=0; ils<2; ils++) {
	  fitresultsfile << 100./ sqrt( pow(1./meanFitError[2][iset+ils*npols][0],2) + pow(1./meanFitError[2][iset+ils*npols][1],2) ) << " ";
	}
      
	fitresultsfile << " \\\\ " << endl;
      }
    }


    fitresultsfile << "\\multicolumn{9}{|c|}{INPUT POL \\%} \\\\" << endl;
    for (int iset=0; iset<npols; iset++) { // pol set
      for (int icheat=0; icheat<ncheat; icheat++) { // cheating
	int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
	for (int ils=0; ils<ILS; ils++) { // L/S
	  TString hh; hh.Form(" & $ %5.1f $", 100*inputPol[icheat][iset+ils*npols][idec] );
	  fitresultsfile << hh;
	}
      }
      fitresultsfile << " \\\\ " << endl;
    }

    fitresultsfile << "\\multicolumn{9}{|c|}{NEVENTS} \\\\" << endl;
    for (int iset=0; iset<npols; iset++) { // pol set
      for (int icheat=0; icheat<ncheat; icheat++) { // cheating
	int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
	for (int ils=0; ils<ILS; ils++) { // L/S
	  TString hh; hh.Form(" & $ %7.0f $", nevents[icheat][iset+ils*npols][idec] );
	  fitresultsfile << hh;
	}
      }
      fitresultsfile << " \\\\ " << endl;
    }
    
    fitresultsfile << "\\multicolumn{9}{|c|}{FITTED POL \\%} \\\\" << endl;
    for (int iset=0; iset<npols; iset++) {
      for (int icheat=0; icheat<ncheat; icheat++) {
	int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
	for (int ils=0; ils<ILS; ils++) {
	  TString hh; hh.Form(" & $ %5.1f \\pm  %5.1f $", 100*meanFitValue[icheat][iset+ils*npols][idec], 100*meanFitValueErr[icheat][iset+ils*npols][idec] );
	  fitresultsfile << hh;
	}
      }
      fitresultsfile << " \\\\ " << endl;
    }

    fitresultsfile << "\\multicolumn{9}{|c|}{MEAN FITTED POL Error \\%} \\\\" << endl;
    for (int iset=0; iset<npols; iset++) {
      for (int icheat=0; icheat<ncheat; icheat++) {
	int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
	for (int ils=0; ils<ILS; ils++) {
	  TString hh; hh.Form(" & $ %4.2f $", 100*meanFitError[icheat][iset+ils*npols][idec] );
	  fitresultsfile << hh;
	}
      }
      fitresultsfile << " \\\\ " << endl;
    }

    fitresultsfile << "\\multicolumn{9}{|c|}{PULL MEAN}\\\\" << endl;
    for (int iset=0; iset<npols; iset++) {
      for (int icheat=0; icheat<ncheat; icheat++) {
	int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
	for (int ils=0; ils<ILS; ils++) {
	  TString hh; hh.Form(" & $ %5.2f \\pm  %5.2f $", pullMean    [icheat][iset+ils*npols][idec], pullMeanErr [icheat][iset+ils*npols][idec] );
	  fitresultsfile << hh;
	}
      }
      fitresultsfile << " \\\\ " << endl;
    }

    fitresultsfile << "\\multicolumn{9}{|c|}{PULL WIDTH}\\\\" << endl;
    for (int iset=0; iset<npols; iset++) {
      for (int icheat=0; icheat<ncheat; icheat++) {
	int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
	for (int ils=0; ils<ILS; ils++) {
	  TString hh; hh.Form(" & $ %5.2f \\pm %5.2f $", pullWidth    [icheat][iset+ils*npols][idec], pullWidthErr [icheat][iset+ils*npols][idec] );
	  fitresultsfile << hh;
	}
      }
      fitresultsfile << " \\\\ " << endl;
    }
  }

  fitresultsfile.close();


  // draw plots

  int idec=0;
  float pi_prec[npols][8];
  for (int iset=0; iset<npols; iset++) {
    int jj=0;
    for (int icheat=0; icheat<ncheat; icheat++) {
      int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
      for (int ils=0; ils<ILS; ils++) {
	pi_prec[iset][jj++]=100*meanFitError[icheat][iset+ils*npols][idec];
      }
    }
  }

  idec=1;
  float rho_prec[npols][10];
  for (int iset=0; iset<npols; iset++) {
    int jj=0;
    for (int icheat=0; icheat<ncheat; icheat++) {
      int ILS = getNLS(cheatStrs[icheat], decStrs[idec] );
      for (int ils=0; ils<ILS; ils++) {
	rho_prec[iset][jj++]=100*meanFitError[icheat][iset+ils*npols][idec];
      }
    }
  }

  gStyle->SetMarkerSize(1.5);
  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetPadRightMargin(0.23);

  const int npP=5;
  const int npR=6;

  float pi_prec_L[npols][npP];
  float pi_prec_S[npols][npP];

  float rho_prec_L[npols][npR];
  float rho_prec_S[npols][npR];

  TGraph* gr_pi_L[npols];
  TGraph* gr_pi_S[npols];

  TGraph* gr_rho_L[npols];
  TGraph* gr_rho_S[npols];

  float indexP[npP] = {1,2,3,5,6};
  float indexR[npR] = {1,2,3,4,5,6};

  TString axLabs[npR] = {"MC-optimal", "MC-approx (no Sel/no Bg)", "MC-approx (Sel/no Bg)", "cheat ECAL (Sel/no Bg)", "no cheat (Sel/no Bg)", "no cheat (Sel/Bg)"};

  for (int i=0; i<npols; i++) {

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

  float mm(0);
  for (int j=0; j<npols; j++) {
    if ( maxmeanfiterror[j] > mm ) mm = maxmeanfiterror[j] ;
  }


  TH2F* hh = new TH2F("prec","",npR,0.5,npR+0.5,10,0, int(mm*100*1.5)+1. );
  hh->GetYaxis()->SetTitle("polarisation precision [%]");

  for (int i=1; i<=npR; i++) {
    hh->GetXaxis()->SetBinLabel(i, axLabs[i-1]);
  }


  cc = new TCanvas("c2","c2",400,400);
  //  cc->Print(indir+"drawpexpresults.pdf[");

  TText ttf;
  ttf.SetTextFont(62);
  ttf.SetTextSize(0.06);



  TLegend* tl=0;

  for (int i=0; i<npols; i++) {
    cc->Clear();

    hh->Draw();

    TString ttx = pollab[i];

    float ymax(2.);
    if ( pollab[i]=="eLpL" || pollab[i]=="eRpR" ) ymax=4.;

    hh->GetYaxis()->SetRangeUser(0, ymax );

    // if      (i==0 ) {hh->GetYaxis()->SetRangeUser(0, maxmeanfiterror*100 ); ttx="eLpR";}
    // else if (i==1 ) {hh->GetYaxis()->SetRangeUser(0, maxmeanfiterror*100. ); ttx="eRpL";}
    // else if (i==2 ) {hh->GetYaxis()->SetRangeUser(0, maxmeanfiterror*100. ); ttx="eLpL";}
    // else if (i==3 ) {hh->GetYaxis()->SetRangeUser(0, maxmeanfiterror*100. ); ttx="eRpR";}

    gr_rho_L[i]->Draw("pl");
    gr_rho_S[i]->Draw("pl");
    gr_pi_L[i]->Draw("pl");
    gr_pi_S[i]->Draw("pl");

    tt.DrawTextNDC(0.6, 0.25, ttx);

    if (!tl) {
      tl=new TLegend(0.2, 0.6, 0.5, 0.85);
      tl->AddEntry( gr_pi_L[i] , "#tau #rightarrow #pi IDR-L", "pl" );
      tl->AddEntry( gr_pi_S[i] , "#tau #rightarrow #pi IDR-S", "pl" );
      tl->AddEntry( gr_rho_L[i], "#tau #rightarrow #rho IDR-L", "pl" );
      tl->AddEntry( gr_rho_S[i], "#tau #rightarrow #rho IDR-S", "pl" );
    }

    tl->Draw();

    ttf.DrawTextNDC(0.38, 0.95, "ILD preliminary");


    TString dd="drawpexpresults"; dd+=i;
    cc->Print(outdir+dd+".eps");
    cc->Print(outdir+dd+".C");

    cc->Print(plotfile);

  }

  cc->Print(plotfile+"]");

  

}
