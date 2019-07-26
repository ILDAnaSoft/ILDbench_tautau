#include "TFile.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <vector>
#include <iostream>

using std::cout;
using std::endl;

TFile* fin[2][2];

int color_size[2]={4,2};
int ltype_size[2]={1,2};
int ptype_size[2]={21,22};

TString ILDlab="ILD preliminary";

TString figdir;

TText ttf;
TText ttt;
TLatex ttltx;

void init() {

  if ( 0 ) {

    fin[0][0]=new TFile("tauFind_ALLEVT_ILD_l5_o1_v02_2f_Z_leptonic_I250106_eL.pR.root","read");
    fin[0][1]=new TFile("tauFind_ALLEVT_ILD_l5_o1_v02_2f_Z_leptonic_I250108_eR.pL.root","read");
    fin[1][0]=new TFile("tauFind_ALLEVT_ILD_s5_o1_v02_2f_Z_leptonic_I250106_eL.pR.root","read");
    fin[1][1]=new TFile("tauFind_ALLEVT_ILD_s5_o1_v02_2f_Z_leptonic_I250108_eR.pL.root","read");

    figdir="figs_ALLEVT/";

  } else if (0) {

    fin[0][0]=new TFile("tauFind_ILD_l5_o1_v02_2f_Z_leptonic_I250106_eL.pR.root","read");
    fin[0][1]=new TFile("tauFind_ILD_l5_o1_v02_2f_Z_leptonic_I250108_eR.pL.root","read");
    fin[1][0]=new TFile("tauFind_ILD_s5_o1_v02_2f_Z_leptonic_I250106_eL.pR.root","read");
    fin[1][1]=new TFile("tauFind_ILD_s5_o1_v02_2f_Z_leptonic_I250108_eR.pL.root","read");

    figdir="figs/";

  } else if ( 0 ) {

    fin[0][0]=new TFile("tauFind_ILD_l5_o1_v02_2f_Z_leptonic_I250106_eL.pR_oneHad.root","read");
    fin[0][1]=new TFile("tauFind_ILD_l5_o1_v02_2f_Z_leptonic_I250108_eR.pL_oneHad.root","read");
    fin[1][0]=new TFile("tauFind_ILD_s5_o1_v02_2f_Z_leptonic_I250106_eL.pR_oneHad.root","read");
    fin[1][1]=new TFile("tauFind_ILD_s5_o1_v02_2f_Z_leptonic_I250108_eR.pL_oneHad.root","read");

    figdir="figs_onehad/";

  } else if ( 1 ) {

    fin[0][0]=new TFile("data_ALLEVToneHad/tauFind_ALLEVToneHad_ILD_l5_o1_v02_2f_Z_leptonic_I250106_eL.pR.root","read");
    fin[0][1]=new TFile("data_ALLEVToneHad/tauFind_ALLEVToneHad_ILD_l5_o1_v02_2f_Z_leptonic_I250108_eR.pL.root","read");
    fin[1][0]=new TFile("data_ALLEVToneHad/tauFind_ALLEVToneHad_ILD_s5_o1_v02_2f_Z_leptonic_I250106_eL.pR.root","read");
    fin[1][1]=new TFile("data_ALLEVToneHad/tauFind_ALLEVToneHad_ILD_s5_o1_v02_2f_Z_leptonic_I250108_eR.pL.root","read");

    figdir="data_ALLEVToneHad/figs_ALLEVToneHad/";

  }

  ttltx.SetTextSize(0.06);

  ttf.SetTextFont(62);
  ttf.SetTextSize(0.06);

  ttt.SetTextSize(0.06);

}

void mcplots() {

  int idet=0; // LS doesnt matter for mc plts

  gStyle->SetOptStat(0);


  TString pln=figdir+"allMCPlots.pdf";

  TCanvas* cc = new TCanvas("cc","cc",400,400);
  cc->Print(pln+"[");

  TString pols[2]={"eLpR", "eRpL"};
  TString pollabs[2]={"eLpR(100%)", "eRpL(100%)"};

  //  TText tt;

  for (int l=0; l<2; l++) {
    TString type=l==0?"Exact":"Approx";

    for (int k=0; k<2; k++) {
      TString dec=k==0?"pi":"rho";
      for (int j=0; j<2; j++) {
	TH1F* h3[3];
	h3[0] = (TH1F*) fin[idet][j]->Get("sample0_"+dec+"_mcPolar"+type);
	h3[1] = (TH1F*) fin[idet][j]->Get("sample0_"+dec+"_mcPolar"+type+"_helNeg");
	h3[2] = (TH1F*) fin[idet][j]->Get("sample0_"+dec+"_mcPolar"+type+"_helPos");

	TLegend tl(0.2, 0.75, 0.6, 0.9);
	tl.AddEntry(h3[1], "#minus ve helicity #tau","l");
	tl.AddEntry(h3[2], "#plus ve helicity #tau","l");

	cc->Clear();
	h3[0]->SetLineColor(1);
	h3[1]->SetLineColor(6);
	h3[2]->SetLineColor(8);

	TString ttyy=type;
	if ( ttyy == "Exact" ) ttyy = "Optimal";

	ttyy.ToLower();

	cout << type << " " << ttyy << endl;

	h3[0]->GetXaxis()->SetTitle("#" + dec + " " + ttyy + " polarimeter");
	h3[0]->SetMaximum(  1.45*h3[0]->GetMaximum() );
	h3[0]->SetTitle("");

	h3[0]->Draw("hist");
	h3[1]->Draw("same; hist");
	h3[2]->Draw("same; hist");

	tl.Draw();

	ttt.DrawTextNDC(0.65, 0.94, pollabs[j]);
	ttt.DrawTextNDC(0.2, 0.94, "MC level");

	cc->Print(figdir+"mcplot_"+pols[j]+"_"+dec+"_mcPolar"+type+".eps");
	cc->Print(figdir+"mcplot_"+pols[j]+"_"+dec+"_mcPolar"+type+".C");
	cc->Print(pln);
      }
    }
  }

  std::vector <TString> hnames1d;
  hnames1d.push_back("mcttmass");
  hnames1d.push_back("sample0_mc_tauMinus_costh");

  std::map < TString, TString > axlabs;
  axlabs["mcttmass"]="#tau-#tau mass [GeV/c^{2}]";
  axlabs["sample0_mc_tauMinus_costh"]="cos #theta_{#tau-}";

  TLegend* tl = new TLegend(0.35, 0.7, 0.75, 0.9);

  TH1F* h[2];
  for ( int i=0; i<hnames1d.size(); i++) {
    for (int j=0; j<2; j++) {
      h[j] = (TH1F*) fin[idet][j]->Get(hnames1d[i]);
      h[j]->SetLineColor( 6 + 2*j );

      if ( i == 1 ) h[j]->Rebin(5);
      h[j]->Scale(1./h[j]->GetEntries());
      h[j]->SetMinimum(0);
    }
    if ( i==0 ) {
      tl->AddEntry(h[0], pollabs[0], "l");
      tl->AddEntry(h[1], pollabs[1], "l");
    }
    cc->Clear();

    h[0]->GetXaxis()->SetTitle(axlabs[hnames1d[i]]);

    h[0]->SetTitle("");
    h[0]->Draw("hist");
    
    for (int j=1; j<2; j++) {
      h[j]->Draw("same;hist");
    }
    tl->Draw();
    ttt.DrawTextNDC(0.2, 0.94, "MC level");

    cc->Print(figdir+"mcplot_"+hnames1d[i]+".eps");
    cc->Print(figdir+"mcplot_"+hnames1d[i]+".C");
    cc->Print(pln);
  }

  

  cc->Print(pln+"]");

  //for (int j=0; j<2; j++) {
  //  fin[j]->Close();
  //}

}


void getDecayModeEff() {
  // get the tau decay mode rec eff matrix
//  TFile* fin[2];
//  fin[0]=new TFile("tauFind_ILD_l5_o1_v02_2f_Z_leptonic_I250106_eL.pR.root","read");
//  fin[1]=new TFile("tauFind_ILD_s5_o1_v02_2f_Z_leptonic_I250106_eL.pR.root","read");


  TCanvas* cc = new TCanvas("cc","cc",400,400);
  cc->Clear();
  cc->Divide(2,2);


  float geffNum[2][3];
  float geffDen[2][3];
  float gpurDen[2][3];

  cout << std::setprecision(2) << std::setw(8) << std::fixed;

  for (int idet=0; idet<2; idet++) {// large/small
    for (int k=0; k<3; k++) {
      geffNum[idet][k]=0;
      geffDen[idet][k]=0;
      gpurDen[idet][k]=0;
    }

    float effMatr[2][4][4];
    float totJet[2][4];


    for (int j=0; j<2; j++) {
    
      cout << "-----------" << endl;
      cout << "***** " << fin[idet][j]->GetName() << endl;

      cc->cd(2*idet + j + 1);
      TH1F* h1 = (TH1F*) fin[idet][j]->Get("sample0_SEL_mcdec");
      TH1F* h2 = (TH1F*) fin[idet][j]->Get("sample0_SEL_rec_pi_mcdec");
      TH1F* h3 = (TH1F*) fin[idet][j]->Get("sample0_SEL_rec_rho_mcdec");
      TH1F* h4 = (TH1F*) fin[idet][j]->Get("sample0_SEL_rec_a1p_mcdec");
      h1->SetLineColor(1);
      h2->SetLineColor(2);
      h3->SetLineColor(4);
      h4->SetLineColor(6);
      h1->Draw();
      h2->Draw("same");
      h3->Draw("same");
      h4->Draw("same");
      
      effMatr[j][0][0]=h1->GetBinContent(2); // all selected cones, matched to MC pi
      effMatr[j][0][1]=h1->GetBinContent(3); // all selected cones, matched to MC rho
      effMatr[j][0][2]=h1->GetBinContent(4); // all selected cones, matched to MC a1p
      effMatr[j][0][3]=h1->GetEntries()-effMatr[j][0][0]-effMatr[j][0][1]-effMatr[j][0][2]; // all selected cones, mathced to neither MC pi nor rho
      totJet[j][0]    =h1->GetEntries(); // all selected cones
      
      effMatr[j][1][0]=h2->GetBinContent(2); // all cones ID as rho, matched to MC pi
      effMatr[j][1][1]=h2->GetBinContent(3); // all cones ID as rho, matched to MC rho
      effMatr[j][1][2]=h2->GetBinContent(4); // all cones ID as rho, matched to MC a1p
      effMatr[j][1][3]=h2->GetEntries()-effMatr[j][1][0]-effMatr[j][1][1]-effMatr[j][1][2];
      totJet[j][1]    =h2->GetEntries();// all cones ID as rho
      
      effMatr[j][2][0]=h3->GetBinContent(2);// all cones ID as pi, matched to MC pi
      effMatr[j][2][1]=h3->GetBinContent(3);// all cones ID as pi, matched to MC rho
      effMatr[j][2][2]=h3->GetBinContent(4);// all cones ID as pi, matched to MC a1p
      effMatr[j][2][3]=h3->GetEntries()-effMatr[j][2][0]-effMatr[j][2][1]-effMatr[j][2][2];
      totJet[j][2]    =h3->GetEntries();// all cones ID as pi
      
      effMatr[j][3][0]=h4->GetBinContent(2);// all cones ID as a1p, matched to MC pi
      effMatr[j][3][1]=h4->GetBinContent(3);// all cones ID as a1p, matched to MC rho
      effMatr[j][3][2]=h4->GetBinContent(4);// all cones ID as a1p, matched to MC a1p
      effMatr[j][3][3]=h4->GetEntries()-effMatr[j][3][0]-effMatr[j][3][1]-effMatr[j][3][2];
      totJet[j][3]    =h4->GetEntries();// all cones ID as pi

      cout << "-----------   matched to      MCpi    MCrho   MCa1p       MCother      MCall" << endl;
      for ( int ii=0; ii<4; ii++) {
	if      ( ii == 0 ) cout << "ALL CONES                 ";
	else if ( ii == 1 ) cout << "ALL CONES RECO as PI      ";
	else if ( ii == 2 ) cout << "ALL CONES RECO as RHO     ";
	else if ( ii == 3 ) cout << "ALL CONES RECO as A1P     ";
	for ( int jj=0; jj<4; jj++) {
	  cout << std::setw(8) << effMatr[j][ii][jj] << " " ;
	}
	cout << " : " << std::setw(8) << totJet[j][ii] << endl;
      }
      
      cout << "-----------" << endl;
      for ( int ii=1; ii<4; ii++) {
	if      ( ii == 0 ) cout << "ALL SELECTED JETS   ";
	else if ( ii == 1 ) cout << "SELECTED AS PI      ";
	else if ( ii == 2 ) cout << "SELECTED AS RHO     ";
	else if ( ii == 3 ) cout << "SELECTED AS A1P     ";
	for ( int jj=0; jj<4; jj++) {
	  float eff = effMatr[j][ii][jj]/effMatr[j][0][jj];
	  cout << 100*eff << " +- " <<  100*sqrt(eff*(1.-eff)/effMatr[j][0][jj]) <<  "    ";
	}
	if ( ii>0 ) {
	  float pur = effMatr[j][ii][ii-1]/totJet[j][ii];
	  cout << " purity " << 100*pur << " +- " << 100*sqrt ( pur*(1.-pur)/totJet[j][ii] );
	}
	cout << endl;
      }

      for (int k=0; k<3; k++) {
	geffNum[idet][k]+=effMatr[j][k+1][k];
	geffDen[idet][k]+=effMatr[j][0][k];
	gpurDen[idet][k]+=totJet[j][k+1];
      }

    }


    cout << "-----------" << endl;
    cout << " --------- AVE OVER SAMPLES ----------------------------" << endl;
    cout << "-----------" << endl;
    cout << "-----------   matched to      MCpi    MCrho   MCa1p       MCother      MCall" << endl;
    for ( int ii=0; ii<4; ii++) {
      if      ( ii == 0 ) cout << "ALL CONES                 ";
      else if ( ii == 1 ) cout << "ALL CONES RECO as PI      ";
      else if ( ii == 2 ) cout << "ALL CONES RECO as RHO     ";
      else if ( ii == 3 ) cout << "ALL CONES RECO as A1P     ";
      for ( int jj=0; jj<4; jj++) {
	cout << std::setw(8) << effMatr[0][ii][jj]+effMatr[1][ii][jj] << " " ;
      }
      cout << " : " << std::setw(8) << totJet[0][ii]+totJet[1][ii] << endl;
    }
    
    cout << "-----------" << endl;
    for ( int ii=1; ii<4; ii++) {
      if      ( ii == 0 ) cout << "ALL SELECTED JETS  & $ ";
      else if ( ii == 1 ) cout << "SELECTED AS PI     & $ ";
      else if ( ii == 2 ) cout << "SELECTED AS RHO    & $ ";
      else if ( ii == 3 ) cout << "SELECTED AS A1P    & $ ";
      for ( int jj=0; jj<4; jj++) {
	float eff = (effMatr[0][ii][jj]+effMatr[1][ii][jj])/(effMatr[0][0][jj]+effMatr[1][0][jj]);
	cout << std::setw(8) << 100*eff << " \\pm " <<  std::setw(8) << 100*sqrt(eff*(1.-eff)/(effMatr[0][0][jj]+effMatr[1][0][jj])) <<  " $ & $    ";
      }
      if ( ii>0 ) {
	float pur = (effMatr[0][ii][ii-1]+effMatr[1][ii][ii-1])/(totJet[0][ii]+totJet[1][ii]);
	// cout << " purity " << 100*pur << " +- " << 100*sqrt ( pur*(1.-pur)/(totJet[0][ii]+totJet[1][ii]) );
	cout << std::setw(8) << 100*pur << " \\pm " << std::setw(8) << 100*sqrt ( pur*(1.-pur)/(totJet[0][ii]+totJet[1][ii]) );
      }
      cout << " $ \\\\ " << endl;
    }




  }
  cc->Print(figdir+"mcplot_decay.eps");
  cc->Print(figdir+"mcplot_decay.C");

  float geff[2][3];
  float gdeff[2][3];
  float gpur[2][3];
  float gdpur[2][3];
  TGraphErrors* gre[2];

  float dec[3]={0.5,1.5,2.5};
  float zero[3]={0,0,0};
  float geffpur[2][3];
  float gdeffpur[2][3];
  TGraphErrors* gre_effPur[2];

  for (int idet=0; idet<2; idet++) {// large/small
    for (int k=0; k<3; k++) {

      float feff = geffNum[idet][k]/geffDen[idet][k];
      float fpur = geffNum[idet][k]/gpurDen[idet][k];

      float deff = sqrt ( feff*(1.-feff) / geffDen[idet][k] );
      float dpur = sqrt ( fpur*(1.-fpur) / gpurDen[idet][k] );

      geff[idet][k] = feff; // geffNum[idet][k]/geffDen[idet][k];
      gpur[idet][k] = fpur; // geffNum[idet][k]/gpurDen[idet][k];
      
      gdeff[idet][k] = deff; // sqrt( geff[idet][k]*(1.-geff[idet][k] ) / geffDen[idet][k] );
      gdpur[idet][k] = dpur; // sqrt( gpur[idet][k]*(1.-gpur[idet][k] ) / gpurDen[idet][k] );

      geffpur[idet][k] = feff*fpur;
      gdeffpur[idet][k] =  feff*fpur * sqrt( pow( deff/feff, 2 ) + pow( dpur/fpur, 2 ) );


    }

    gre[idet] = new TGraphErrors(3, &(geff[idet][0]), &(gpur[idet][0]), &(gdeff[idet][0]), &(gdpur[idet][0]) );

    gre_effPur[idet] = new TGraphErrors(3, dec, &(geffpur[idet][0]), zero, &(gdeffpur[idet][0]) );
    
  }

  TLegend* tl = new TLegend(0.23, 0.75, 0.5, 0.88);
  tl->SetFillStyle(0);
  tl->AddEntry( gre[0], "IDR-L", "l" );
  tl->AddEntry( gre[1], "IDR-S", "l" );


  cc->Clear();
  TH2F* hs = new TH2F("effPur","",10,0.45, 1., 10, 0.45, 1.);
  hs->GetXaxis()->SetTitle("efficiency");
  hs->GetYaxis()->SetTitle("purity");
  hs->Draw();
  for (int idet=0; idet<2; idet++) {// large/small
    gre[idet]->SetMarkerColor( color_size[idet] );
    gre[idet]->SetLineColor( color_size[idet] );
    gre[idet]->SetLineStyle( ltype_size[idet] );
    gre[idet]->Draw("same;pl");
  }
  tl->Draw();

  ttltx.DrawLatexNDC(0.82, 0.64, "#pi");
  ttltx.DrawLatexNDC(0.58, 0.78, "#rho");
  ttltx.DrawLatexNDC(0.33, 0.3 , "a_{1}");

  ttf.DrawTextNDC(0.5, 0.95, ILDlab);
  ttltx.DrawLatexNDC(0.55, 0.35, "e^{+}e^{#minus} #rightarrow #tau^{+}#tau^{#minus}");
  ttltx.DrawLatexNDC(0.55, 0.25, "#sqrt{s} = 500 GeV");


  cc->Print(figdir+"decmodeEffPur.eps");
  cc->Print(figdir+"decmodeEffPur.C");



  cc->Clear();
  TH2F* hseq = new TH2F("effPur","",3,0., 3., 10, 0.2, 0.85);
  hseq->GetXaxis()->SetBinLabel(1, "#pi");
  hseq->GetXaxis()->SetBinLabel(2, "#rho");
  hseq->GetXaxis()->SetBinLabel(3, "a_{1}");
  hseq->GetYaxis()->SetTitle("efficiency #times purity");

  hseq->GetXaxis()->SetLabelSize(0.13);

  hseq->Draw();
  for (int idet=0; idet<2; idet++) {// large/small
    gre_effPur[idet]->SetMarkerColor( color_size[idet] );
    gre_effPur[idet]->SetLineColor( color_size[idet] );
    gre_effPur[idet]->SetLineStyle( ltype_size[idet] );
    gre_effPur[idet]->Draw("same;pl");
  }
  tl->SetX1NDC(0.65);
  tl->SetX2NDC(0.85);
  tl->Draw();

  ttf.DrawTextNDC(0.5, 0.95, ILDlab);
  ttltx.DrawLatexNDC(0.22, 0.35, "e^{+}e^{#minus} #rightarrow #tau^{+}#tau^{#minus}");
  ttltx.DrawLatexNDC(0.22, 0.25, "#sqrt{s} = 500 GeV");

  cc->Print(figdir+"decmodeEffTimesPur.eps");
  cc->Print(figdir+"decmodeEffTimesPur.C");




  tl->SetX1NDC(0.65);
  tl->SetX2NDC(0.85);


  TH1F* hdpol[2][2];
  TH1D* hpold[2][2];
  for (int idet=0; idet<2; idet++) {// large/small
    for (int j=0; j<2; j++) {
      cc->cd(2*idet+j+1);
      hdpol[idet][j] = (TH1F*) fin[idet][j]->Get("sample0_SEL_recrho_mcrho_dpol");
      hdpol[idet][j]->SetLineColor( color_size[idet] );
      hdpol[idet][j]->SetLineStyle( ltype_size[idet] );
      hdpol[idet][j]->Rebin(2);
    }
    hdpol[idet][0]->Add( hdpol[idet][1] );
  }

  cc->Clear();
  //  hdpol[0][0]->SetMaximum(1.2*hdpol[0][0]->GetMaximum());
  hdpol[0][0]->Draw();
  hdpol[0][0]->SetTitle("");
  hdpol[0][0]->GetXaxis()->SetRangeUser(-0.5, 0.5);
  hdpol[0][0]->GetXaxis()->SetTitle("#rho polarimeter (rec - MC)");
  hdpol[1][0]->Draw("same");

  tl->Draw();
  ttf.DrawTextNDC(0.5, 0.95, ILDlab);

  cc->Print(figdir+"dPolRho.eps");
  cc->Print(figdir+"dPolRho.C");


  for (int idet=0; idet<2; idet++) {// large/small
    for (int j=0; j<2; j++) {
      cc->cd(2*idet+j+1);
      hdpol[idet][j] = (TH1F*) fin[idet][j]->Get("sample0_SEL_recpi_mcpi_dpol");
      hdpol[idet][j]->SetLineColor( color_size[idet] );
      hdpol[idet][j]->SetLineStyle( ltype_size[idet] );
      //      hdpol[idet][j]->Rebin(2);
    }
    hdpol[idet][0]->Add( hdpol[idet][1] );
  }

  cc->Clear();
  hdpol[0][0]->SetMaximum(1.25 * hdpol[0][0]->GetMaximum() );
  hdpol[0][0]->Draw();
  hdpol[0][0]->SetTitle("");
  hdpol[0][0]->GetXaxis()->SetRangeUser(-0.1, 0.1);
  hdpol[0][0]->GetXaxis()->SetTitle("#pi polarimeter (rec - MC)");
  hdpol[1][0]->Draw("same");

  tl->Draw();
  ttf.DrawTextNDC(0.5, 0.95, ILDlab);

  cc->Print(figdir+"dPolPi.eps");
  cc->Print(figdir+"dPolPi.C");


  for (int idet=0; idet<2; idet++) {// large/small
    for (int j=0; j<2; j++) {
      cc->cd(2*idet+j+1);
      TString hhg="xx"; hhg+=idet; hhg+=j;
      hpold[idet][j] = ( (TH2F*) fin[idet][j]->Get("sample0_SEL_recrho_mcrho_pol") )->ProjectionY(hhg);
      hpold[idet][j]->SetLineColor( color_size[idet] );
      hpold[idet][j]->SetLineStyle( ltype_size[idet] );
      hpold[idet][j]->SetMarkerColor( color_size[idet] );
      hpold[idet][j]->SetMarkerStyle( ptype_size[idet] );
      hpold[idet][j]->Rebin(5);
    }
  }

  cc->Clear();
  hpold[0][0]->SetMaximum(1.25*hpold[0][0]->GetMaximum());
  hpold[0][0]->Draw();
  hpold[0][0]->SetTitle("");
  hpold[0][0]->GetXaxis()->SetTitle("#rho polarimeter (rec)");
  hpold[1][0]->Draw("same");
  hpold[0][1]->Draw("same;p;e");
  hpold[1][1]->Draw("same;p;e");

  tl->Draw();
  ttf.DrawTextNDC(0.5, 0.95, ILDlab);

  cc->Print(figdir+"PolRho.eps");
  cc->Print(figdir+"PolRho.C");

  for (int idet=0; idet<2; idet++) {// large/small
    for (int j=0; j<2; j++) {
      cc->cd(2*idet+j+1);
      TString hhg="xx"; hhg+=idet; hhg+=j;
      hpold[idet][j] = ( (TH2F*) fin[idet][j]->Get("sample0_SEL_recpi_mcpi_pol") )->ProjectionY(hhg);
      hpold[idet][j]->SetLineColor( color_size[idet] );
      hpold[idet][j]->SetLineStyle( ltype_size[idet] );
      hpold[idet][j]->SetMarkerColor( color_size[idet] );
      hpold[idet][j]->SetMarkerStyle( ptype_size[idet] );
      //      hpold[idet][j]->Rebin(10);
      hpold[idet][j]->Rebin(5);
    }
  }

  cc->Clear();

  cout << hpold[0][0]->GetMaximum() << endl;

  hpold[0][0]->SetMaximum(1.25*hpold[1][0]->GetMaximum());
  hpold[0][0]->Draw();
  hpold[0][0]->SetTitle("");
  hpold[0][0]->GetXaxis()->SetTitle("#pi polarimeter (rec)");
  hpold[1][0]->Draw("same");
  hpold[0][1]->Draw("same;p;e");
  hpold[1][1]->Draw("same;p;e");

  tl->Draw();
  ttf.DrawTextNDC(0.5, 0.95, ILDlab);

  cc->Print(figdir+"PolPi.eps");
  cc->Print(figdir+"PolPi.C");



  for (int idet=0; idet<2; idet++) {// large/small
    for (int j=0; j<2; j++) {
      cc->cd(2*idet+j+1);
      TString hhg="xx"; hhg+=idet; hhg+=j;
      hpold[idet][j] = ( (TH2F*) fin[idet][j]->Get("sample0_SEL_recrho_mcrho_pol") )->ProjectionX(hhg);
      hpold[idet][j]->SetLineColor( color_size[idet] );
      hpold[idet][j]->SetLineStyle( ltype_size[idet] );
      hpold[idet][j]->SetMarkerColor( color_size[idet] );
      hpold[idet][j]->SetMarkerStyle( ptype_size[idet] );
      hpold[idet][j]->Rebin(5);
    }
  }

  cc->Clear();
  hpold[0][0]->SetMaximum(1.25*hpold[0][0]->GetMaximum());
  hpold[0][0]->Draw();
  hpold[0][0]->SetTitle("");
  hpold[0][0]->GetXaxis()->SetTitle("#rho polarimeter (MC)");
  hpold[1][0]->Draw("same");
  hpold[0][1]->Draw("same;p;e");
  hpold[1][1]->Draw("same;p;e");

  tl->Draw();
  ttf.DrawTextNDC(0.5, 0.95, ILDlab);

  cc->Print(figdir+"PolRhoMC.eps");
  cc->Print(figdir+"PolRhoMC.C");

  for (int idet=0; idet<2; idet++) {// large/small
    for (int j=0; j<2; j++) {
      cc->cd(2*idet+j+1);
      TString hhg="xx"; hhg+=idet; hhg+=j;
      hpold[idet][j] = ( (TH2F*) fin[idet][j]->Get("sample0_SEL_recpi_mcpi_pol") )->ProjectionX(hhg);
      hpold[idet][j]->SetLineColor( color_size[idet] );
      hpold[idet][j]->SetLineStyle( ltype_size[idet] );
      hpold[idet][j]->SetMarkerColor( color_size[idet] );
      hpold[idet][j]->SetMarkerStyle( ptype_size[idet] );
      //      hpold[idet][j]->Rebin(10);
      hpold[idet][j]->Rebin(5);
    }
  }

  cc->Clear();
  //  hpold[0][0]->SetMaximum(1.2*hpold[0][0]->GetMaximum());

  cout << hpold[0][0]->GetMaximum() << endl;

  hpold[0][0]->SetMaximum( hpold[0][1]->GetMaximum() * 1.45 );

  cout << hpold[0][0]->GetMaximum() << endl;

  hpold[0][0]->Draw();

  cout << hpold[0][0]->GetMaximum() << endl;

  hpold[0][0]->SetTitle("");
  hpold[0][0]->GetXaxis()->SetTitle("#pi polarimeter (mc)");
  hpold[1][0]->Draw("same");
  hpold[0][1]->Draw("same;p;e");
  hpold[1][1]->Draw("same;p;e");

  tl->Draw();
  ttf.DrawTextNDC(0.5, 0.95, ILDlab);

  cc->Print(figdir+"PolPiMC.eps");
  cc->Print(figdir+"PolPiMC.C");


//  for (int j=0; j<2; j++) {
//    fin[j]->Close();
//  }


}


void simpleLargeSmallPlots() {

  gStyle->SetOptStat(0);

//  TFile* fin[2];
//  fin[0]=new TFile("tauFind_ALLEVT_ILD_l5_o1_v02_2f_Z_leptonic_I250106_eL.pR.root","read");
//  fin[1]=new TFile("tauFind_ALLEVT_ILD_s5_o1_v02_2f_Z_leptonic_I250106_eL.pR.root","read");

  std::vector <TString> hnames1d;
  hnames1d.push_back("sample0__PI_cone_ngammapfo");
  hnames1d.push_back("sample0__RHO_cone_ngammapfo");
  hnames1d.push_back("sample0__A1P_cone_ngammapfo");
  hnames1d.push_back("sample0__RHO_cone_visMass");
  hnames1d.push_back("sample0__A1P_cone_visMass");
  hnames1d.push_back("sample0__RHO_SELcone_visEnergyDiff");
  hnames1d.push_back("sample0__A1P_SELcone_visEnergyDiff");
  //  hnames1d.push_back("sample0__RHO_SELcone_neutralvisEnergyDiff");

  std::map < TString, TString > xaxtitle;
  xaxtitle["sample0__PI_cone_ngammapfo"]="# reconstructed photons";
  xaxtitle["sample0__RHO_cone_ngammapfo"]="# reconstructed photons";
  xaxtitle["sample0__A1P_cone_ngammapfo"]="# reconstructed photons";
  xaxtitle["sample0__RHO_cone_visMass"]="visible mass [GeV/c^{2}]";
  xaxtitle["sample0__A1P_cone_visMass"]="visible mass [GeV/c^{2}]";
  xaxtitle["sample0__RHO_SELcone_visEnergyDiff"]="en. diff. (rec-true) [GeV]";
  xaxtitle["sample0__A1P_SELcone_visEnergyDiff"]="en. diff. (rec-true) [GeV]";
  //  xaxtitle["sample0__RHO_SELcone_neutralvisEnergyDiff"]="neutral en. diff. (rec-true) [GeV]";

  std::map < TString, int > linlog;
  linlog["sample0__PI_cone_ngammapfo"]               =0;
  linlog["sample0__RHO_cone_ngammapfo"]              =0;
  linlog["sample0__A1P_cone_ngammapfo"]              =0;
  linlog["sample0__RHO_cone_visMass"]                =0;
  linlog["sample0__A1P_cone_visMass"]                =0;
  linlog["sample0__RHO_SELcone_visEnergyDiff"]       =0;
  linlog["sample0__A1P_SELcone_visEnergyDiff"]       =0;
  //  linlog["sample0__RHO_SELcone_neutralvisEnergyDiff"]=0;

  std::map < TString, int > rebin;
  rebin["sample0__PI_cone_ngammapfo"]               =1;
  rebin["sample0__RHO_cone_ngammapfo"]              =1;
  rebin["sample0__A1P_cone_ngammapfo"]              =1;
  rebin["sample0__RHO_cone_visMass"]                =1;
  rebin["sample0__A1P_cone_visMass"]                =1;
  rebin["sample0__RHO_SELcone_visEnergyDiff"]       =2;
  rebin["sample0__A1P_SELcone_visEnergyDiff"]       =4;
  //  rebin["sample0__RHO_SELcone_neutralvisEnergyDiff"]=2;

  std::map < TString, float > hmaxX;
  hmaxX["sample0__PI_cone_gammaen"]=20.;
  hmaxX["sample0__PI_cone_nhaden"]=40;
  hmaxX["sample0__RHO_SELconeTRIM_neutralvisMass"]=1.;

  TCanvas* cc = new TCanvas("cc","cc",400,400);
  cc->Print(figdir+"allSimpleLargeSmallPlots.pdf[");

  TLegend* tl = new TLegend( 0.6, 0.7, 0.9, 0.9);
  TH1F* h[2];
  for ( int i=0; i<hnames1d.size(); i++) {

    int ilog = 0;
    if ( linlog.find( hnames1d[i] ) != linlog.end() ) {
      ilog=linlog[ hnames1d[i] ];
    }

    float hmax(0);
    float maxx=-1000;
    if ( hmaxX.find( hnames1d[i] ) != hmaxX.end() ) maxx = hmaxX[hnames1d[i]];

    int ireb = rebin.find(  hnames1d[i] ) != rebin.end() ? rebin[hnames1d[i]]: 1;


    for (int j=0; j<2; j++) {
      h[j] = (TH1F*) fin[j][0]->Get(hnames1d[i]);
      
      if (ireb>1) h[j] ->Rebin(ireb);

      if ( maxx>-999 ) {
	h[j]->GetXaxis()->SetRangeUser( h[j]->GetXaxis()->GetBinLowEdge(1), maxx );
      }

      h[j]->SetLineColor( color_size[j] );
      h[j]->SetLineStyle( ltype_size[j] );

      h[j]->Scale(1./h[j]->GetEntries());
      if ( ilog==1 ) h[j]->SetMinimum(0.001);
      
      if ( h[j]->GetMaximum()>hmax) hmax= h[j]->GetMaximum();

    }
    if ( i==0 ) {
      tl->AddEntry(h[0], "IDR-L", "l");
      tl->AddEntry(h[1], "IDR-S", "l");
    }
    cc->Clear();


    cc->SetLogy(ilog);
    if ( xaxtitle.find(hnames1d[i])!=xaxtitle.end() ) {
      h[0]->SetTitle("");
      h[0]->GetXaxis()->SetTitle(xaxtitle[hnames1d[i]]);
    }
    h[0]->SetMaximum(hmax*1.3);
    h[0]->Draw("hist");
    for (int j=1; j<2; j++) {
      h[j]->Draw("same;hist");
    }
    tl->Draw();

    if ( hnames1d[i].Contains("_RHO_") ) {
      ttltx.DrawLatexNDC(0.23, 0.82, "#tau #rightarrow #rho");
    } else if ( hnames1d[i].Contains("_PI_") ) {
      ttltx.DrawLatexNDC(0.23, 0.82, "#tau #rightarrow #pi");
    } else if ( hnames1d[i].Contains("_A1P_") ) {
      ttltx.DrawLatexNDC(0.23, 0.82 , "#tau #rightarrow a_{1}");
    }
    ttf.DrawTextNDC(0.5, 0.95, ILDlab);


    cc->Print(figdir+"compareplot_"+hnames1d[i]+".eps");
    cc->Print(figdir+"compareplot_"+hnames1d[i]+".C");
    cc->Print(figdir+"allSimpleLargeSmallPlots.pdf");
  }

  cc->SetLogy(0);

  TH2F* hj = (TH2F*) fin[0][0]->Get("neutralHadronPFO_mainMCcontrib_pdgEn");
  hj->RebinX(2);
  hj->GetYaxis()->SetRangeUser(1,5);
  hj->SetTitle("");
  hj->GetXaxis()->SetTitle("PFO energy [GeV]");
  hj->Draw("box");
  ttf.DrawTextNDC(0.5, 0.95, ILDlab);


  cc->Print(figdir+"neutralHadPFOs.eps");
  cc->Print(figdir+"neutralHadPFOs.C");
  cc->Print(figdir+"allSimpleLargeSmallPlots.pdf");

  TH1F* jj = (TH1F*) fin[0][0]->Get("rhoDecaySinglePhoClus_reason");
  jj->SetTitle("");
  jj->Draw();
  ttf.DrawTextNDC(0.5, 0.95, ILDlab);
  cc->Print(figdir+"onlyOneClusReason.eps");
  cc->Print(figdir+"onlyOneClusReason.C");
  cc->Print(figdir+"allSimpleLargeSmallPlots.pdf");

  TH2F* jj2 = (TH2F*) fin[0][0]->Get("neutronClus_eOnP_beforeAfter");
  jj2->SetTitle("");
  jj2->GetXaxis()->SetTitle("E/p [orig]");
  jj2->GetYaxis()->SetTitle("E/p [merged]");
  jj2->Draw("box");
  ttf.DrawTextNDC(0.5, 0.95, ILDlab);
  cc->Print(figdir+"EonPmergingNHAD.eps");
  cc->Print(figdir+"EonPmergingNHAD.C");
  cc->Print(figdir+"allSimpleLargeSmallPlots.pdf");


  cc->Clear();
  cc->SetLogy();
  TH1D* l1 = jj2->ProjectionX();
  TH1D* l2 = jj2->ProjectionY();

  TLegend ddl(0.4, 0.7, 0.8, 0.9);
  ddl.AddEntry(l1, "original", "l");
  ddl.AddEntry(l2, "with neutral cluster", "l");

  l1->GetXaxis()->SetTitle("E/p");
  
  l1->Draw();
  l2->SetLineColor(2);
  l2->Draw("same");
  ttf.DrawTextNDC(0.5, 0.95, ILDlab);

  ddl.Draw();

  cc->Print(figdir+"EonPmergingNHAD_1d.eps");
  cc->Print(figdir+"EonPmergingNHAD_1d.C");
  cc->Print(figdir+"allSimpleLargeSmallPlots.pdf");
  cc->SetLogy(0);


  TH1F* jj3 = (TH1F*) fin[0][0]->Get("neutronClus_eOnPdiff_beforeAfter");
  cc->Clear();
  jj3->SetTitle("");
  jj3->GetXaxis()->SetTitle("|E/p|_{orig} - |E/p|_{merge}");
  //  jj3->GetYaxis()->SetTitle("events");
  jj3->Draw();
  ttf.DrawTextNDC(0.5, 0.95, ILDlab);
  
  cc->Print(figdir+"EonPmergingNHAD_DIFF.eps");
  cc->Print(figdir+"EonPmergingNHAD_DIFF.C");
  cc->Print(figdir+"allSimpleLargeSmallPlots.pdf");




  TH2F* jjk[2][4];
  cc->Clear();
  jjk[0][0] = (TH2F*) fin[0][0]->Get("rhoDecaySinglePhoClus_mergedGammaClusterEval1");
  jjk[0][1] = (TH2F*) fin[0][0]->Get("rhoDecaySinglePhoClus_mergedGammaClusterEval2");
  jjk[0][2] = (TH2F*) fin[0][0]->Get("rhoDecaySinglePhoClus_mergedGammaClusterEvalRatio");

  jjk[1][0] = (TH2F*) fin[0][0]->Get("rhoDecaySinglePhoClus_singleGammaClusterEval1");
  jjk[1][1] = (TH2F*) fin[0][0]->Get("rhoDecaySinglePhoClus_singleGammaClusterEval2");
  jjk[1][2] = (TH2F*) fin[0][0]->Get("rhoDecaySinglePhoClus_singleGammaClusterEvalRatio");

  for (int j=0; j<3; j++) {
    jjk[0][j]->SetLineColor(2);

    jjk[0][j]->SetTitle("");
    jjk[0][j]->GetXaxis()->SetTitle("cluster energy [GeV]");

    jjk[0][j]->Rebin2D(2,2);
    jjk[1][j]->Rebin2D(2,2);

    jjk[0][j]->Draw("box");
    jjk[1][j]->Draw("box;same");

    ttt.SetTextColor(2);
    ttt.DrawTextNDC(0.65, 0.2, "merged");

    ttt.SetTextColor(1);
    ttt.DrawTextNDC(0.4, 0.2, "single");


    if ( j==0 ) {
      jjk[0][j]->GetYaxis()->SetTitle("1st cluster e-val [mm]");
      cc->Print(figdir+"eval1.eps");
      cc->Print(figdir+"eval1.C");
    } else if ( j==1 ) {
      jjk[0][j]->GetYaxis()->SetTitle("2nd cluster e-val [mm]");
      cc->Print(figdir+"eval2.eps");
      cc->Print(figdir+"eval2.C");
    } else if ( j==2 ) {
      jjk[0][j]->GetYaxis()->SetTitle("1st/2nd cluster e-val ratio");
      cc->Print(figdir+"evalRatio.eps");
      cc->Print(figdir+"evalRatio.C");
    }



    cc->Print(figdir+"allSimpleLargeSmallPlots.pdf");
  }


  for ( int j=0; j<4; j++) {
    for (int i=0; i<2; i++) {
      TString dec = i==0? "RHO" : "A1P";
      int col=i==0?1:2;
      TString nn; nn+=j+1;
      jjk[i][j] = (TH2F*) fin[0][0]->Get("sample0__"+dec+"_SELcone_vis_neutral_Mass_"+nn+"g");
      jjk[i][j]->Rebin2D(4,4);
      jjk[i][j]->SetLineColor(col);
    }

    jjk[0][j]->Draw("box");
    jjk[1][j]->Draw("box; same");
    ttf.DrawTextNDC(0.5, 0.95, ILDlab);

    cc->Print(figdir+"allSimpleLargeSmallPlots.pdf");

  }

  TH1F* ggg[3][3];
  ggg[0][0] = (TH1F*) fin[0][0]->Get("sample0__PI_seed_logclusterWidth1");
  ggg[0][1] = (TH1F*) fin[0][0]->Get("sample1__E_seed_logclusterWidth1");
  ggg[0][2] = (TH1F*) fin[0][0]->Get("sample1__M_seed_logclusterWidth1");
  ggg[1][0] = (TH1F*) fin[0][0]->Get("sample0__PI_seed_logclusterWidth2");
  ggg[1][1] = (TH1F*) fin[0][0]->Get("sample1__E_seed_logclusterWidth2");
  ggg[1][2] = (TH1F*) fin[0][0]->Get("sample1__M_seed_logclusterWidth2");
  ggg[2][0] = (TH1F*) fin[0][0]->Get("sample0__PI_seed_logclusterLength");
  ggg[2][1] = (TH1F*) fin[0][0]->Get("sample1__E_seed_logclusterLength");
  ggg[2][2] = (TH1F*) fin[0][0]->Get("sample1__M_seed_logclusterLength");

  for (int i=0; i<3; i++) {
    cc->Clear();
    float hmm(0);
    for (int j=0; j<3; j++) {
      ggg[i][j]->SetLineColor(j+1);
      if ( ggg[i][j]->GetMaximum()>hmm ) hmm=ggg[i][j]->GetMaximum();
    }

    for (int j=0; j<3; j++) {
      ggg[i][j]->SetMaximum( hmm*1.2);
      TString opt=j==0?"":"same";
      ggg[i][j]->Draw(opt);
    }
    cc->Print(figdir+"allSimpleLargeSmallPlots.pdf");
  }

  TH1F* hEff[2][2];

  for (int itdec=0; itdec<2; itdec++) {

    TString dd = itdec==0 ? "pi" : "rho";

    for ( int idet=0; idet<2; idet++) {


      TH1F* heff[2];

      for (int ipol=0; ipol<2; ipol++) {

	if ( ipol==0 ) {
	  heff[0]=(TH1F*) fin[idet][ipol] -> Get("sample0_SEL_"+dd+"_mcPolarExact_helNeg");
	  heff[0]->Add( (TH1F*) fin[idet][ipol] -> Get("sample0_SEL_"+dd+"_mcPolarExact_helPos") ); 

	  heff[1]=(TH1F*) fin[idet][ipol] -> Get("sample0_"+dd+"_mcPolarExact");
	} else {
	  heff[0]->Add((TH1F*) fin[idet][ipol] -> Get("sample0_SEL_"+dd+"_mcPolarExact_helNeg") );
	  heff[0]->Add( (TH1F*) fin[idet][ipol] -> Get("sample0_SEL_"+dd+"_mcPolarExact_helPos") ); 

	  heff[1]->Add((TH1F*) fin[idet][ipol] -> Get("sample0_"+dd+"_mcPolarExact"));
	}

      }

      heff[0]->Rebin(5);
      heff[1]->Rebin(5);

      //cc->Clear();
      //heff[0]->SetLineColor(4);
      //heff[1]->Draw();
      //heff[0]->Draw("same");
      //
      //cc->Print(figdir+"allSimpleLargeSmallPlots.pdf");

      hEff[idet][itdec] = (TH1F*) heff[0]->Clone("eff");
      hEff[idet][itdec]->Sumw2();
      hEff[idet][itdec]->SetLineColor( color_size[idet] );
      hEff[idet][itdec]->SetMarkerColor( color_size[idet] );
      hEff[idet][itdec]->SetLineStyle( ltype_size[idet] );

      hEff[idet][itdec]->Divide(heff[0],heff[1], 1, 1, "B");

      // hEff[idet][itdec]->Draw("e");
      //      hEff[idet][itdec]->Fit("pol1");
      //cc->Print(figdir+"allSimpleLargeSmallPlots.pdf");

    }

  }

  cc->Clear();
  cc->Divide(1,2);
  for (int itdec=0; itdec<2; itdec++) {
    cc->cd(itdec+1);

    hEff[0][itdec]->SetTitle("");
    hEff[0][itdec]->GetYaxis()->SetRangeUser(0.25, 0.75);

    hEff[0][itdec]->GetXaxis()->SetTitle("exact polarimeter (MC)");
    hEff[0][itdec]->GetYaxis()->SetTitle("selection efficiency");

    hEff[0][itdec]->Draw("e");
    hEff[1][itdec]->Draw("same;e");

    ttltx.SetTextSize(0.12);
    if ( itdec==0 ) {
      ttltx.DrawLatexNDC(0.45, 0.42, "#tau #rightarrow #pi");
    } else if (  itdec==1 ) {
      ttltx.DrawLatexNDC(0.45, 0.42, "#tau #rightarrow #rho");
    }

    tl->SetX1NDC(0.7);
    tl->SetX2NDC(0.9);
    tl->SetY1NDC(0.2);
    tl->SetY2NDC(0.4);
    tl->Draw();

    ttf.DrawTextNDC(0.7, 0.95, ILDlab);

  }
  cc->Print(figdir+"selectionEffsPols.eps");
  cc->Print(figdir+"selectionEffsPols.C");

  cc->Print(figdir+"allSimpleLargeSmallPlots.pdf");


//  TH2F* hcl[2][2];
//  for ( int idet=0; idet<2; idet++) {
//    for (int ipol=0; ipol<2; ipol++) {
//
//      TH2F* h0=(TH2F*) fin[idet][ipol] -> Get("rhoDecaySinglePhoClus_mergedGammaClusterEvals");
//      TH2F* h1=(TH2F*) fin[idet][ipol] -> Get("rhoDecaySinglePhoClus_singleGammaClusterEvals");
//
//      cout <<  hcl[idet][0] << " " <<  hcl[idet][1] << " " << h0 << " " << h1 << endl;
//
//      if (ipol==0) {
//	hcl[idet][0]=h0;
//	hcl[idet][1]=h1;
//      } else {
//	hcl[idet][0]->Add(h0);
//	hcl[idet][1]->Add(h1);
//      }
//
//      hcl[idet][0]->SetLineColor(2);
//
//    }
//
//  }

//  cc->Clear();
//  cc->Divide(2,2);
//  int ic=1;
//  for ( int idet=0; idet<2; idet++) {
//
//    for (int i=0; i<2; i++) hcl[idet][i]->Scale(1./ hcl[idet][i]->GetEntries());
//
//    cc->cd(ic++);
//    hcl[idet][0]->ProjectionX()->Draw("hist");
//    hcl[idet][1]->ProjectionX()->Draw("same;hist");
//
//    cc->cd(ic++);
//    hcl[idet][0]->ProjectionY()->Draw("hist");
//    hcl[idet][1]->ProjectionY()->Draw("hist;same");
//
//
//  }
//  cc->Print(figdir+"allSimpleLargeSmallPlots.pdf");

  cc->Print(figdir+"allSimpleLargeSmallPlots.pdf]");

}
