#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
using std::cout; using std::endl;

void photonsInCones() {

  TFile* fin[2];
  fin[0] = new TFile("tauFind_ILD_l5_o1_v02_2f_Z_leptonic_I250106_eL.pR.root","read");
  fin[1] = new TFile("tauFind_ILD_s5_o1_v02_2f_Z_leptonic_I250106_eL.pR.root","read");

  TString hnames[4];
  hnames[0] = "sample0__PI_cone_ngammapfo";
  hnames[1] = "sample0__RHO_cone_ngammapfo";
  hnames[2] = "sample0__A1P_cone_ngammapfo";
  hnames[3] = "sample0__A3P_cone_ngammapfo";

  TCanvas* cc = new TCanvas();
  cc->Divide(2,2);

  TString allname="photonsInCones.pdf";
  cc->Print(allname+"[");

  for (int i=0; i<4; i++) {
    cc->cd(i+1);
    TH1F* h0 = (TH1F*) fin[0]->Get(hnames[i]);
    TH1F* h1 = (TH1F*) fin[1]->Get(hnames[i]);
    h0->SetLineColor(2);
    h1->SetLineColor(4);
    h1->Draw();
    h0->Draw("same");
  }

  cc->Print("nRecNgammaPfos.eps");
  cc->Print(allname);


  cc->Clear();
  cc->Divide(2,2);

  hnames[0] =  "sample0__PI_cone_gammaen";
  hnames[1] = "sample0__RHO_cone_gammaen";
  hnames[2] = "sample0__A1P_cone_gammaen";
  hnames[3] = "sample0__A3P_cone_gammaen";

  for (int i=0; i<4; i++) {
    cc->cd(i+1);
    TH1F* h0 = (TH1F*) fin[0]->Get(hnames[i]);
    TH1F* h1 = (TH1F*) fin[1]->Get(hnames[i]);
    h0->SetLineColor(2);
    h1->SetLineColor(4);
    h1->Draw();
    h0->Draw("same");
  }

  cc->Print("totGammaEn.eps");
  cc->Print(allname);



  hnames[0] =  "sample0__PI_cone_nnhadpfo";
  hnames[1] = "sample0__RHO_cone_nnhadpfo";
  hnames[2] = "sample0__A1P_cone_nnhadpfo";
  hnames[3] = "sample0__A3P_cone_nnhadpfo";

  cc->Clear();
  cc->Divide(2,2);

  for (int i=0; i<4; i++) {
    cc->cd(i+1);
    TH1F* h0 = (TH1F*) fin[0]->Get(hnames[i]);
    TH1F* h1 = (TH1F*) fin[1]->Get(hnames[i]);
    h0->SetLineColor(2);
    h1->SetLineColor(4);
    h1->SetMaximum(1.1*h1->GetMaximum());
    h1->Draw();
    h0->Draw("same");
  }

  cc->Print("nRecNnhadPfos.eps");
  cc->Print(allname);


  hnames[0] =  "sample0__PI_cone_nchpfo";
  hnames[1] = "sample0__RHO_cone_nchpfo";
  hnames[2] = "sample0__A1P_cone_nchpfo";
  hnames[3] = "sample0__A3P_cone_nchpfo";

  cc->Clear();
  cc->Divide(2,2);

  for (int i=0; i<4; i++) {
    cc->cd(i+1);
    TH1F* h0 = (TH1F*) fin[0]->Get(hnames[i]);
    TH1F* h1 = (TH1F*) fin[1]->Get(hnames[i]);
    h0->SetLineColor(2);
    h1->SetLineColor(4);
    h1->SetMaximum(1.1*h1->GetMaximum());
    h1->Draw();
    h0->Draw("same");
  }

  cc->Print("nChgPfos.eps");
  cc->Print(allname);


  hnames[0] =  "sample0__PI_cone_ntracks";
  hnames[1] = "sample0__RHO_cone_ntracks";
  hnames[2] = "sample0__A1P_cone_ntracks";
  hnames[3] = "sample0__A3P_cone_ntracks";

  cc->Clear();
  cc->Divide(2,2);

  for (int i=0; i<4; i++) {
    cc->cd(i+1);
    TH1F* h0 = (TH1F*) fin[0]->Get(hnames[i]);
    TH1F* h1 = (TH1F*) fin[1]->Get(hnames[i]);
    h0->SetLineColor(2);
    h1->SetLineColor(4);
    h1->SetMaximum(1.1*h1->GetMaximum());
    h1->Draw();
    h0->Draw("same");
  }

  cc->Print("nTracks.eps");
  cc->Print(allname);



  hnames[0] =  "sample0__PI_cone_trackNTPC";
  hnames[1] = "sample0__RHO_cone_trackNTPC";
  hnames[2] = "sample0__A1P_cone_trackNTPC";
  hnames[3] = "sample0__A3P_cone_trackNTPC";

  cc->Clear();
  cc->Divide(2,2);

  for (int i=0; i<4; i++) {
    cc->cd(i+1);
    TH1F* h0 = (TH1F*) fin[0]->Get(hnames[i]);
    TH1F* h1 = (TH1F*) fin[1]->Get(hnames[i]);
    h0->SetLineColor(2);
    h1->SetLineColor(4);
    h1->SetMaximum(1.1*h1->GetMaximum());
    h1->Draw();
    h0->Draw("same");
  }

  cc->Print("trackNTPC.eps");
  cc->Print(allname);


  cc->Clear();
  cc->Divide(2,2);

  hnames[0] = "sample0__RHO_SELcone_visMass";
  hnames[1] = "sample0__RHO_SELcone_neutralvisMass";
  hnames[2] = "sample0__RHO_SELconeTRIM_visMass";
  hnames[3] = "sample0__RHO_SELconeTRIM_neutralvisMass";

  for (int i=0; i<4; i++) {
    cc->cd(i+1);
    TH1F* h0 = (TH1F*) fin[0]->Get(hnames[i]);
    TH1F* h1 = (TH1F*) fin[1]->Get(hnames[i]);
    h0->SetLineColor(2);
    h1->SetLineColor(4);
    h1->Draw();
    h0->Draw("same");
  }

  cc->Print("tauJetMass.eps");
  cc->Print(allname);





  cc->Clear();
  cc->Divide(2,2);

  hnames[0] = "sample0__RHO_SELconeTRIM_visMass";
  hnames[1] = "sample0__RHO_SELconeTRIM_neutralvisMass";
  hnames[2] = "sample0__A1P_SELconeTRIM_visMass";
  hnames[3] = "sample0__A1P_SELconeTRIM_neutralvisMass";

  for (int i=0; i<4; i++) {
    cc->cd(i+1);
    TH1F* h0 = (TH1F*) fin[0]->Get(hnames[i]);
    TH1F* h1 = (TH1F*) fin[1]->Get(hnames[i]);
    h0->SetLineColor(2);
    h1->SetLineColor(4);
    h1->Draw();
    h0->Draw("same");
  }

  cc->Print("tauJetMass_a1_rho.eps");
  cc->Print(allname);

  cc->Clear();
  cc->Divide(2,2);

  hnames[0] = "sample0__PI_SELconeTRIM_visMass";
  hnames[1] = "sample0__PI_SELconeTRIM_neutralvisMass";
  hnames[2] = "sample0__RHO_SELconeTRIM_visMass";
  hnames[3] = "sample0__RHO_SELconeTRIM_neutralvisMass";

  for (int i=0; i<4; i++) {
    cc->cd(i+1);
    TH1F* h0 = (TH1F*) fin[0]->Get(hnames[i]);
    TH1F* h1 = (TH1F*) fin[1]->Get(hnames[i]);
    h0->SetLineColor(2);
    h1->SetLineColor(4);
    h1->Draw();
    h0->Draw("same");
  }

  cc->Print("tauJetMass_rho_pi.eps");
  cc->Print(allname);

  cc->Clear();
  cc->Divide(2,2);

  hnames[0] = "sample0__RHO_cone_npi0pfo";
  hnames[1] = "sample0__A1P_cone_npi0pfo";
  hnames[2] = "sample0__OTHER_cone_npi0pfo";
  hnames[3] = "sample0__PI_cone_npi0pfo";

  for (int i=0; i<4; i++) {
    cc->cd(i+1);
    TH1F* h0 = (TH1F*) fin[0]->Get(hnames[i]);
    TH1F* h1 = (TH1F*) fin[1]->Get(hnames[i]);
    h0->SetLineColor(2);
    h1->SetLineColor(4);
    h1->Draw();
    h0->Draw("same");
  }

  cc->Print("tauJetMass_npi0.eps");
  cc->Print(allname);

  TString ppp[3]={"PI","RHO","A1P"};

  cc->Clear();
  cc->Divide(2,2);

  for  (int ip=0; ip<3; ip++) {

    cc->cd(ip+1);

    hnames[0] = "sample0__"+ppp[ip]+"_mc_ngamma";
    hnames[1] = "sample0__"+ppp[ip]+"_mc_cone_ngamma";
    hnames[2] = "sample0__"+ppp[ip]+"_mcall_cone_ngamma";
    hnames[3] = "sample0__"+ppp[ip]+"_cone_ngammapfo";


    TH1F* h0 = (TH1F*) fin[0]->Get(hnames[0]);
    TH1F* h1 = (TH1F*) fin[0]->Get(hnames[1]);
    TH1F* h2 = (TH1F*) fin[0]->Get(hnames[2]);
    TH1F* h3 = (TH1F*) fin[0]->Get(hnames[3]);
    TH1F* h4 = (TH1F*) fin[1]->Get(hnames[3]);

    h3->SetLineColor(2);
    h4->SetLineColor(4);

    h0->SetLineColor(1);  h0->SetLineStyle(1);
    h1->SetLineColor(1);  h1->SetLineStyle(2);
    h2->SetLineColor(1);  h2->SetLineStyle(3);

    h0->Draw();
    h1->Draw("same");
    h2->Draw("same");
    h3->Draw("same");
    h4->Draw("same");

  }

  cc->Print("mcphotons_diffDec.eps");
  cc->Print(allname);



  TString pp[2]={"RHO","A1P"};

  for (int ip=0; ip<2; ip++) {
    cc->Clear();

    TH1F* h0 = (TH1F*) fin[0]->Get("sample0__"+pp[ip]+"_mc_vismass");
    TH1F* h1 = (TH1F*) fin[0]->Get("sample0__"+pp[ip]+"_SELcone_visMass");
    TH1F* h2 = (TH1F*) fin[0]->Get("sample0__"+pp[ip]+"_SELconeTRIM_visMass");
    TH1F* h3 = (TH1F*) fin[1]->Get("sample0__"+pp[ip]+"_SELcone_visMass");
    TH1F* h4 = (TH1F*) fin[1]->Get("sample0__"+pp[ip]+"_SELconeTRIM_visMass" );

    h0->SetLineColor(1);
    h0->Scale( 1.*h1->GetEntries()/h0->GetEntries() );

    h1->SetLineColor(2);
    h2->SetLineColor(2);  h2->SetLineStyle(2);

    h3->SetLineColor(4);
    h4->SetLineColor(4);  h4->SetLineStyle(2);


    cc->Clear();
    cc->Divide(2,2);

    cc->cd(1);

    h0->Draw("hist");
    h1->Draw("same");
    h2->Draw("same");

    cc->cd(2);

    h0->Draw("hist");
    h3->Draw("same");
    h4->Draw("same");

    cc->cd(3);

    h0->Draw("hist");
    h1->Draw("same");
    h3->Draw("same");

    cc->cd(4);

    h0->Draw("hist");
    h2->Draw("same");
    h4->Draw("same");


    cc->Print("vismass_"+pp[ip]+".eps");
    cc->Print(allname);

  }

  for (int ip=0; ip<2; ip++) {
    cc->Clear();

    TH1F* h0 = (TH1F*) fin[0]->Get("sample0__"+pp[ip]+"_mc_visneutralmass");
    TH1F* h1 = (TH1F*) fin[0]->Get("sample0__"+pp[ip]+"_SELcone_neutralvisMass");
    TH1F* h2 = (TH1F*) fin[0]->Get("sample0__"+pp[ip]+"_SELconeTRIM_neutralvisMass");
    TH1F* h3 = (TH1F*) fin[1]->Get("sample0__"+pp[ip]+"_SELcone_neutralvisMass");
    TH1F* h4 = (TH1F*) fin[1]->Get("sample0__"+pp[ip]+"_SELconeTRIM_neutralvisMass");


    h0->SetLineColor(1);
    
    h0->Scale( 1.*h1->GetEntries()/h0->GetEntries() );

    h1->SetLineColor(2);
    h2->SetLineColor(2);  h2->SetLineStyle(2);

    h3->SetLineColor(4);
    h4->SetLineColor(4);  h4->SetLineStyle(2);


    cc->Clear();
    cc->Divide(2,2);

    cc->cd(1);

    h0->Draw("hist");
    h1->Draw("same");
    h2->Draw("same");

    cc->cd(2);

    h0->Draw("hist");
    h3->Draw("same");
    h4->Draw("same");

    cc->cd(3);

    h0->Draw("hist");
    h1->Draw("same");
    h3->Draw("same");

    cc->cd(4);

    h0->Draw("hist");
    h2->Draw("same");
    h4->Draw("same");


    cc->Print("visneutralmass_"+pp[ip]+".eps");
    cc->Print(allname);

  }




  for (int ip=0; ip<2; ip++) {
    cc->Clear();

    TH1F* h1 = (TH1F*) fin[0]->Get("sample0__"+pp[ip]+"_SELcone_visEnergyDiff");
    TH1F* h2 = (TH1F*) fin[0]->Get("sample0__"+pp[ip]+"_SELconeTRIM_visEnergyDiff");
    TH1F* h3 = (TH1F*) fin[1]->Get("sample0__"+pp[ip]+"_SELcone_visEnergyDiff");
    TH1F* h4 = (TH1F*) fin[1]->Get("sample0__"+pp[ip]+"_SELconeTRIM_visEnergyDiff");

    h1->SetLineColor(2);
    h2->SetLineColor(2);  h2->SetLineStyle(2);

    h3->SetLineColor(4);
    h4->SetLineColor(4);  h4->SetLineStyle(2);

    cc->Clear();
    cc->Divide(2,2);

    cc->cd(1);

    h1->Draw("");
    h2->Draw("same");

    cc->cd(2);

    h3->Draw("");
    h4->Draw("same");

    cc->cd(3);

    h1->Draw("");
    h3->Draw("same");

    cc->cd(4);

    h2->Draw("");
    h4->Draw("same");


    cc->Print("visEnergyDiff_"+pp[ip]+".eps");
    cc->Print(allname);

  }


  for (int ip=0; ip<2; ip++) {
    cc->Clear();

    TH1F* h1 = (TH1F*) fin[0]->Get("sample0__"+pp[ip]+"_SELcone_neutralvisEnergyDiff");
    TH1F* h2 = (TH1F*) fin[0]->Get("sample0__"+pp[ip]+"_SELconeTRIM_neutralvisEnergyDiff");
    TH1F* h3 = (TH1F*) fin[1]->Get("sample0__"+pp[ip]+"_SELcone_neutralvisEnergyDiff");
    TH1F* h4 = (TH1F*) fin[1]->Get("sample0__"+pp[ip]+"_SELconeTRIM_neutralvisEnergyDiff");

    h1->SetLineColor(2);
    h2->SetLineColor(2);  h2->SetLineStyle(2);

    h3->SetLineColor(4);
    h4->SetLineColor(4);  h4->SetLineStyle(2);

    cc->Clear();
    cc->Divide(2,2);

    cc->cd(1);

    h1->Draw("");
    h2->Draw("same");

    cc->cd(2);

    h3->Draw("");
    h4->Draw("same");

    cc->cd(3);

    h1->Draw("");
    h3->Draw("same");

    cc->cd(4);

    h2->Draw("");
    h4->Draw("same");


    cc->Print("neutralvisEnergyDiff_"+pp[ip]+".eps");
    cc->Print(allname);

  }



  for (int ij=0; ij<6; ij++) {

    cc->Clear();
    cc->Divide(2,2);

    TString suff="";
    if (ij>0) {
      suff="_"; suff+=ij-1; suff+="g";
    }

    hnames[0] = "sample0__ALL_SELconeTRIM_vis_neutral_Mass" + suff;
    hnames[1] =  "sample0__PI_SELconeTRIM_vis_neutral_Mass" + suff;
    hnames[2] = "sample0__RHO_SELconeTRIM_vis_neutral_Mass" + suff;
    hnames[3] = "sample0__A1P_SELconeTRIM_vis_neutral_Mass" + suff;

    float hmax=0;
    for (int i=0; i<4; i++) {
      cc->cd(i+1)->SetLogz();
      TH2F* h0 = (TH2F*) fin[0]->Get(hnames[i]);
      if (i==0) hmax = h0->GetMaximum();
      else h0->SetMaximum(hmax);

      h0->SetMinimum(1);

      h0->GetXaxis()->SetRangeUser(0,2);
      h0->GetYaxis()->SetRangeUser(0,2);

      //    TH1F* h1 = (TH1F*) fin[1]->Get(hnames[i]);
      //h0->SetLineColor(2);
      //h1->SetLineColor(4);
      h0->Draw("zcol");
      //h0->Draw("same");
    }

    cc->Print("tauJetVisNeuTrimMass"+suff+".eps");
    cc->Print(allname);

  }

  for (int ij=0; ij<6; ij++) {

    cc->Clear();
    cc->Divide(2,2);

    TString suff="";
    if (ij>0) {
      suff="_"; suff+=ij-1; suff+="g";
    }

    hnames[0] = "sample0__ALL_SELcone_vis_neutral_Mass" + suff;
    hnames[1] = "sample0__PI_SELcone_vis_neutral_Mass"  + suff;
    hnames[2] = "sample0__RHO_SELcone_vis_neutral_Mass" + suff;
    hnames[3] = "sample0__A1P_SELcone_vis_neutral_Mass" + suff;

    float hmax=0;
    for (int i=0; i<4; i++) {
      cc->cd(i+1)->SetLogz();
      TH2F* h0 = (TH2F*) fin[0]->Get(hnames[i]);
      if (i==0) hmax = h0->GetMaximum();
      else h0->SetMaximum(hmax);

      h0->SetMinimum(1);

      h0->GetXaxis()->SetRangeUser(0,2);
      h0->GetYaxis()->SetRangeUser(0,2);

      //    TH1F* h1 = (TH1F*) fin[1]->Get(hnames[i]);
      //h0->SetLineColor(2);
      //h1->SetLineColor(4);
      h0->Draw("zcol");
      //h0->Draw("same");
    }

    cc->Print("tauJetVisNeuMass"+suff+".eps");
    cc->Print(allname);

  }


  cc->Clear();

  TH1F* h0 = (TH1F*) fin[0]->Get("sample0_pi_coneFit_polRecMcDiff");
  TH1F* h1 = (TH1F*) fin[0]->Get("sample0_pi_ipFit_polRecMcDiff");
  TH1F* h2 = (TH1F*) fin[0]->Get("sample0_rho_coneFit_polRecMcDiff");
  TH1F* h3 = (TH1F*) fin[0]->Get("sample0_rho_ipFit_polRecMcDiff");

  TH1F* h4 = (TH1F*) fin[1]->Get("sample0_pi_coneFit_polRecMcDiff");
  TH1F* h5 = (TH1F*) fin[1]->Get("sample0_pi_ipFit_polRecMcDiff");
  TH1F* h6 = (TH1F*) fin[1]->Get("sample0_rho_coneFit_polRecMcDiff");
  TH1F* h7 = (TH1F*) fin[1]->Get("sample0_rho_ipFit_polRecMcDiff");

  h0->SetLineColor(2);  h0->SetLineStyle(2);
  h1->SetLineColor(2);
  h2->SetLineColor(2);  h2->SetLineStyle(2);
  h3->SetLineColor(2);

  h4->SetLineColor(4);  h4->SetLineStyle(2);
  h5->SetLineColor(4);
  h6->SetLineColor(4);  h6->SetLineStyle(2);
  h7->SetLineColor(4);


  h0->Scale(0.5);
  h2->Scale(0.5);
  h4->Scale(0.5);
  h6->Scale(0.5);

  cc->Divide(3,2);
  cc->cd(1); h0->Draw("hist"); h4->Draw("hist;same");
  cc->cd(2); h1->Draw("hist"); h5->Draw("hist;same");
  cc->cd(3); h1->Draw("hist"); h0->Draw("hist;same");
  cc->cd(4); h2->Draw("hist"); h6->Draw("hist;same");
  cc->cd(5); h3->Draw("hist"); h7->Draw("hist;same");
  cc->cd(6); h3->Draw("hist"); h2->Draw("hist;same");

  cc->Print("solutions.eps");
  cc->Print(allname);



  cc->Print(allname+"]");


}
