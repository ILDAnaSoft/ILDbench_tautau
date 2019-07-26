void effPlotByCut() {

  TFile* finput[4];

  
  finput[0] = new TFile("tauFind_ALLEVT_ILD_l5_o1_v02_2f_Z_leptonic_I250106_eL.pR_oneHadSS.root","read");
  finput[1] = new TFile("tauFind_ALLEVT_ILD_s5_o1_v02_2f_Z_leptonic_I250106_eL.pR_oneHadSS.root","read");
  finput[2] = new TFile("tauFind_ALLEVT_ILD_l5_o1_v02_2f_Z_leptonic_I250108_eR.pL_oneHadSS.root","read");
  finput[3] = new TFile("tauFind_ALLEVT_ILD_s5_o1_v02_2f_Z_leptonic_I250108_eR.pL_oneHadSS.root","read");

  TH1F* hPiRho[4][2][9];

  int icol[9]={1,2,3,4,5,6,7,8,9};

  TString label[9] = {"all", "presel", "2ndSeedEn", "outsideEnergy", "acoLin", "acoPlan", "ISR", "lepton", "clusterShape"};
  TString pname = "effPlotByCut.pdf";
  TCanvas* cc = new TCanvas();
  cc->Print( pname+"[" );
  TLegend* tl = new TLegend(0.45, 0.1, 0.65, 0.4);

  TString pollab[2]={"eLpR(100%)","eRpL(100%)"};
  TText tt;

  for ( int i=0; i<4; i++) {

    TFile* fin = finput[i];

    hPiRho[i][0][0] = (TH1F*) fin->Get("sample0_pi_mcPolarApprox");
    hPiRho[i][0][1] = (TH1F*) fin->Get("sample0_selA_pi_mcPolarApprox");
    hPiRho[i][0][2] = (TH1F*) fin->Get("sample0_selB_pi_mcPolarApprox");
    hPiRho[i][0][3] = (TH1F*) fin->Get("sample0_selC_pi_mcPolarApprox");
    hPiRho[i][0][4] = (TH1F*) fin->Get("sample0_selD_pi_mcPolarApprox");
    hPiRho[i][0][5] = (TH1F*) fin->Get("sample0_selE_pi_mcPolarApprox");
    hPiRho[i][0][6] = (TH1F*) fin->Get("sample0_selF_pi_mcPolarApprox");
    hPiRho[i][0][7] = (TH1F*) fin->Get("sample0_selG_pi_mcPolarApprox");
    hPiRho[i][0][8] = (TH1F*) fin->Get("sample0_selH_pi_mcPolarApprox");

    hPiRho[i][1][0] = (TH1F*) fin->Get("sample0_rho_mcPolarApprox");
    hPiRho[i][1][1] = (TH1F*) fin->Get("sample0_selA_rho_mcPolarApprox");
    hPiRho[i][1][2] = (TH1F*) fin->Get("sample0_selB_rho_mcPolarApprox");
    hPiRho[i][1][3] = (TH1F*) fin->Get("sample0_selC_rho_mcPolarApprox");
    hPiRho[i][1][4] = (TH1F*) fin->Get("sample0_selD_rho_mcPolarApprox");
    hPiRho[i][1][5] = (TH1F*) fin->Get("sample0_selE_rho_mcPolarApprox");
    hPiRho[i][1][6] = (TH1F*) fin->Get("sample0_selF_rho_mcPolarApprox");
    hPiRho[i][1][7] = (TH1F*) fin->Get("sample0_selG_rho_mcPolarApprox");
    hPiRho[i][1][8] = (TH1F*) fin->Get("sample0_selH_rho_mcPolarApprox");


    for (int idec=0; idec<2; idec++) {
      for (int isel=0; isel<9; isel++) {
	hPiRho[i][idec][isel]->Rebin(5);
	hPiRho[i][idec][isel]->SetLineColor(icol[isel]);
	hPiRho[i][idec][isel]->SetMarkerColor(icol[isel]);
	if ( i==0 && idec==0 ) {
	  tl->AddEntry(  hPiRho[i][idec][isel], label[isel], "l" );
	}
      }
    }

  }


  for (int ieff=0 ; ieff<2; ieff++ ) {

    if ( ieff==1 ) {
      for (int ifile=0; ifile<4; ifile++) {
	for (int idec=0; idec<2; idec++) {
	  for (int isel=1; isel<9; isel++) {
	    hPiRho[ifile][idec][isel]->Sumw2();
	    hPiRho[ifile][idec][isel]->Divide( hPiRho[ifile][idec][isel],  hPiRho[ifile][idec][0], 1, 1, "B" );
	  }
	}
      }
    }

    int iselmin= ieff==0 ? 0 : 1;

    for (int ipol=0; ipol<2; ipol++) {
      for (int idec=0; idec<2; idec++) {
	cc->Clear();
	for (int isel=iselmin; isel<9; isel++) {
	  for ( int idet=0; idet<2; idet++) {
	    TString opt = "";
	    if ( isel>iselmin || idet>0 ) opt+="same";
	    if ( ieff==0 ) {
	      if ( idet>0 ) opt+=";ple";
	    } else {
	      if ( idet==0 ) opt+=";hist";
	    }

	    hPiRho[2*ipol+idet][idec][isel]->GetYaxis()->SetTitle("efficiency");

	    if (idec==0 ) {
	      hPiRho[2*ipol+idet][idec][isel]->GetXaxis()->SetTitle("MC approx polarimeter (#pi)");
	    } else {
	      hPiRho[2*ipol+idet][idec][isel]->GetXaxis()->SetTitle("MC approx polarimeter (#rho)");
	    }

	    if ( ieff==1 ) {
	      hPiRho[2*ipol+idet][idec][isel]->GetYaxis()->SetRangeUser(0.5, 1.0);
	    }

	    hPiRho[2*ipol+idet][idec][isel]->Draw(opt);
	  }
	}
	tl->Draw();
	tt.DrawTextNDC( 0.6, 0.93, pollab[ipol]);

	tt.DrawTextNDC( 0.25, 0.25, "points: IDR-S");
	tt.DrawTextNDC( 0.25, 0.2, "line: IDR-L");
	cc->Print( pname);
      }
    }
    
  }


  cc->Print( pname+"]" );


}
