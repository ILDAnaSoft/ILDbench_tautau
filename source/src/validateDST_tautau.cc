#include "validateDST_tautau.h"
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <fstream>
#include <algorithm>

#include <marlin/Global.h>
#include "lcio.h"

#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"

#include "TLorentzVector.h"
#include "TVector3.h"

#include "UTIL/LCRelationNavigator.h"

#include "tauUtils.h"

#include "TauSpinner/tau_reweight_lib.h"
#include "Tauola/Tauola.h"

using namespace lcio ;
using namespace marlin ;
using std::cout;
using std::endl;

validateDST_TauTauProcessor avalidateDST_TauTauProcessor ;


validateDST_TauTauProcessor::validateDST_TauTauProcessor() : Processor("validateDST_TauTauProcessor") {

  // modify processor description
  _description = "validateDST_TauTauProcessor makes histograms of basic properties for tau pair events" ;

  registerProcessorParameter("outputFilename",
                             "name of output file",
                             _outfile,
                             std::string( "validateTauTau.root") );

  registerProcessorParameter("dbdsample",
			     "is this a dbd sample?",
			     _isDBD,
			     int (0) );

  registerProcessorParameter("pythiaTauDecays",
			     "taus decayed by pythia? (no intermediate states in event record)",
			     _isPythiaTauDecay,
			     int (0) );

  registerProcessorParameter("mccolskimmed",
			     "skimmed mclist?",
			     _mcskim,
			     int (0) );

  _eventFitter = new tautau2fEventFitter();
  _eventFitter->setVerbose();




}


void validateDST_TauTauProcessor::init() {
  cout << "hello from validateDST_TauTauProcessor::init" << endl;

  iplot=0;
  // the output root file for histograms
  _fout = new TFile(_outfile.c_str(),"recreate");

  TString tauDecStrings[NTAUDEC]={"mu", "el", "a1_3p", "a1_1p", "rho", "hadW", "kchg", "pi"};
  TString pfoTypeStrings[NTYPE]={"mu", "el", "gam", "chad", "nhad", "kshort", "lambda"};

  _h_eventCounter = new TH1F("eventCounter","eventCounter",3,-0.5,2.5);

  _h_mcDecMode = new TH1F("mcDecMode","mcDecMode",NTAUDEC,-0.5, NTAUDEC-0.5);
  _h_mcDecMode2d = new TH2F("mcDecMode2d","mcDecMode2d",NTAUDEC,-0.5, NTAUDEC-0.5, NTAUDEC,-0.5, NTAUDEC-0.5);
  for (int i=0; i<NTAUDEC; i++) {
    _h_mcDecMode->GetXaxis()->SetBinLabel(i+1, tauDecStrings[i] );
    _h_mcDecMode2d->GetXaxis()->SetBinLabel(i+1, tauDecStrings[i] );
    _h_mcDecMode2d->GetYaxis()->SetBinLabel(i+1, tauDecStrings[i] );
  }

  _h_mcHas94 = new TH1F( "mcHas94",  "mcHas94", 4,-1,3 );

  _h_mctaunu = new TH1F( "mctaunu",  "mctaunu", 10,0,10 );

  _h_mcTauPlusVtxRZ  = new TH2F( "mcTauPlusVtxRZ", "mcTauPlusVtxRZ", 100, 0, 100, 100, -50, 50 );
  _h_mcTauPlusVtxXY  = new TH2F( "mcTauPlusVtxXY", "mcTauPlusVtxXY", 100, -50, 50, 100, -50, 50 );
  _h_mcTauMinusVtxRZ = new TH2F( "mcTauMinusVtxRZ", "mcTauMinusVtxRZ", 100, 0, 100, 100, -50, 50 );
  _h_mcTauMinusVtxXY = new TH2F( "mcTauMinusVtxXY", "mcTauMinusVtxXY", 100, -50, 50, 100, -50, 50 );

  _h_mcTauMinusDecl  = new TH1F( "mcTauMinusDecl",  "mcTauMinusDecl", 100,0,50);
  _h_mcTauPlusDecl   = new TH1F( "mcTauPlusDecl",  "mcTauPlusDecl", 100,0,50);

  _h_mcTauMinusLife  = new TH1F( "mcTauMinusLife",  "mcTauMinusLife", 250,0,5);
  _h_mcTauPlusLife   = new TH1F( "mcTauPlusLife",  "mcTauPlusLife", 250,0,5);

  _h_mcTauTauMass    = new TH1F( "mcTauTauMass", "mcTauTauMass", 200, 0, 510 );
  _h_mcMuMuMass      = new TH1F( "mcMuMuMass", "mcMuMuMass", 200, 0, 510 );

  _h_mctoten_mumu    = new TH1F("mctoten_mumu","mctoten_mumu",500,0,2000) ;
  _h_mctoten_tautau  = new TH1F("mctoten_tautau","mctoten_tautau",500,0,2000) ;

  _h_mcVis_aColin  = new TH2F( "mcVis_aColin",  "mcVis_aColin",  250, 0, 500, 400, 0, 3.14159);
  _h_mcVis_aCoplan = new TH2F( "mcVis_aCoplan", "mcVis_aCoplan", 250, 0, 500, 400, 0, 3.14159);
  _h_mcVis_invMass = new TH1F( "mcVis_invMass", "mcVis_invMass", 300, 0, 3 );

  _h_mc_aColin     = new TH2F( "mc_aColin" ,    "mc_aColin",     250, 0, 500, 400, 0, 3.14159);
  _h_mc_aCoplan    = new TH2F( "mc_aCoplan",    "mc_aCoplan",    250, 0, 500, 400, 0, 3.14159);

  _h_mcCheckTauMass = new TH1F( "mcCheckTauMass", "mcCheckTauMass", 200,1,3 );


  _h_mcPfoMatch_tauPi_nChNeu = new TH2F( "mcPfoMatch_tauPi_nChNeu", "mcPfoMatch_tauPi_nChNeu", 10,-0.5, 9.5, 10, -0.5, 9.5);
  _h_mcPfoMatch_tauPi_nChpi = new TH1F( "mcPfoMatch_tauPi_nChpi", "mcPfoMatch_tauPi_nChpi",  10,-0.5, 9.5 );
  _h_mcPfoMatch_tauPi_nNeuh = new TH1F( "mcPfoMatch_tauPi_nNeuh", "mcPfoMatch_tauPi_nNeuh",  10,-0.5, 9.5 );
  _h_mcPfoMatch_tauPi_nGamma = new TH1F( "mcPfoMatch_tauPi_nGamma", "mcPfoMatch_tauPi_nGamma",  10,-0.5, 9.5 );
  _h_mcPfoMatch_tauPi_nEl = new TH1F( "mcPfoMatch_tauPi_nEl", "mcPfoMatch_tauPi_nEl",  10,-0.5, 9.5 );
  _h_mcPfoMatch_tauPi_nMu = new TH1F( "mcPfoMatch_tauPi_nMu", "mcPfoMatch_tauPi_nMu",  10,-0.5, 9.5 );
  _h_mcPfoMatch_tauPi_nV0 = new TH1F( "mcPfoMatch_tauPi_nV0", "mcPfoMatch_tauPi_nV0",  10,-0.5, 9.5 );

  _h_mcPfoMatch_tauRho_nChNeu = new TH2F( "mcPfoMatch_tauRho_nChNeu", "mcPfoMatch_tauRho_nChNeu", 10,-0.5, 9.5, 10, -0.5, 9.5);
  _h_mcPfoMatch_tauRho_nChpi = new TH1F( "mcPfoMatch_tauRho_nChpi", "mcPfoMatch_tauRho_nChpi",  10,-0.5, 9.5 );
  _h_mcPfoMatch_tauRho_nNeuh = new TH1F( "mcPfoMatch_tauRho_nNeuh", "mcPfoMatch_tauRho_nNeuh",  10,-0.5, 9.5 );
  _h_mcPfoMatch_tauRho_nGamma = new TH1F( "mcPfoMatch_tauRho_nGamma", "mcPfoMatch_tauRho_nGamma",  10,-0.5, 9.5 );
  _h_mcPfoMatch_tauRho_nEl = new TH1F( "mcPfoMatch_tauRho_nEl", "mcPfoMatch_tauRho_nEl",  10,-0.5, 9.5 );
  _h_mcPfoMatch_tauRho_nMu = new TH1F( "mcPfoMatch_tauRho_nMu", "mcPfoMatch_tauRho_nMu",  10,-0.5, 9.5 );
  _h_mcPfoMatch_tauRho_nV0 = new TH1F( "mcPfoMatch_tauRho_nV0", "mcPfoMatch_tauRho_nV0",  10,-0.5, 9.5 );

  _h_mcPfoMatch_chgz0 = new TH1F("mcPfoMatch_chgz0","mcPfoMatch_chgz0",400,-2,2);
  _h_mcPfoMatch_chgz0WrtPv = new TH1F("mcPfoMatch_chgz0WrtPv", "mcPfoMatch_chgz0WrtPv" ,400,-2,2);


  _h_mcfitall_ttpz = new TH1F("mcfitall_ttpz","mcfitall_ttpz", 500,0,1000);
  _h_mcfitall_ttE  = new TH1F("mcfitall_ttE" ,"mcfitall_ttE" , 500,0,1000);
  _h_mcfitall_ttpt = new TH1F("mcfitall_ttpt","mcfitall_ttpt", 500,0,1000);



  _h_pfototen_mumu    = new TH1F("pfototen_mumu","pfototen_mumu",500,0,2000) ;
  _h_pfototen_tautau  = new TH1F("pfototen_tautau","pfototen_tautau",500,0,2000) ;

  for (int ien=0; ien<2; ien++) {
    TString hnen= ien==0 ? "_zm" : "_highm";

    _h_mcTauTauMassSel[ien] = new TH1F( "mcTauTauMass"+hnen, "mcTauTauMass"+hnen, 200, 0, 510 );

    _h_mcTauEn[ien] = new TH2F("mcTauEn"+hnen, "mcTauEn"+hnen, 100, 0, 260,100,0,260);
    _h_mcTauCosth[ien] = new TH2F("mcTauCosth"+hnen, "mcTauCosth"+hnen, 110, -1.05, 1.05, 110, -1.05, 1.05);

    _h_mcPiPiDecAng[ien]  = new TH2F("mcPiPiDecAng"+hnen,"mcPiPiDecAng"+hnen,20,0,3.15,20,0.,3.15);

    _h_mcPiPi_PiEn[ien]  = new TH2F("mcPiPi_PiEn"+hnen, "mcPiPi_PiEn"+hnen,30,0.,300.,30,0.,300.);
    _h_mcPiPi_PiEnFrac[ien]  = new TH2F("mcPiPi_PiEnFrac"+hnen, "mcPiPi_PiEnFrac"+hnen,20,0.,1.,20,0.,1.);

    _h_mcPiPi_polThTh[ien] = new TH2F("mcPiPi_polCosTh"+hnen, "mcPiPi_polCosTh"+hnen, 20,-1,1,20,-1,1);
    _h_mcPiRho_polThTh[ien] = new TH2F("mcPiRho_polCosTh"+hnen, "mcPiRho_polCosTh"+hnen, 20,-1,1,20,-1,1);
    _h_mcRhoPi_polThTh[ien] = new TH2F("mcRhoPi_polCosTh"+hnen, "mcRhoPi_polCosTh"+hnen, 20,-1,1,20,-1,1);
    _h_mcRhoRho_polThTh[ien] = new TH2F("mcRhoRho_polCosTh"+hnen, "mcRhoRho_polCosTh"+hnen, 20,-1,1,20,-1,1);

    _h_mcRhoRho_colin_polThTh[ien] = new TH2F("mcRhoRho_colin_polCosTh"+hnen, "mcRhoRho_colin_polCosTh"+hnen, 20,-1,1,20,-1,1);

    _h_mcPiPi_polComb1[ien] = new TH1F("mcPiPi_polComb1"+hnen, "mcPiPi_polComb1"+hnen, 40, -2, 2);
    _h_mcPiRho_polComb1[ien] = new TH1F("mcPiRho_polComb1"+hnen, "mcPiRho_polComb1"+hnen, 40, -2, 2);
    _h_mcRhoPi_polComb1[ien] = new TH1F("mcRhoPi_polComb1"+hnen, "mcRhoPi_polComb1"+hnen, 40, -2, 2);
    _h_mcRhoRho_polComb1[ien] = new TH1F("mcRhoRho_polComb1"+hnen, "mcRhoRho_polComb1"+hnen, 40, -2, 2);

    _h_mcPiPi_polComb2[ien] = new TH1F("mcPiPi_polComb2"+hnen, "mcPiPi_polComb2"+hnen, 40, -2, 2);
    _h_mcPiRho_polComb2[ien] = new TH1F("mcPiRho_polComb2"+hnen, "mcPiRho_polComb2"+hnen, 40, -2, 2);
    _h_mcRhoPi_polComb2[ien] = new TH1F("mcRhoPi_polComb2"+hnen, "mcRhoPi_polComb2"+hnen, 40, -2, 2);
    _h_mcRhoRho_polComb2[ien] = new TH1F("mcRhoRho_polComb2"+hnen, "mcRhoRho_polComb2"+hnen, 40, -2, 2);


    _h_pfo_id[ien] = new TH1F("pfo_id"+hnen,"pfo_id"+hnen,NTYPE, -0.5, NTYPE-0.5);
    for (int i=0; i<NTYPE; i++) {
      _h_pfo_id[ien]->GetXaxis()->SetBinLabel(i+1, pfoTypeStrings[i] );
    }

    _h_nPFO[ien]   = new TH1F("nPFO"+hnen,"nPFO"+hnen,100,0,100);

    _h_pfos_by_type[ien] = new TH2F("pfos_by_type"+hnen,"pfos_by_type"+hnen,NTAUDEC,-0.5, NTAUDEC-0.5, NTYPE, -0.5, NTYPE-0.5);
    for (int i=0; i<NTAUDEC; i++) {
      _h_pfos_by_type[ien]->GetXaxis()->SetBinLabel(i+1, tauDecStrings[i] );
    }
    for (int i=0; i<NTYPE; i++) {
      _h_pfos_by_type[ien]->GetYaxis()->SetBinLabel(i+1, pfoTypeStrings[i] );
    }

    for (int j=0; j<NTAUDEC; j++) {
      for (int i=0; i<NTYPE; i++) {
	TString hn = "TAUDEC_"+tauDecStrings[j]+"_PFOTYPE_"+pfoTypeStrings[i]+hnen;
	_h_pfos_by_dec_type[ien][j][i]=new TH1F(hn,hn,10,-0.5,9.5);
      }

      TString hn = "TAUDEC_"+tauDecStrings[j]+"_visMass"+hnen;
      _h_visMass[ien][j] = new TH1F(hn,hn,1000,0,4);
      hn = "TAUDEC_"+tauDecStrings[j]+"_ggMass"+hnen;
      _h_ggMass[ien][j] = new TH1F(hn,hn,1000,0,4);

      hn = "TAUDEC_"+tauDecStrings[j]+"_npfo"+hnen;
      _h_npfos[ien][j] = new TH1F(hn,hn,20,0,20);

    }

    _h_pfotauang[ien] = new TH2F("pfotauang"+hnen,"pfotauang"+hnen,50,0,3.2,50,0,3.2);
    _h_ntrk    [ien] = new TH1F("pfontrk"+hnen,"pfontrk"+hnen,5,0,5);
    _h_trkD0   [ien] = new TH1F("d0"+hnen,"d0"+hnen,1000,-5,5);
    _h_trkD0err[ien] = new TH1F("d0err"+hnen,"d0err"+hnen,1000,0,1);
    _h_trkD0sig[ien] = new TH1F("d0sig"+hnen,"d0sig"+hnen,1000,-60,60);


  }

  _h_mumu_nPFO         = new TH1F("h_mumu_nPFO"        , "h_mumu_nPFO"        , 100,0,100);
  _h_mumu_nPFOChg	 = new TH1F("h_mumu_nPFOChg"     , "h_mumu_nPFOChg"     , 40,0,40);
  _h_mumu_nPFONeu      = new TH1F("h_mumu_nPFONeu"     , "h_mumu_nPFONeu"     , 40,0,40);
  _h_mumu_pfoNeu_costh = new TH1F("h_mumu_pfoNeu_costh", "h_mumu_pfoNeu_costh", 100,-1,1);
  _h_mumu_pfoChg_costh = new TH1F("h_mumu_pfoChg_costh", "h_mumu_pfoChg_costh", 100,-1,1);
  _h_mumu_pfoNeu_en    = new TH1F("h_mumu_pfoNeu_en",    "h_mumu_pfoNeu_en"   , 50,0,25);
  _h_mumu_pfoChg_en    = new TH1F("h_mumu_pfoChg_en",    "h_mumu_pfoChg_en"   , 50,0,25);

  _h_mumu_pfoMCchad_costh = new TH1F("h_mumu_pfoMCchad_costh", "h_mumu_pfoMCchad_costh", 100,-1,1);
  _h_mumu_pfoMCnhad_costh = new TH1F("h_mumu_pfoMCnhad_costh", "h_mumu_pfoMCnhad_costh", 100,-1,1);
  _h_mumu_pfoMCgam_costh = new TH1F("h_mumu_pfoMCgam_costh", "h_mumu_pfoMCgam_costh", 100,-1,1);
  _h_mumu_pfoMCele_costh = new TH1F("h_mumu_pfoMCele_costh", "h_mumu_pfoMCele_costh", 100,-1,1);
  _h_mumu_pfoMCmuo_costh = new TH1F("h_mumu_pfoMCmuo_costh", "h_mumu_pfoMCmuo_costh", 100,-1,1);
  _h_mumu_pfoMCoth_costh = new TH1F("h_mumu_pfoMCoth_costh", "h_mumu_pfoMCoth_costh", 100,-1,1);

  _h_mumu_trk_costh = new TH1F("h_mumu_trk_costh","h_mumu_trk_costh",100,-1,1);
  _h_mumu_mcele_costh = new TH1F("h_mumu_mcele_costh","h_mumu_mcele_costh",100,-1,1);


  _mc_taujet_np = new TH1F("taujet_n","taujet_n",10,0,10);
  _mc_taujet_m  = new TH1F("taujet_m","taujet_m",100,0,10);
  _mc_taujet_e  = new TH1F("taujet_e","taujet_e",100,0,500);


  _h_cpcheck = new TH1F("cpangcheck","cpangcheck",100,-3.14159,3.14159) ;

  _h_cpcheck_leplep = new TH1F("cpangcheck_leplep","cpangcheck_leplep",20,-3.14159,3.14159) ;
  _h_cpcheck_pipi = new TH1F("cpangcheck_pipi","cpangcheck_pipi",20,-3.14159,3.14159) ;
  _h_cpcheck_rhorho = new TH1F("cpangcheck_rhorho","cpangcheck_rhorho",20,-3.14159,3.14159) ;
  _h_cpcheck_a3a3 = new TH1F("cpangcheck_a3a3","cpangcheck_a3a3",20,-3.14159,3.14159) ;
  _h_cpcheck_a1a1 = new TH1F("cpangcheck_a1a1","cpangcheck_a1a1",20,-3.14159,3.14159) ;

  _h_cpcheck2_leplep = new TH2F("cpangcheck2_leplep","cpangcheck2_leplep",20,0,1,20,-3.14159,3.14159) ;
  _h_cpcheck2_pipi   = new TH2F("cpangcheck2_pipi"  ,"cpangcheck2_pipi",  20,0,1,20,-3.14159,3.14159) ;
  _h_cpcheck2_rhorho = new TH2F("cpangcheck2_rhorho","cpangcheck2_rhorho",20,0,1,20,-3.14159,3.14159) ;
  _h_cpcheck2_a3a3   = new TH2F("cpangcheck2_a3a3"  ,"cpangcheck2_a3a3",  20,0,1,20,-3.14159,3.14159) ;
  _h_cpcheck2_a1a1   = new TH2F("cpangcheck2_a1a1"  ,"cpangcheck2_a1a1",  20,0,1,20,-3.14159,3.14159) ;

  _h_cpcheckPolAngle_a3a3 = new TH2F("cpcheckPolAngle_a3a3","cpcheckPolAngle_a3a3",20,-1,1,20,-1,1);
  _h_cpcheckDPhiDecPlane_a3a3 = new TH1F( "cpcheckDPhiDecPlane_a3a3", "cpcheckDPhiDecPlane_a3a3", 20,-3.14159,3.14159);
  _h_cpcheckDPhiDecPlane_a3a3_r1 = new TH1F( "cpcheckDPhiDecPlane_a3a3_r1", "cpcheckDPhiDecPlane_a3a3_r1", 20,-3.14159,3.14159);
  _h_cpcheckDPhiDecPlane_a3a3_r2 = new TH1F( "cpcheckDPhiDecPlane_a3a3_r2", "cpcheckDPhiDecPlane_a3a3_r2", 20,-3.14159,3.14159);
  _h_cpcheckDPhiDecPlane_a3a3_r3 = new TH1F( "cpcheckDPhiDecPlane_a3a3_r3", "cpcheckDPhiDecPlane_a3a3_r3", 20,-3.14159,3.14159);

  _h_cpcheckPolAngle_pipi = new TH2F("cpcheckPolAngle_pipi","cpcheckPolAngle_pipi",20,-1,1,20,-1,1);
  _h_cpcheckDPhiDecPlane_pipi = new TH1F( "cpcheckDPhiDecPlane_pipi", "cpcheckDPhiDecPlane_pipi", 20,-3.14159,3.14159);
  _h_cpcheckDPhiDecPlane_pipi_r1 = new TH1F( "cpcheckDPhiDecPlane_pipi_r1", "cpcheckDPhiDecPlane_pipi_r1", 20,-3.14159,3.14159);
  _h_cpcheckDPhiDecPlane_pipi_r2 = new TH1F( "cpcheckDPhiDecPlane_pipi_r2", "cpcheckDPhiDecPlane_pipi_r2", 20,-3.14159,3.14159);
  _h_cpcheckDPhiDecPlane_pipi_r3 = new TH1F( "cpcheckDPhiDecPlane_pipi_r3", "cpcheckDPhiDecPlane_pipi_r3", 20,-3.14159,3.14159);

  _h_cpcheckPolAngle_rhorho = new TH2F("cpcheckPolAngle_rhorho","cpcheckPolAngle_rhorho",20,-1,1,20,-1,1);
  _h_cpcheckDPhiDecPlane_rhorho = new TH1F( "cpcheckDPhiDecPlane_rhorho", "cpcheckDPhiDecPlane_rhorho", 20,-3.14159,3.14159);
  _h_cpcheckDPhiDecPlane_rhorho_r1 = new TH1F( "cpcheckDPhiDecPlane_rhorho_r1", "cpcheckDPhiDecPlane_rhorho_r1", 20,-3.14159,3.14159);
  _h_cpcheckDPhiDecPlane_rhorho_r2 = new TH1F( "cpcheckDPhiDecPlane_rhorho_r2", "cpcheckDPhiDecPlane_rhorho_r2", 20,-3.14159,3.14159);
  _h_cpcheckDPhiDecPlane_rhorho_r3 = new TH1F( "cpcheckDPhiDecPlane_rhorho_r3", "cpcheckDPhiDecPlane_rhorho_r3", 20,-3.14159,3.14159);

  Tauolapp::Tauola::initialize();
  string name="cteq6ll.LHpdf";
  LHAPDF::initPDFSetByName(name);

  // Next 3 lines are used to initialize TauSpinner
  //  CMSENE should be adjusted, Ipp    = true should be kept and
  // examples are set to work with unpolarized sammples only, Ipol is not
  // fully checked yet.

  // DANIEL thinks not important for this application (we use tauspinner just to calculate polarimeters)
  double CMSENE = 14000.; // center of mass system energy used in PDF calculation
  bool   Ipp    = true;   // for pp collisions
  int    Ipol   = 1;      // are input samples polarized?

  // Initialize TauSpinner (flags nonSM and nonSMN are 0 for these tests)
  TauSpinner::initialize_spinner(Ipp, Ipol, 0, 0, CMSENE);


  return;
}

void validateDST_TauTauProcessor::processRunHeader( LCRunHeader* ) {
  //  cout << "hello from validateDST_TauTauProcessor::processRunHeader" << endl;
}

TVector3 validateDST_TauTauProcessor::getTauolaPolMC(MCParticle* mctau) {
  assert(mctau);

  std::vector <std::pair <int, TLorentzVector> > tlv_daughters;

  //  double* TauSpinner::calculateHH(int tau_pdgid, vector<Particle> &tau_daughters, double phi, double theta);
  // TauSpinner::Particle ts_tau ( mctau->getMomentum()[0],  mctau->getMomentum()[1],  mctau->getMomentum()[2], 
  // 				mctau->getEnergy(), mctau->getPDG() );
  // TauSpinner::Particle ts_nu_tau;
  
  vector<TauSpinner::Particle> ts_tau_daughters;
  for ( size_t idau=0; idau<mctau->getDaughters().size(); idau++) {
    MCParticle* mcdau = mctau->getDaughters()[idau];

    TLorentzVector tlv( mcdau->getMomentum()[0],
			mcdau->getMomentum()[1],
			mcdau->getMomentum()[2],
			mcdau->getEnergy() );
    std::pair <int, TLorentzVector> tlvp( mcdau->getPDG(), tlv );
    
    //TauSpinner::Particle pp( mcdau->getMomentum()[0],mcdau->getMomentum()[1],mcdau->getMomentum()[2],
    // mcdau->getEnergy(),mcdau->getPDG() );

    // if ( abs(mcdau->getPDG())==16 ) {
    //   ts_nu_tau = pp;
    //   ts_tau_daughters.push_back(pp);
    //   cout << "got tau nu" << endl;
    // } else 
    if (  abs(mcdau->getPDG())!=15 &&  abs(mcdau->getPDG())!=24 && abs(mcdau->getPDG())!=20213  && abs(mcdau->getPDG())!=213) {
      cout << "daughter " << mcdau->getPDG() << endl;

      tlv_daughters.push_back( tlvp );

      // ts_tau_daughters.push_back(pp);
    } else {
      for ( size_t jdau=0; jdau<mcdau->getDaughters().size(); jdau++) {
	MCParticle* mcgdau = mcdau->getDaughters()[jdau];
	
	TLorentzVector tlv2(  mcgdau->getMomentum()[0], 
			      mcgdau->getMomentum()[1], 
			      mcgdau->getMomentum()[2],
			      mcgdau->getEnergy() );

	std::pair <int, TLorentzVector> tlvp2( mcgdau->getPDG(), tlv2 );

	// TauSpinner::Particle pp2(  mcgdau->getMomentum()[0], mcgdau->getMomentum()[1], mcgdau->getMomentum()[2],
	//    mcgdau->getEnergy(),      mcgdau->getPDG() );
	if (  abs(mcgdau->getPDG())!=15 &&  abs(mcgdau->getPDG())!=24 ) {
	  
	  cout << "grand daughter " << mcgdau->getPDG() << endl;
	  
	  tlv_daughters.push_back( tlvp2 );
	  
	  // ts_tau_daughters.push_back(pp2);
	} else {
	  cout << "WARNING, weird decay products" << mcgdau->getPDG() << endl;
	}
      }
    }
  }


  TLorentzVector tlv_total(0,0,0,0);
  for ( size_t i = 0; i<tlv_daughters.size(); i++) {
    tlv_total+=tlv_daughters[i].second;
  }
  for ( size_t i = 0; i<tlv_daughters.size(); i++) {
    tlv_daughters[i].second.Boost( - tlv_total.BoostVector() );

    TauSpinner::Particle pp( tlv_daughters[i].second.Px(),  
			     tlv_daughters[i].second.Py(), 
			     tlv_daughters[i].second.Pz(), 
			     tlv_daughters[i].second.E(), 
			     tlv_daughters[i].first ); 

    ts_tau_daughters.push_back( pp );

  }


  cout << "inputs:" << endl;
  for ( size_t i = 0; i<ts_tau_daughters.size(); i++) {
    cout << i << " " << ts_tau_daughters[i].pdgid() << " " << ts_tau_daughters[i].px()  << " " << ts_tau_daughters[i].py()  << " " << ts_tau_daughters[i].pz()  << endl;
  }


  double phi = 0.0, theta = 0.0;
  //  TauSpinner::prepareKinematicForHH   (ts_tau, ts_nu_tau, ts_tau_daughters, &phi, &theta);

  cout << "rotated" << endl;
  for ( size_t i = 0; i<ts_tau_daughters.size(); i++) {
    cout << i << " " << ts_tau_daughters[i].pdgid() << " " << ts_tau_daughters[i].px()  << " " << ts_tau_daughters[i].py()  << " " << ts_tau_daughters[i].pz()  << endl;
  }
  cout << "angles " << phi << " " << theta << endl;

  double* polarimeter = TauSpinner::calculateHH(  mctau->getPDG() , ts_tau_daughters, phi, theta);

  cout << "getTauolaPolMC polarimeter " << polarimeter[0] << " " << polarimeter[1] << " " << polarimeter[2] << endl;

  return TVector3( polarimeter[0] , polarimeter[1] , polarimeter[2]);
}



float validateDST_TauTauProcessor::getTauPiDecayAngle_MC( MCParticle* mctau ) {
  assert(mctau);

  MCParticle* mcpi(0);

  for ( size_t i=0; i<mctau->getDaughters().size(); i++) {
    if ( abs(mctau->getDaughters()[i]->getPDG())==211 ) {
      mcpi=mctau->getDaughters()[i];
      break;
    }
  }

  assert (mcpi);

  TLorentzVector tltau( mctau->getMomentum()[0], mctau->getMomentum()[1], mctau->getMomentum()[2], mctau->getEnergy() );
  TLorentzVector tlpi( mcpi->getMomentum()[0], mcpi->getMomentum()[1], mcpi->getMomentum()[2], mcpi->getEnergy() );

  tlpi.Boost( -tltau.BoostVector() );

  float decayangle = tlpi.Vect().Angle( tltau.Vect() );

  return decayangle;
}

std::vector <MCParticle*> validateDST_TauTauProcessor::getTauDaughters_byPdg_MC( MCParticle* mctau, int pdg, bool bothCharges ) {
  assert(mctau);
  std::vector <MCParticle*> daughters;

  for ( size_t i=0; i<mctau->getDaughters().size(); i++) {
    if ( (  bothCharges && abs(mctau->getDaughters()[i]->getPDG())==pdg ) ||
	 ( !bothCharges && mctau->getDaughters()[i]->getPDG()==pdg ) ||
	 pdg==0 ) {
      daughters.push_back( mctau->getDaughters()[i] );
    }
  }

  return daughters;
}


std::vector <MCParticle*> validateDST_TauTauProcessor::getTauGrandDaughters_byPdg_MC( MCParticle* mctau, int pdg, bool bothCharges ) {
  assert(mctau);
  std::vector <MCParticle*> gdaughters;

  for ( size_t i=0; i<mctau->getDaughters().size(); i++) {

    MCParticle* dau = mctau->getDaughters()[i];

    if ( dau->getGeneratorStatus()==2 ) { // unstable in generator, look for granddaughters

      for ( size_t j=0; j<dau->getDaughters().size(); j++) {

	MCParticle* gdau = dau->getDaughters()[j];

	if ( (  bothCharges && abs(gdau->getPDG())==pdg ) ||
	     ( !bothCharges && gdau->getPDG()==pdg ) ||
	     pdg==0 ) {
	  gdaughters.push_back( gdau );
	}
      }
    }
  }

  return gdaughters;
}


float validateDST_TauTauProcessor::getTauPiEnFrac_MC( MCParticle* mctau ) {
  assert(mctau);
  MCParticle* mcpi(0);
  for ( size_t i=0; i<mctau->getDaughters().size(); i++) {
    if ( abs(mctau->getDaughters()[i]->getPDG())==211 ) {
      mcpi=mctau->getDaughters()[i];
      break;
    }
  }
  assert (mcpi);
  return mcpi->getEnergy() / mctau->getEnergy() ;
}

float validateDST_TauTauProcessor::getTauPiEn_MC( MCParticle* mctau ) {
  assert(mctau);
  MCParticle* mcpi(0);
  for ( size_t i=0; i<mctau->getDaughters().size(); i++) {
    if ( abs(mctau->getDaughters()[i]->getPDG())==211 ) {
      mcpi=mctau->getDaughters()[i];
      break;
    }
  }
  assert (mcpi);
  return mcpi->getEnergy() ;
}

int validateDST_TauTauProcessor::getTauDecayMode_MC( MCParticle* mctau ) {


  assert ( abs( mctau->getPDG() ) ==15 );

  int a1chg    = 0;
  int rho      = 0;
  int k0Long   = 0;
  int k0Short  = 0;
  int kStarChg = 0;
  int kChg     = 0;
  int k0       = 0;
  int pichg    = 0;
  int pizero   = 0;
  int el       = 0;
  int mu       = 0;
  int nu       = 0;
  int W        = 0;
  int gamma    = 0;

  int nChStableHad=0;
  int nNeuStableHad=0;

  for ( size_t i=0; i<mctau->getDaughters().size(); i++) {

    MCParticle* mcdau=mctau->getDaughters()[i];

    cout << "tau daughter:: " << mcdau->getPDG() << " genstat " << mcdau->getGeneratorStatus() << " simstat " << mcdau->getSimulatorStatus () 
	 << " crInSim " << mcdau->isCreatedInSimulation () << " vtxNotParEnd " << mcdau->vertexIsNotEndpointOfParent () 
	 << " decInTrk " << mcdau->isDecayedInTracker () 
	 << " ; vtx " << mcdau->getVertex()[0] << " "<< mcdau->getVertex()[1] << " " << mcdau->getVertex()[2] << endl;

    cout << i << "/" << 
      mctau->getDaughters().size() << " " << 
      mcdau->getPDG() << " " << 
      mcdau->getGeneratorStatus() << " " << 
      mcdau->getDaughters().size() << endl;

    int pdg=abs(mcdau->getPDG());
    if ( pdg == 20213 ) {
      a1chg++;
    } else if ( pdg == 213 ) {
      rho++;
    } else if ( pdg == 323 ) {
      kStarChg++;
      nChStableHad++;
    } else if ( pdg == 321 ) {
      kChg++;
      nChStableHad++;
    } else if ( pdg == 130 ) {
      k0Long++;
      nNeuStableHad++;
    } else if ( pdg == 310 ) {
      k0Short++;
      nNeuStableHad++;
    } else if ( pdg == 211 ) {
      pichg++;
      nChStableHad++;
    } else if ( pdg == 111 ) {
      pizero++;
      nNeuStableHad++;
    } else if ( pdg == 11 ) {
      el++;
    } else if ( pdg == 13 ) {
      mu++;
    } else if ( pdg == 12 || pdg == 14 || pdg == 16 ) {
      nu++;
    } else if ( pdg == 24 ) {
      W++;
    } else if ( pdg == 22 ) {
      gamma++;
    }
    //else {cout << "unrecognised pdg..." << pdg << endl;}
    

  }


  int decayMode(-1);

  if        ( mu    >0                ) {
    //    cout << "tauDec_mu;   " << endl;  
    decayMode = tauDec_mu;    
  } else if ( el    >0                ) {
    //    cout << "tauDec_el;   " << endl;  
    decayMode = tauDec_el;
  } else if ( a1chg >0                ) {
    if ( pichg==3 ) {
      //      cout << "tauDec_a1chg; 3 prong" << endl;  
      decayMode = tauDec_a1chg_3p;
    } else if ( pichg==1 ) {
      //      cout << "tauDec_a1chg; 1 prong" << endl;  
      decayMode = tauDec_a1chg_1p;
    } else {

      if ( _isDBD ) {

	int nchpi(0);
	for ( size_t i=0; i<mctau->getDaughters().size(); i++) {
	  MCParticle* mcdau=mctau->getDaughters()[i];
	  if ( abs(  mcdau->getPDG() ) == 20213 ) {
	    for ( size_t j=0; j<mcdau->getDaughters().size(); j++) {
	      MCParticle* mcgdau = mcdau->getDaughters()[j];
	      if ( abs(mcgdau->getPDG()) == 211 ) nchpi++;
	    }
	  }
	}


	if ( nchpi==3 ) {
	  decayMode = tauDec_a1chg_3p;
	} else if  ( nchpi==1 ) {
	  decayMode = tauDec_a1chg_1p;
	} else {
	  cout << "weird a1 decay (dbd) with nchpi = " << nchpi << endl;
	  assert(0); return -1;
	}

      } else {

	cout << "weird a1 decay " << pichg << endl;  
      
	assert(0); return -1;
      }
    }
  } else if ( rho   >0                ) {
    //    cout << "tauDec_rho;  " << endl;  
    decayMode = tauDec_rho;
  } else if ( W     >0                ) {
    //    cout << "tauDec_hadW; " << endl;  
    decayMode = tauDec_hadW;
  } else if ( kStarChg >0 || kChg >0  ) {
    // cout << "tauDec_kchg; " << endl;  
    decayMode = tauDec_kchg;
  } else if ( pichg >0                ) {

    if ( pichg == 1 && nChStableHad==1 && nNeuStableHad==0 ) {
      decayMode = tauDec_pi;
    } else if ( pichg == 1 && nChStableHad==1 && nNeuStableHad==1 && pizero == 1 ) {
      decayMode = tauDec_rho;
    } else if ( pichg == 1 && nChStableHad==1 && nNeuStableHad==2 && pizero == 2 ) {
      decayMode = tauDec_a1chg_1p;
    } else if ( pichg == 3 && nChStableHad==3 && nNeuStableHad==0 ) {
      decayMode = tauDec_a1chg_3p;
    } else {
      decayMode = tauDec_hadW; // general hadronic decay
    }
  }

  if ( decayMode < 0 ) {
    cout << "un recognised tau decay mode" << endl;
    assert(0);
  }

  cout << "decay mode: " << decayMode << endl;
  printDecayMode(decayMode);
  return decayMode; 

}

void validateDST_TauTauProcessor::printDecayMode( int imode ) {
  switch ( imode ) {
  case  tauDec_mu       : {cout << imode << " tauDec_mu       " << endl; break;}
  case  tauDec_el       : {cout << imode << " tauDec_el	" << endl; break;}
  case  tauDec_a1chg_3p : {cout << imode << " tauDec_a1chg_3p	" << endl; break;}
  case  tauDec_a1chg_1p : {cout << imode << " tauDec_a1chg_1p	" << endl; break;}
  case  tauDec_rho      : {cout << imode << " tauDec_rho	" << endl; break;}
  case  tauDec_hadW     : {cout << imode << " tauDec_hadW	" << endl; break;}
  case  tauDec_kchg     : {cout << imode << " tauDec_kchg	" << endl; break;}
  case  tauDec_pi       : {cout << imode << " tauDec_pi       " << endl; break;}
  default               : { cout << "undefined mode " << imode << endl;}
  }
}

int validateDST_TauTauProcessor::gettype( int pdg ) {
  int apdg=abs(pdg);
  if      ( apdg==13 ) return type_mu;
  else if ( apdg==11 ) return type_el;
  else if ( apdg== 22 ) return type_gam;
  else if ( apdg== 211 ) return type_chad;
  else if ( apdg== 2112 ) return type_nhad;
  else if ( apdg== 310 ) return type_kshort;
  else if ( apdg== 3122 ) return type_lambda;
  else cout << "un recognised pfo type (pdg) " << pdg << endl;
  assert(0);
}

void validateDST_TauTauProcessor::mcana( std::vector <MCParticle*> cheatPart ) {

  //cout << "hello from mcana : " << endl;
  //cout << cheatPart.size() << endl;

  std::vector <MCParticle*> fidPart;

  for ( size_t i=0; i<cheatPart.size(); i++) {
    float mom(0);
    for (int j=0; j<3; j++) {
      mom += pow( cheatPart[i]->getMomentum()[j], 2 );
    }
    mom=sqrt(mom);
    float costh=mom>0 ? cheatPart[i]->getMomentum()[2]/mom : 0;
    if ( mom>0.1 && fabs(costh)<0.96 ) {
      fidPart.push_back( cheatPart[i] );
    }
  }

  //  cout << "--------------------" << endl;

//  for ( size_t i=0; i<fidPart.size(); i++) {
//    cout << i << " " << fidPart[i]->getPDG() << " " << fidPart[i]->getEnergy() << endl;
//  }

  // find highest momentum chg fid part
  float maxE(-1);
  float maxIndx(-1);
  for ( size_t i=0; i<fidPart.size(); i++) {
    if ( abs(fidPart[i]->getCharge())>0 ) {
      if ( fidPart[i]->getEnergy() > maxE ) {
	maxE=fidPart[i]->getEnergy();
	maxIndx=i;
      }
    }
  }

  //  cout << "MAX particle: " << maxIndx << " " << maxE << endl;

  if ( maxIndx>=0 ) {

    TVector3 maxmom( fidPart[maxIndx]->getMomentum()[0],
		     fidPart[maxIndx]->getMomentum()[1],
		     fidPart[maxIndx]->getMomentum()[2]);

    // now find highest energy in other phi hemisphere
    float nextE(-1);
    float nextIndx(-1);

    TVector3 nextmom;
    for ( size_t i=0; i<fidPart.size(); i++) {
      if (i==maxIndx) continue;
      if ( abs(fidPart[i]->getCharge())>0 ) {
	
	nextmom.SetXYZ( fidPart[i]->getMomentum()[0],
			fidPart[i]->getMomentum()[1],
			fidPart[i]->getMomentum()[2]);

	float dphi = fabs( nextmom.DeltaPhi( maxmom ) );
	if ( dphi > acos(-1)/2. ) {
	  if ( fidPart[i]->getEnergy() > nextE ) {
	    nextE=fidPart[i]->getEnergy();
	    nextIndx=i;
	  }
	}
      }
    }

    if (nextIndx>=0) {
      nextmom.SetXYZ( fidPart[nextIndx]->getMomentum()[0],
		      fidPart[nextIndx]->getMomentum()[1],
		      fidPart[nextIndx]->getMomentum()[2]);

      // now associate other particles to one (or none) of these
      std::vector < MCParticle* > tau1;
      std::vector < MCParticle* > tau2;

      tau1.push_back(fidPart[maxIndx]);
      tau2.push_back(fidPart[nextIndx]);

      TLorentzVector tl1;
      TLorentzVector tl2;
      TLorentzVector tltemp;

      tl1.SetVectM(maxmom,  fidPart[maxIndx]->getMass() );
      tl2.SetVectM(nextmom, fidPart[nextIndx]->getMass() );

      fidPart.erase( find(  fidPart.begin(),  fidPart.end(), tau1[0] ) );
      fidPart.erase( find(  fidPart.begin(),  fidPart.end(), tau2[0] ) );

      if ( fidPart.size()>0 ) {
	TVector3 thismom;
	for ( size_t i=0; i<fidPart.size(); i++) {
	  thismom.SetXYZ( fidPart[i]->getMomentum()[0],
			  fidPart[i]->getMomentum()[1],
			  fidPart[i]->getMomentum()[2]);
	  float ang1 = thismom.Angle( tl1.Vect() );
	  float ang2 = thismom.Angle( tl2.Vect() );
	  //	  cout << "angles " << ang1 << " " << ang2 << endl;
	  int itau(0);
	  if        ( ang1 < ang2 && ang1 < acos(-1)/4. ) {
	    itau=1;
	  } else if ( ang2 < ang1 && ang2 < acos(-1)/4. ) {
	    itau=2;
	  }
	  float mass = fidPart[i]->getMass();
	  tltemp.SetVectM(thismom, mass);
	  if ( itau==1 ) {
	    tl1 += tltemp;
	    tau1.push_back( fidPart[i] );
	    fidPart.erase(fidPart.begin() + i);
	  } else if ( itau==2 ) {
	    tl2 += tltemp;
	    tau2.push_back( fidPart[i] );
	    fidPart.erase(fidPart.begin() + i);
	  }
	}
      }

      //      cout << "jet nparts " << tau1.size() << " " << tau2.size() << endl;
      //      cout << "jet masses " << tl1.M() << " " << tl2.M() << endl;

      _mc_taujet_np->Fill(tau1.size());
      _mc_taujet_np->Fill(tau2.size());

      _mc_taujet_m->Fill(tl1.M());
      _mc_taujet_m->Fill(tl2.M());

      _mc_taujet_e->Fill(tl1.E());
      _mc_taujet_e->Fill(tl2.E());

    }

  }

  return;
}

bool validateDST_TauTauProcessor::dau_created_in_sim(  const MCParticle* mcp ) {
  bool interactedWithMaterial(false);
  for ( int k=0; k<mcp->getDaughters().size(); k++) {
    const MCParticle* dau = mcp->getDaughters()[k];
    if ( dau->isCreatedInSimulation() ) {
      interactedWithMaterial=true;
      break;
    }
  }
  return interactedWithMaterial;
}


std::vector <const MCParticle*> validateDST_TauTauProcessor::getStableTauDescendents( const MCParticle* mctau ) {

  bool weird=false;

  std::vector <const MCParticle*> tauStableDescendents;

  for ( int j=0; j<mctau->getDaughters().size(); j++) {
    const MCParticle* dau=mctau->getDaughters()[j];


    if ( dau->getGeneratorStatus()==1 || dau_created_in_sim(dau) ) tauStableDescendents.push_back( dau );
    else {
      for ( int k=0; k<dau->getDaughters().size(); k++) {
	const MCParticle* gdau = dau->getDaughters()[k];
	if ( gdau->getGeneratorStatus()==1 || dau_created_in_sim(gdau) ) tauStableDescendents.push_back( gdau );
	else {
	  for ( int l=0; l<gdau->getDaughters().size(); l++) {
	    const MCParticle* ggdau = gdau->getDaughters()[l];

	    if ( ggdau->getGeneratorStatus()==1 || dau_created_in_sim(ggdau) ) tauStableDescendents.push_back( ggdau );
	    else {
	      for ( int m=0; m<ggdau->getDaughters().size(); m++) {
		const MCParticle* gggdau = ggdau->getDaughters()[m];
		if ( gggdau->getGeneratorStatus()==1 || dau_created_in_sim(gggdau) ) tauStableDescendents.push_back( gggdau );
		else {
		  cout << "getStableTauDescendents::WARNING need more generations!!!! " <<
		    gggdau->getPDG() << " " << gggdau->getGeneratorStatus() << endl;
		  weird = true;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  if ( weird ) tauStableDescendents.clear();

  return tauStableDescendents;
}


void validateDST_TauTauProcessor::processEvent( LCEvent * evt ) {

  float mccosth(-99);
  float mcphi(-99);

  //  if ( evt->getEventNumber()%10000 == 0 )
    cout << "new event " << evt->getEventNumber() << endl;


  _h_eventCounter->Fill(0);

  //---------------------------------
  // MC information
  //---------------------------------

  MCParticle* tauMC[2]={0,0};
  //  MCParticle* tauPlus(0);
  MCParticle* muMC[2]={0,0};
  //MCParticle* muPlus(0);

  LCRelationNavigator* relNavi(0);
  try {
    LCCollection* linkcol = evt->getCollection( "RecoMCTruthLink" );
    relNavi = new LCRelationNavigator(linkcol);
  } catch(DataNotAvailableException &e) {};


  int has94=0;
  int ntaunu=0;

  float mctoten(0);
  float pfototen(0);

  std::vector <MCParticle*> cheatedParticles_nooverlay;
  std::vector <MCParticle*> cheatedParticles;


  std::string colname = _mcskim==0 ? "MCParticle" : "MCParticlesSkimmed";

  try {


    //    LCCollection* mccol = evt->getCollection( "MCParticlesSkimmed" );
    LCCollection* mccol = evt->getCollection( colname );
    assert ( mccol );
    for (int j=0; j<mccol->getNumberOfElements(); j++) {
      MCParticle* mcp = dynamic_cast<MCParticle*> (mccol->getElementAt(j));

      if ( mcp->getPDG()==94 ) {
	has94=1;
      }

      if ( abs(  mcp->getPDG() )==16 ) {
	ntaunu++;
      }

      if ( mcp->isCreatedInSimulation()==0 && mcp->getGeneratorStatus() == 1 ) { 
	cheatedParticles.push_back(mcp);
	if ( mcp->isOverlay()==0 ) {
	  mctoten+=mcp->getEnergy();
	  cheatedParticles_nooverlay.push_back(mcp);
	}
      }

      if ( mcp->getPDG()==-13 && !muMC[0] ) {
	muMC[0]=mcp;
      }
      if ( mcp->getPDG()==13 && !muMC[1] ) {
	muMC[1]=mcp;
      }

      // find the taus before they decay
      if ( mcp->getDaughters().size()>1 ) {
	if ( mcp->getPDG()==-15 ) {
	  tauMC[0]=mcp;
	} else if ( mcp->getPDG()==15 ) {
	  tauMC[1]=mcp;
	}
      }
    }
  } catch(DataNotAvailableException &e) {};



  mcana( cheatedParticles_nooverlay );


  _h_mctaunu->Fill(ntaunu);

  if ( tauMC[0] && tauMC[1] ) {
    _h_mctoten_tautau ->Fill(mctoten);
  } else {
    _h_mctoten_mumu   ->Fill(mctoten);
  }

  if ( !tauMC[0] && !tauMC[1] ) {  // probably mu-mu+

    if ( !muMC[0] || !muMC[1] ) {
      cout << "ERROR no pair of taus or muons!" << endl;
      assert(0);
    }

    _h_eventCounter->Fill(1);


    TVector3 mcMuMom[2] = {muMC[0]->getMomentum(), muMC[1]->getMomentum() };
    TLorentzVector mcMuTLV[2];
    mcMuTLV[0].SetVectM(mcMuMom[0], 0.105);
    mcMuTLV[1].SetVectM(mcMuMom[1], 0.105);
    float mmmass = (mcMuTLV[0]+mcMuTLV[1]).M();
    _h_mcMuMuMass->Fill( mmmass );

    try {
      LCCollection* pfocol = evt->getCollection( "PandoraPFOs" );
      _h_mumu_nPFO->Fill(pfocol->getNumberOfElements());

      int nPfoNeu(0);
      int nPfoChg(0);
      for (int j=0; j<pfocol->getNumberOfElements(); j++) {
	ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (pfocol->getElementAt(j));
	float costh=pfo->getMomentum()[2]/sqrt( 
					       pow( pfo->getMomentum()[0],2) + 
					       pow( pfo->getMomentum()[1],2) + 
					       pow( pfo->getMomentum()[2],2) );
	if ( pfo->getCharge()==0 ) {
	  nPfoNeu++;
	  _h_mumu_pfoNeu_costh->Fill(costh);
	  _h_mumu_pfoNeu_en->Fill(pfo->getEnergy());
	} else {
	  nPfoChg++;
	  _h_mumu_pfoChg_costh->Fill(costh);
	  _h_mumu_pfoChg_en->Fill(pfo->getEnergy());
	}


	int bestTrkPdg=0;
	float bestTrkWt=0;
	int bestCalPdg=0;
	float bestCalWt=0;

	for (size_t k=0; k<relNavi->getRelatedToObjects(pfo).size(); k++) {
	  MCParticle* mcp = dynamic_cast <MCParticle*> (relNavi->getRelatedToObjects(pfo)[k]);
	  float wgt = relNavi->getRelatedToWeights(pfo)[k];
	  int pdg = abs( mcp->getPDG() );

	  float trackwgt = (int(wgt)%10000)/1000.;
	  float clusterwgt = (int(wgt)/10000)/1000. ;

	  //	  cout << k << " : " << pdg << " " << trackwgt << " "<< clusterwgt << endl;

	  if ( trackwgt>bestTrkWt ) {
	    bestTrkWt=trackwgt;
	    bestTrkPdg=pdg;
	  }
	  if ( clusterwgt > bestCalWt ) {
	    bestCalWt = clusterwgt;
	    bestCalPdg=pdg;
	  }
	}

	int bestpdg(0);
	if ( bestTrkWt > 0.5 ) {
	  bestpdg=bestTrkPdg;
	} else if (bestTrkWt+bestCalWt>0) {
	  bestpdg = bestTrkWt>bestCalWt ? bestTrkPdg : bestCalPdg;
	}

	//	cout << "best pdg " << bestpdg << endl;

	if ( bestpdg==22 ) {
	  _h_mumu_pfoMCgam_costh->Fill(costh);
	} else if ( bestpdg==211 || bestpdg==321 ||  bestpdg==2212 ) {
	  _h_mumu_pfoMCchad_costh->Fill(costh);
	} else if ( bestpdg==310 || bestpdg==2112 ) {
	  _h_mumu_pfoMCnhad_costh->Fill(costh);
	} else if ( bestpdg==11 ) {
	  _h_mumu_pfoMCele_costh->Fill(costh);
	} else if ( bestpdg==13 ) {
	  _h_mumu_pfoMCmuo_costh->Fill(costh);
	} else {
	  //cout << "did not know how to match to " << bestpdg << endl;
	  _h_mumu_pfoMCoth_costh->Fill(costh);
	}

      }
      
      _h_mumu_nPFOChg->Fill(nPfoChg);
      _h_mumu_nPFONeu->Fill(nPfoNeu);

    } catch(DataNotAvailableException &e) {};


    try {
      LCCollection* trkcol = evt->getCollection( "MarlinTrkTracks" );
      for (int j=0; j<trkcol->getNumberOfElements(); j++) {
	Track* trk = dynamic_cast<Track*> (trkcol->getElementAt(j));

	float tanl = trk->getTanLambda();
	float costh = tanl/sqrt( 1 + tanl*tanl);

	_h_mumu_trk_costh->Fill(costh);

      }
    } catch(DataNotAvailableException &e) {};


    try {
      LCCollection* mccol = evt->getCollection( "MCParticlesSkimmed" );
      assert ( mccol );
      for (int j=0; j<mccol->getNumberOfElements(); j++) {
	MCParticle* mcp = dynamic_cast<MCParticle*> (mccol->getElementAt(j));
	if ( !mcp->isCreatedInSimulation() && abs( mcp->getPDG() ) == 11 ) {

	  float costh = mcp->getMomentum()[2]/sqrt( 
						   pow( mcp->getMomentum()[0],2) + 
						   pow( mcp->getMomentum()[1],2) + 
						   pow( mcp->getMomentum()[2],2) 
						    );

	  _h_mumu_mcele_costh->Fill(costh);
	}
      }
    } catch(DataNotAvailableException &e) {};




  } else if ( tauMC[0] && tauMC[1] ) {

    _h_eventCounter->Fill(2);

    // this is the tau tau analysis

    const float lightSpeed = 2.99e8 * 1.e3 * 1.e-12; // in mm / ps
    float vtxr[2];
    float vtxdecl[2];
    float life[2];
    for (int i=0; i<2; i++) {
      vtxr[i]= sqrt( tauMC[i]->getEndpoint()[0]* tauMC[i]->getEndpoint()[0] + tauMC[i]->getEndpoint()[1]* tauMC[i]->getEndpoint()[1] );
      vtxdecl[i]=  sqrt( tauMC[i]->getEndpoint()[0]* tauMC[i]->getEndpoint()[0] + 
			  tauMC[i]->getEndpoint()[1]* tauMC[i]->getEndpoint()[1] +  
			  tauMC[i]->getEndpoint()[2]* tauMC[i]->getEndpoint()[2] );

      float gamma = tauMC[i]->getEnergy() / tauMC[i]->getMass();
      life[i]= vtxdecl[i] / ( gamma*lightSpeed );
    }


    _h_mcTauPlusVtxRZ ->Fill(vtxr[1],  tauMC[1]->getEndpoint()[2] );
    _h_mcTauPlusVtxXY ->Fill(tauMC[1]->getEndpoint()[0], tauMC[1]->getEndpoint()[1]);
    _h_mcTauMinusVtxRZ->Fill(vtxr[0], tauMC[0]->getEndpoint()[2]);
    _h_mcTauMinusVtxXY->Fill(tauMC[0]->getEndpoint()[0],  tauMC[0]->getEndpoint()[1] );

    _h_mcTauPlusDecl->Fill(vtxdecl[1]);
    _h_mcTauMinusDecl->Fill(vtxdecl[0]);

    _h_mcTauPlusLife->Fill(life[1]);
    _h_mcTauMinusLife->Fill(life[0]);

    _h_mcHas94->Fill( has94 );

    int decMode[2];
    for (int i=0; i<2; i++) {
      decMode[i]=getTauDecayMode_MC(tauMC[i]);
      _h_mcDecMode->Fill(decMode[i]);
    }
    _h_mcDecMode2d->Fill(decMode[0], decMode[1]);


    // try getting TauSpinner polarimeter from MC info

    TVector3 tauola_hel_MC[2];
    for (int i=0; i<2; i++) {
      tauola_hel_MC[i] = getTauolaPolMC(tauMC[i]);
    }

    // choose a reference direction (to calculate the trans spin correls)
    // tau0 direction in tau1 restframe
    TLorentzVector temp0(  tauMC[0]->getMomentum()[0],  tauMC[0]->getMomentum()[1],  tauMC[0]->getMomentum()[2], tauMC[0]->getEnergy() );
    TLorentzVector temp1(  tauMC[1]->getMomentum()[0],  tauMC[1]->getMomentum()[1],  tauMC[1]->getMomentum()[2], tauMC[1]->getEnergy() );
    temp1.Boost( - temp0.BoostVector() );
    TVector3 refDir = temp1.Vect()*(1./temp1.Vect().Mag());
    TVector3 transComp[2];

    refDir.Print();

      // get transverse components
    for (int itau=0; itau<2; itau++) {
      transComp[itau] = tauola_hel_MC[itau];
      transComp[itau] -= refDir*(transComp[itau].Dot(refDir));
      cout << "check dot " << transComp[itau].Dot(refDir) << endl;
    }    
    float cross = transComp[0].Cross( transComp[1] ).Dot( refDir );
    float theangle=transComp[0].Angle(transComp[1]);
    cout << "angle between trans comps " << transComp[0].Angle(transComp[1]) << " " << transComp[0].Cross( transComp[1] ).Dot( refDir )  << endl;
    if ( cross<0 ) theangle*=-1;
    _h_cpcheck->Fill(theangle);

    // get the contrast fn
    float ang0 = tauola_hel_MC[0].Angle( refDir );
    float ang1 = tauola_hel_MC[1].Angle( refDir );
    float constrast = sin(ang0)*sin(ang1)/( 1 + cos(ang0)*cos(ang1) );

    //  enum {tauDec_mu=0, tauDec_el, tauDec_a1chg_3p, tauDec_a1chg_1p, tauDec_rho, tauDec_hadW, tauDec_kchg, tauDec_pi, NTAUDEC};

    if ( tauola_hel_MC[0].Mag()>0 && tauola_hel_MC[1].Mag()>0 ) {

      if ( ( decMode[0] == tauDec_mu || decMode[0] == tauDec_el ) && 
	   ( decMode[1] == tauDec_mu || decMode[1] == tauDec_el ) ) {
	_h_cpcheck_leplep->Fill(theangle);
	_h_cpcheck2_leplep->Fill(constrast, theangle);
      } else if ( decMode[0] == tauDec_pi && decMode[1] == tauDec_pi ) {
	_h_cpcheck_pipi->Fill(theangle);
	_h_cpcheck2_pipi->Fill(constrast, theangle);
      } else if ( decMode[0] == tauDec_rho && decMode[1] == tauDec_rho ) {
	_h_cpcheck_rhorho->Fill(theangle);
	_h_cpcheck2_rhorho->Fill(constrast, theangle);
      } else if ( decMode[0] == tauDec_a1chg_3p && decMode[1] == tauDec_a1chg_3p ) {
	_h_cpcheck_a3a3->Fill(theangle);
	_h_cpcheck2_a3a3->Fill(constrast, theangle);
      } else if ( decMode[0] == tauDec_a1chg_1p && decMode[1] == tauDec_a1chg_1p ) {
	_h_cpcheck_a1a1->Fill(theangle);
	_h_cpcheck2_a1a1->Fill(constrast, theangle);
      }

    }


    if ( 
	( decMode[0] == tauDec_a1chg_3p && decMode[1] == tauDec_a1chg_3p ) || 
	( decMode[0] == tauDec_pi && decMode[1] == tauDec_pi ) ||
	( decMode[0] == tauDec_rho && decMode[1] == tauDec_rho ) 
	 ) {

      TLorentzVector tltau[2];
      std::vector <const MCParticle*> stable_tau_desc[2];
      TLorentzVector tlNuTau[2];
      for (int itau=0; itau<2; itau++) {
	tltau[itau].SetXYZT( tauMC[itau]->getMomentum()[0],
			     tauMC[itau]->getMomentum()[1],
			     tauMC[itau]->getMomentum()[2],
			     tauMC[itau]->getEnergy() );

	stable_tau_desc[itau] = getStableTauDescendents( tauMC[itau] );

        for ( size_t j=0; j< stable_tau_desc[itau].size(); j++) {
	  if ( fabs(  stable_tau_desc[itau][j]->getPDG() )==16 ) {
	    tlNuTau[itau].SetXYZT(  stable_tau_desc[itau][j]->getMomentum()[0],
				    stable_tau_desc[itau][j]->getMomentum()[1],
				    stable_tau_desc[itau][j]->getMomentum()[2],
				    stable_tau_desc[itau][j]->getEnergy() );
	  }
	}
      }
      // now boost to higgs restframe
      TVector3 boostVector = -(tltau[0]+tltau[1]).BoostVector();

      TLorentzVector tlNuTau_higgsRF[2];
      TLorentzVector tlTau_higgsRF[2];
      TLorentzVector tlNuTau_tauRF[2];

      float polAngle[2];
      TVector3 decayPlane[2];

      for (int itau=0; itau<2; itau++) {
	tlTau_higgsRF[itau]=tltau[itau];	tlTau_higgsRF[itau].Boost( boostVector );
	tlNuTau_higgsRF[itau]=tlNuTau[itau];	tlNuTau_higgsRF[itau].Boost( boostVector );
	tlNuTau_tauRF[itau]=tlNuTau[itau];	tlNuTau_tauRF[itau].Boost(  -tltau[itau].BoostVector() );

	polAngle[itau] = cos(tlNuTau_tauRF[itau].Angle( tltau[itau].Vect() ) );
	decayPlane[itau] = tlNuTau_higgsRF[itau].Vect().Cross( tlTau_higgsRF[itau].Vect() );
      }

      float deltaPhi =  decayPlane[0].Angle( decayPlane[1] ) ;
      float deltaPhiSign =  decayPlane[0].Cross( decayPlane[1] ).Dot( tlTau_higgsRF[0].Vect() );
      if ( deltaPhiSign < 0 ) {
	deltaPhi*=-1.;
      }




      if ( decMode[0] == tauDec_a1chg_3p && decMode[1] == tauDec_a1chg_3p ) {

	_h_cpcheckPolAngle_a3a3->Fill( polAngle[0], polAngle[1] );
	_h_cpcheckDPhiDecPlane_a3a3 -> Fill( deltaPhi );

	if ( polAngle[0]>0 && polAngle[1]>0 )
	  _h_cpcheckDPhiDecPlane_a3a3_r1 -> Fill( deltaPhi );
	else if ( polAngle[0]<0 && polAngle[1]<0 )
	  _h_cpcheckDPhiDecPlane_a3a3_r2 -> Fill( deltaPhi );
	else
	  _h_cpcheckDPhiDecPlane_a3a3_r3 -> Fill( deltaPhi );

      } else if ( decMode[0] == tauDec_pi && decMode[1] == tauDec_pi ) {

	_h_cpcheckPolAngle_pipi->Fill( polAngle[0], polAngle[1] );
	_h_cpcheckDPhiDecPlane_pipi -> Fill( deltaPhi );

	if ( polAngle[0]>0 && polAngle[1]>0 )
	  _h_cpcheckDPhiDecPlane_pipi_r1 -> Fill( deltaPhi );
	else if ( polAngle[0]<0 && polAngle[1]<0 )
	  _h_cpcheckDPhiDecPlane_pipi_r2 -> Fill( deltaPhi );
	else
	  _h_cpcheckDPhiDecPlane_pipi_r3 -> Fill( deltaPhi );

      } else if ( decMode[0] == tauDec_rho && decMode[1] == tauDec_rho ) {

	_h_cpcheckPolAngle_rhorho->Fill( polAngle[0], polAngle[1] );
	_h_cpcheckDPhiDecPlane_rhorho -> Fill( deltaPhi );

	if ( polAngle[0]>0 && polAngle[1]>0 )
	  _h_cpcheckDPhiDecPlane_rhorho_r1 -> Fill( deltaPhi );
	else if ( polAngle[0]<0 && polAngle[1]<0 )
	  _h_cpcheckDPhiDecPlane_rhorho_r2 -> Fill( deltaPhi );
	else
	  _h_cpcheckDPhiDecPlane_rhorho_r3 -> Fill( deltaPhi );

      }



    }


    //--------------------------


    TVector3 mcTauMom[2] = {tauMC[0]->getMomentum(), tauMC[1]->getMomentum() };

    const float tauMass = 1.777;

    TLorentzVector mcTauTLV[2];
    mcTauTLV[0].SetVectM(mcTauMom[0], tauMass);
    mcTauTLV[1].SetVectM(mcTauMom[1], tauMass);

    float ttmass = (mcTauTLV[0]+mcTauTLV[1]).M();
    _h_mcTauTauMass->Fill( ttmass );

    //    cout << "***** TRUE tau tau mass " << ttmass << endl;


    TLorentzVector mcTauVisTLV[2];

    std::vector <const MCParticle*> stable_tau_desc[2];

    for (int i=0; i<2; i++) {
      mcTauVisTLV[i].SetXYZT(0,0,0,0);
      stable_tau_desc[i] = getStableTauDescendents( tauMC[i] );
      if (  stable_tau_desc[i].size()==0 ) { // weird event....
	cout << "WEIRD EVENT (MC) ?? " << evt->getEventNumber() << endl;
	LCCollection* mccol = evt->getCollection( colname );
	assert ( mccol );
	for (int j=0; j<mccol->getNumberOfElements(); j++) {
	  MCParticle* mcp = dynamic_cast<MCParticle*> (mccol->getElementAt(j));
	  cout << j << " " << mcp << " " << mcp->getPDG() << " " << mcp->getGeneratorStatus() << endl;
	}
      }
      TLorentzVector mcTauCheckTLV(0,0,0,0);
      for ( size_t j=0; j< stable_tau_desc[i].size(); j++) {
	const MCParticle* mm =  stable_tau_desc[i][j];
	TLorentzVector aa( mm->getMomentum()[0],
			   mm->getMomentum()[1],
			   mm->getMomentum()[2],
			   mm->getEnergy() );
	mcTauCheckTLV+=aa;
	if ( abs( mm->getPDG() ) !=12 &&  abs( mm->getPDG() ) !=14 &&  abs( mm->getPDG() ) != 16 ) {
	  mcTauVisTLV[i]+=aa;
	} else {
	  //	  cout << "MC neutrino from tau " << i << " : " << mm->getPDG() << " " << mm->getEnergy() << endl;
	}
      }
      _h_mcVis_invMass->Fill(mcTauVisTLV[i].M());
      _h_mcCheckTauMass->Fill( mcTauCheckTLV.M() );
    }

    _h_mcVis_aColin -> Fill ( ttmass, 3.14159 - mcTauVisTLV[0].Vect().Angle( mcTauVisTLV[1].Vect() ) );
    _h_mcVis_aCoplan -> Fill( ttmass, 3.14159 - mcTauVisTLV[0].DeltaPhi( mcTauVisTLV[1] ) );
    _h_mc_aColin -> Fill ( ttmass, 3.14159 - mcTauTLV[0].Vect().Angle( mcTauTLV[1].Vect() ) );
    _h_mc_aCoplan -> Fill( ttmass, 3.14159 - mcTauTLV[0].DeltaPhi( mcTauTLV[1] ) );




//    for (int itau=0; itau<2; itau++) {
//
//      //  double* TauSpinner::calculateHH(int tau_pdgid, vector<Particle> &tau_daughters, double phi, double theta);
//      TauSpinner::Particle ts_tau ( tauMC[itau]->getMomentum()[0],  tauMC[itau]->getMomentum()[1],  tauMC[itau]->getMomentum()[2], 
//				 tauMC[itau]->getEnergy(), tauMC[itau]->getPDG() );
//      TauSpinner::Particle ts_nu_tau;
//
//      vector<TauSpinner::Particle> ts_tau_daughters;
//      for ( size_t idau=0; idau<tauMC[itau]->getDaughters().size(); idau++) {
//	MCParticle* mcdau = tauMC[itau]->getDaughters()[idau];
//
//
//	TauSpinner::Particle pp( mcdau->getMomentum()[0],mcdau->getMomentum()[1],mcdau->getMomentum()[2],
//		     mcdau->getEnergy(),mcdau->getPDG() );
//
//	if ( abs(mcdau->getPDG())==16 ) {
//	  ts_nu_tau = pp;
//	} else if (  abs(mcdau->getPDG())!=15 &&  abs(mcdau->getPDG())!=24 ) {
//	  ts_tau_daughters.push_back(pp);
//	} else {
//	  for ( size_t jdau=0; jdau<mcdau->getDaughters().size(); jdau++) {
//	    MCParticle* mcgdau = mcdau->getDaughters()[jdau];
//	    TauSpinner::Particle pp2(  mcgdau->getMomentum()[0], mcgdau->getMomentum()[1], mcgdau->getMomentum()[2],
//			   mcgdau->getEnergy(),      mcgdau->getPDG() );
//	    if (  abs(mcgdau->getPDG())!=15 &&  abs(mcgdau->getPDG())!=24 ) {
//	      ts_tau_daughters.push_back(pp2);
//	    } else {
//	      cout << "WARNING, weird decay products" << mcgdau->getPDG() << endl;
//	    }
//	  }
//	}
//      }
//
//      double phi = 0.0, theta = 0.0;
//      //TauSpinner::prepareKinematicForHH   (ts_tau, ts_nu_tau, ts_tau_daughters, &phi, &theta);
//      double* polarimeter = TauSpinner::calculateHH(  tauMC[itau]->getPDG() , ts_tau_daughters, phi, theta);
//
//      cout << "tau " << itau << " polarimeter " << polarimeter[0] << " " << polarimeter[1] << " " << polarimeter[2] << endl;
//
//      // get transverse components
//
//
//      transComp[itau].SetXYZ( polarimeter[0], polarimeter[1], polarimeter[2] );
//
//      transComp[itau] -= refDir*(transComp[itau].Dot(refDir));
//
//
//      cout << "check dot " << transComp[itau].Dot(refDir) << endl;
//
//      delete polarimeter;
//
//    }    


//    cout << "angle between trans comps " << transComp[0].Angle(transComp[1]) << " " << transComp[0].Cross( transComp[1] ).Dot( refDir )  << endl;



    TVector3 polarimeters[2];
    float thetaPol[2];
    float thetaPolColin[2]={-1,-1};
    for (int i=0; i<2; i++) {
      if ( decMode[i] == tauDec_pi ) {
	std::vector <MCParticle*> mcds = getTauDaughters_byPdg_MC( tauMC[i], 211, true );
	if ( mcds.size()!=1 || abs(mcds[0]->getPDG()) != 211 ) {
	  cout << "crazy tau -> pion decay ????" << endl;
	  for (size_t kk=0; kk<mcds.size(); kk++) 
	    cout << kk << " " << mcds[kk]->getPDG() << endl;
	  assert(0);
	}
	// get 4-mom of MC pion 
	TLorentzVector tlpi( mcds[0]->getMomentum()[0], mcds[0]->getMomentum()[1], mcds[0]->getMomentum()[2], mcds[0]->getEnergy() ); 
	// transform to rest fram of tau
	tlpi.Boost( -mcTauTLV[i].BoostVector() );

	//TVector3 polvec=tlpi.Vect();
	//polvec*=1./polvec.Mag();

	polarimeters[i]=tlpi.Vect();
	polarimeters[i]*=1./polarimeters[i].Mag();


	//	cout << "pi mode polarimeter: " << endl; 
	//	polvec.Print();
	//	tauola_hel_MC[i].Print();
	//	cout << "    angle to tauola " << polvec.Angle(	tauola_hel_MC[i] ) << endl;


	thetaPol[i] = polarimeters[i].Angle( mcTauTLV[i].Vect() );
      } else if ( decMode[i] == tauDec_rho ) {

	std::vector <MCParticle*> mcd_chpi;
	std::vector <MCParticle*> mcd_neupi;

	if ( _isPythiaTauDecay ) {
	  mcd_chpi =  getTauDaughters_byPdg_MC( tauMC[i], 211, true );
	  mcd_neupi =  getTauDaughters_byPdg_MC( tauMC[i], 111, false );
	} else {
	  mcd_chpi =  getTauGrandDaughters_byPdg_MC( tauMC[i], 211, true );
	  mcd_neupi =  getTauGrandDaughters_byPdg_MC( tauMC[i], 111, false );
	}

	std::vector <MCParticle*> mcd_neutrino  =  getTauDaughters_byPdg_MC( tauMC[i], 16, true );
	if ( mcd_chpi.size()!=1 || mcd_neupi.size()!=1 || mcd_neutrino.size()!=1 ) {
	  cout << "crazy tau -> rho decay ????" << mcd_chpi.size() << " " << mcd_neupi.size() << " " << mcd_neutrino.size() << endl;
	  assert(0);
        }
	TLorentzVector tlpiChg( mcd_chpi[0]->getMomentum()[0],  mcd_chpi[0]->getMomentum()[1],  mcd_chpi[0]->getMomentum()[2],  mcd_chpi[0]->getEnergy() );
	TLorentzVector tlpiNeu( mcd_neupi[0]->getMomentum()[0],  mcd_neupi[0]->getMomentum()[1],  mcd_neupi[0]->getMomentum()[2],  mcd_neupi[0]->getEnergy() );
	TLorentzVector tlNeutrino( mcd_neutrino[0]->getMomentum()[0],  
				   mcd_neutrino[0]->getMomentum()[1],
				   mcd_neutrino[0]->getMomentum()[2],
				   mcd_neutrino[0]->getEnergy() );
	// transform to tau rest frame
	tlpiChg.Boost( -mcTauTLV[i].BoostVector() );
	tlpiNeu.Boost( -mcTauTLV[i].BoostVector() );
	tlNeutrino.Boost( -mcTauTLV[i].BoostVector() );
	// now calculate the polarimeter vector 

	//	TVector3 polvec = 
	polarimeters[i] = 
	  tauMC[i]->getMass() * (tlpiChg.E() - tlpiNeu.E() ) * ( tlpiChg.Vect() - tlpiNeu.Vect() ) +
	  0.5 * ( tlpiChg + tlpiNeu ).M2() * tlNeutrino.Vect();
	thetaPol[i] = polarimeters[i].Angle( mcTauTLV[i].Vect() );	

	polarimeters[i]*=1./polarimeters[i].Mag();
	//	cout << "rho mode polarimeter: " << endl; 
	//	polvec.Print();
	//	tauola_hel_MC[i].Print();
	//	cout << "    angle to tauola " << polvec.Angle(	tauola_hel_MC[i] ) << endl;

	// approximate neutrino by colin approx

	TLorentzVector tlVis = tlpiChg+tlpiNeu;

	//TVector3 visTauDir( tlpiChg.Vect() + tlpiNeu.Vect() );
	//visTauDir.Scale( 1./visTauDir.Mag() );

	float colinNeuEn = ( pow(tauMass,2) - tlVis.M2() ) / ( 2*( tlVis.E() - tlVis.Vect().Mag() )  );

	TLorentzVector tlNeutrinoColinear;
	tlNeutrinoColinear.SetVectM( colinNeuEn*(tlVis.Vect().Unit()), 0 );

	TVector3 polvecColin =
          tauMass * (tlpiChg.E() - tlpiNeu.E() ) * ( tlpiChg.Vect() - tlpiNeu.Vect() ) +
          0.5 * ( tlpiChg + tlpiNeu ).M2() * tlNeutrinoColinear.Vect();

	thetaPolColin[i] = polvecColin.Angle( tlVis.Vect() );

      } else if ( decMode[i] == tauDec_a1chg_3p ) {

	// hmm not sure what to do here

      }
    }



    //    TVector3 polarimeters[2];
    TVector3 polarimetersTrans[2];
    for (int i=0; i<2; i++) {

      polarimetersTrans[i]=polarimeters[i];

      polarimetersTrans[i] -= refDir*(polarimetersTrans[i].Dot(refDir));

    }


    // from davier paper
    // I think this is not correct....
    float Omega1 = ( cos(thetaPol[0]) + cos(thetaPol[1]) ) / (1. + cos(thetaPol[0]) * cos(thetaPol[1]));
    float Omega2 = ( cos(thetaPol[0]) - cos(thetaPol[1]) ) / (1. - cos(thetaPol[0]) * cos(thetaPol[1]));

    int massrange = -1;
    if      ( ttmass > 450 ) massrange=1;
    else if ( ttmass > 81 && ttmass<101 ) massrange=0;
    else return;

    _h_mcTauTauMassSel[massrange]->Fill( ttmass );

    _h_mcTauEn[massrange]->Fill( tauMC[0]->getEnergy(), tauMC[1]->getEnergy() );
    _h_mcTauCosth[massrange]->Fill( mcTauMom[0].CosTheta(), mcTauMom[1].CosTheta() );

    if ( decMode[0]==tauDec_pi && decMode[1]==tauDec_pi ) {
      _h_mcPiPiDecAng[massrange]->Fill(getTauPiDecayAngle_MC( tauMC[0] ), getTauPiDecayAngle_MC( tauMC[1] ) );

      _h_mcPiPi_PiEnFrac[massrange]->Fill( getTauPiEnFrac_MC( tauMC[0] ), getTauPiEnFrac_MC( tauMC[1] ) );
      _h_mcPiPi_PiEn[massrange]->Fill( getTauPiEn_MC( tauMC[0] ), getTauPiEn_MC( tauMC[1] ) );

      _h_mcPiPi_polThTh[massrange]->Fill( cos(thetaPol[0]), cos(thetaPol[1]) );
      _h_mcPiPi_polComb1[massrange]->Fill( Omega1 );
      _h_mcPiPi_polComb2[massrange]->Fill( Omega2 );
    } else if ( decMode[0]==tauDec_pi && decMode[1]==tauDec_rho ) {
      _h_mcPiRho_polThTh[massrange]->Fill( cos(thetaPol[0]), cos(thetaPol[1]) );
      _h_mcPiRho_polComb1[massrange]->Fill( Omega1 );
      _h_mcPiRho_polComb2[massrange]->Fill( Omega2 );
    } else if ( decMode[0]==tauDec_rho && decMode[1]==tauDec_pi ) {
      _h_mcRhoPi_polThTh[massrange]->Fill( cos(thetaPol[0]), cos(thetaPol[1]) );
      _h_mcRhoPi_polComb1[massrange]->Fill( Omega1 );
      _h_mcRhoPi_polComb2[massrange]->Fill( Omega2 );
    } else if ( decMode[0]==tauDec_rho && decMode[1]==tauDec_rho ) {
      _h_mcRhoRho_polThTh[massrange]->Fill( cos(thetaPol[0]), cos(thetaPol[1]) );
      _h_mcRhoRho_polComb1[massrange]->Fill( Omega1 );
      _h_mcRhoRho_polComb2[massrange]->Fill( Omega2 );

      _h_mcRhoRho_colin_polThTh[massrange]->Fill( cos(thetaPolColin[0]), cos(thetaPolColin[1]) );
    }

    // now try to do a fit using MC information

    if ( 
	false && 
	( decMode[0]==tauDec_pi || decMode[0]==tauDec_rho ) && 
	( decMode[1]==tauDec_pi || decMode[1]==tauDec_rho ) ) {

      const MCParticle* mcchg[2];
      TLorentzVector neuTLV[2];
      bool mcokevent=true;

      for (int i=0; i<2; i++) {

	mcchg[i]=0;
	neuTLV[i].SetXYZT(0,0,0,0);
	int nch(0);
	for ( size_t j=0; j< stable_tau_desc[i].size(); j++) {
	  const MCParticle* mcp=stable_tau_desc[i][j];
	  int pdg=abs( mcp->getPDG() );

	  if ( fabs(mcp->getCharge())>0.5 ) {
	    mcchg[i]=stable_tau_desc[i][j];
	    nch++;
	  } else if (pdg!=12 && pdg!=14 && pdg!=16) {
	    neuTLV[i]+=TLorentzVector( mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2], mcp->getEnergy() );
	  }
	}

	if ( nch>1 ) mcokevent=false;
      }


      if ( mcokevent ) {

	//	cout << "TRY TO FIT!!" << endl;

	_eventFitter->reset();
	_eventFitter->setEcom(500);
	// _eventFitter->setIP( TVector3( tauMC[0]->getVertex()[0], tauMC[0]->getVertex()[1], tauMC[0]->getVertex()[2] ) );
	//	_eventFitter->setVerbose();

	for (int i=0; i<2; i++) {
	  _eventFitter->setNeutralHadronMomentum( i, neuTLV[i] );
	  _eventFitter->setChargedHadronMomentum( i, TLorentzVector( mcchg[i]->getMomentum()[0], mcchg[i]->getMomentum()[1], mcchg[i]->getMomentum()[2], mcchg[i]->getEnergy() ) );
	  // point on trajectory
	  _eventFitter->setChargedHadronImpactParameterVector( i , TVector3(mcchg[i]->getVertex()[0], mcchg[i]->getVertex()[1], mcchg[i]->getVertex()[2] ) );
	}

	int strategy=0;
	cout << "fitting..." << endl;
	_eventFitter->fitIt_single_single(strategy);
	cout << "fitting done" << endl;


//	// let's have a look at something...
//
//	int nb=200;
//	float psirange=0.2;
//
//	iplot++;
//	bool doplot=iplot<20;
//
//	TString hn = "MCFIT_pt_ev";
//	TString hnz = "MCFIT_pz_ev";
//	TString hnm = "MCFIT_isrM_ev";
//	TString hnE = "MCFIT_E_ev";
//
//	hn+=iplot;
//	hnE+=iplot;
//	hnm+=iplot;
//	hnz+=iplot;
//
//	if (doplot) {
//	  _fout->cd();
//	  for ( int iss=0; iss<4; iss++) {
//	    TString ll("_sol"); ll+=iss;
//	    _h_evnPt[iss] = new TH2F( hn+ll, hn+ll, nb, -psirange - psirange/nb , psirange +  psirange/nb,  nb, -psirange -  psirange/nb, psirange +  psirange/nb);
//	    _h_evnPz[iss] = new TH2F( hnz+ll, hnz+ll, nb, -psirange - psirange/nb , psirange +  psirange/nb,  nb, -psirange -  psirange/nb, psirange +  psirange/nb);
//	    _h_evnISRMass[iss] = new TH2F( hnm+ll, hnm+ll, nb, -psirange - psirange/nb , psirange +  psirange/nb,  nb, -psirange -  psirange/nb, psirange +  psirange/nb);
//	    _h_evnE[iss] = new TH2F( hnE+ll, hnE+ll, nb, -psirange - psirange/nb , psirange +  psirange/nb,  nb, -psirange -  psirange/nb, psirange +  psirange/nb);
//	  }
//
//	  for (int i1=0; i1<nb; i1++) {
//
//	    float psi1 = -psirange + i1*(psirange*2/(nb-1));
//
//	    for (int i2=0; i2<nb; i2++) {
//	      
//	      float psi2 = -psirange + i2*(psirange*2/(nb-1));
//
//	      _eventFitter->setPsi( 0, psi1 );
//	      _eventFitter->setPsi( 1, psi2 );
//
//	      // the 2x2 solutions
//	      for (int isol1=0; isol1<2; isol1++) {
//
//		// require positive decay length
//		if ( _eventFitter->getDecayLength( 0, isol1 ) < 0 ) continue;
//
//		for (int isol2=0; isol2<2; isol2++) {
//
//		  // require positive decay length
//		  if ( _eventFitter->getDecayLength( 1, isol2 ) < 0 ) continue;
//
//		  float ne1 = _eventFitter->getNeutrino4Momentum( 0, isol1 ).E();
//		  float ne2 = _eventFitter->getNeutrino4Momentum( 1, isol2 ).E();
//		
//		  if (ne1 > 0 && ne2 > 0 ) {
//
//		    //		  cout << i1 << " " << i2 << " " << psi1 << " " << psi2;
//		    //		  cout << " : SOLUTION " << isol1 << " " << isol2 << " neutrino energies: " << ne1 << ", " << ne2 << endl;
//		
//		    float tautau_pt = _eventFitter->getPt(isol1, isol2);
//
//		    TLorentzVector tautautlv = _eventFitter->getEvent4Momentum(isol1, isol2);
//
//		    //		  cout << " pt = "  << tautau_pt << " " << tautautlv.Pt();
//		    //		  tautautlv.Print();
//
//		    TLorentzVector isrTLV(0,0,0,500);
//		    isrTLV -= tautautlv;
//
//
//		    _h_evnPz[isol1+2*isol2] -> Fill( psi1, psi2, tautautlv.Pz() );
//		    _h_evnPt[isol1+2*isol2] -> Fill( psi1, psi2, tautau_pt);
//		    _h_evnE[isol1+2*isol2] -> Fill( psi1, psi2, tautautlv.E() );
//		    _h_evnISRMass[isol1+2*isol2] -> Fill( psi1, psi2, isrTLV.M() );
//		  
//		    _h_mcfitall_ttpz->Fill( tautautlv.Pz() );
//		    _h_mcfitall_ttE->Fill( tautautlv.E() );
//		    _h_mcfitall_ttpt->Fill( tautautlv.Pt() );
//		    
//		  }
//		}
//	      }
//	    }
//	  }
//	}



	//int strategy=2;
	//_eventFitter->fitIt_single_single(strategy);

	//	cout << "DONE FIT!!" << endl;

      } // OK event


    } // decay mode



    // now let's try to use PFOs...
    if ( 0 ) {


    LCRelationNavigator* relNavi2(0);
    try {
      LCCollection* linkcol = evt->getCollection( "MCTruthRecoLink" );
      relNavi2 = new LCRelationNavigator(linkcol);
    } catch(DataNotAvailableException &e) {};

    // first cheating the assiociation to MC particles
    if ( ( decMode[0]==tauDec_pi || decMode[0]==tauDec_rho ) || 
	 ( decMode[1]==tauDec_pi || decMode[1]==tauDec_rho ) ) {

      cout << "START MATCHING TO PFOS" << endl;

      std::vector < ReconstructedParticle*> matchedPfos[2];
      for (int i=0; i<2; i++) {

	cout << "TAU " << i << " dec mode " << decMode[i] << endl;


	for ( size_t j=0; j<stable_tau_desc[i].size(); j++) {
	  MCParticle* mcp = const_cast <MCParticle*> (stable_tau_desc[i][j]);

	  cout << "MC tau daghter: " << i << " " << j << " " << mcp->getPDG() << " " << mcp->getEnergy() << endl;

	  if ( abs( mcp->getPDG() ) ==12 ||  abs( mcp->getPDG() ) ==14 ||  abs( mcp->getPDG() ) == 16 ) continue;
	  float maxwt(-9999);
	  ReconstructedParticle* bestMatch(0);

	  cout << "size" << relNavi2->getRelatedToObjects( mcp ).size() << endl;

	  for (size_t k=0; k<relNavi2->getRelatedToObjects( mcp ).size(); k++) {
	    ReconstructedParticle* rp = (ReconstructedParticle*) relNavi2->getRelatedToObjects( mcp )[k];
	    float wt = relNavi2->getRelatedToWeights( mcp )[k];
	    if ( wt > maxwt ) {
	      maxwt=wt;
	      bestMatch=rp;
	    }
	  }
	  if ( bestMatch ) {

	    cout << "   best PFO match: " << bestMatch->getType() << " " << bestMatch->getEnergy() << " weight " << maxwt << endl;

	    if ( find (  matchedPfos[i].begin(),  matchedPfos[i].end(), bestMatch ) == matchedPfos[i].end() ) // don't add same PFO twice
	      matchedPfos[i].push_back(bestMatch);
	  } else {
	    cout << " no match " << endl;
	  }
	}

	// enum {pfoType_Photon=0, pfoType_ChHad, pfoType_NeuHad, pfoType_Mu, pfoType_El, pfoType_Lambda, pfoType_KShort, NPFOTYPE };
	int nPfoByType[NPFOTYPE]={0};
	int nch(0), nneu(0);

	cout << "PFOs matched to tau " << i << endl;
	for ( size_t j=0; j< matchedPfos[i].size(); j++) {
	  cout << " -- " << j << " " << matchedPfos[i][j]->getType() << " " << matchedPfos[i][j]->getEnergy() << endl;

	  switch ( abs( matchedPfos[i][j]->getType() ) ) {
	  case 22:
	    nPfoByType[pfoType_Photon]++;
	    nneu++;
	    break;
	  case 211:
	    nPfoByType[pfoType_ChHad]++;
	    nch++;
	    break;
	  case 2112:
	    nPfoByType[pfoType_NeuHad]++;
	    nneu++;
	    break;
	  case 13:
	    nPfoByType[pfoType_Mu]++;
	    nch++;
	    break;
	  case 11:
	    nPfoByType[pfoType_El]++;
	    nch++;
	    break;
	  case 3122:
	    nPfoByType[pfoType_Lambda]++;
	    nneu++;
	    break;
	  case 310:
	    nPfoByType[pfoType_KShort]++;
	    nneu++;
	    break;
	  default:
	    cout << "unknown PFO type: " << matchedPfos[i][j]->getType() << endl;
	  }

	  if ( fabs(matchedPfos[i][j]->getCharge())>0.5 ) {

	    float z0 = matchedPfos[i][j]->getTracks()[0]->getZ0();

	    _h_mcPfoMatch_chgz0->Fill( z0 );

	    _h_mcPfoMatch_chgz0WrtPv->Fill( z0 - tauMC[0]->getVertex()[2] );

	  }


	}

	
	if ( decMode[i] == tauDec_pi ) {
	  _h_mcPfoMatch_tauPi_nChNeu->Fill(nch, nneu);
	  _h_mcPfoMatch_tauPi_nChpi  ->Fill( nPfoByType[pfoType_ChHad] );
	  _h_mcPfoMatch_tauPi_nNeuh  ->Fill( nPfoByType[pfoType_NeuHad] );
	  _h_mcPfoMatch_tauPi_nGamma ->Fill( nPfoByType[pfoType_Photon] );
	  _h_mcPfoMatch_tauPi_nEl    ->Fill( nPfoByType[pfoType_El] );
	  _h_mcPfoMatch_tauPi_nMu    ->Fill( nPfoByType[pfoType_Mu] );
	  _h_mcPfoMatch_tauPi_nV0    ->Fill( nPfoByType[pfoType_Lambda]+nPfoByType[pfoType_KShort] );
	} else if ( decMode[i] == tauDec_rho ) {
	  _h_mcPfoMatch_tauRho_nChNeu->Fill(nch, nneu);
	  _h_mcPfoMatch_tauRho_nChpi  ->Fill( nPfoByType[pfoType_ChHad] );
	  _h_mcPfoMatch_tauRho_nNeuh  ->Fill( nPfoByType[pfoType_NeuHad] );
	  _h_mcPfoMatch_tauRho_nGamma ->Fill( nPfoByType[pfoType_Photon] );
	  _h_mcPfoMatch_tauRho_nEl    ->Fill( nPfoByType[pfoType_El] );
	  _h_mcPfoMatch_tauRho_nMu    ->Fill( nPfoByType[pfoType_Mu] );
	  _h_mcPfoMatch_tauRho_nV0    ->Fill( nPfoByType[pfoType_Lambda]+nPfoByType[pfoType_KShort] );
	}


      }


    }





    int nPfoByType[2][NTYPE];
    for (int i=0; i<2; i++) {
      for (int j=0; j<NTYPE; j++) {
	nPfoByType[i][j]=0;
      }
    }

    try {
      LCCollection* pfocol = evt->getCollection( "PandoraPFOs" );
      _h_nPFO[massrange]->Fill(pfocol->getNumberOfElements());

      pfototen=0;

      TLorentzVector tot4mom[2];
      for (int i=0; i<2; i++) tot4mom[i].SetXYZT(0,0,0,0);

      std::vector < TLorentzVector > photon4moms[2];

      for (int j=0; j<pfocol->getNumberOfElements(); j++) {
	ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (pfocol->getElementAt(j));

	pfototen+=pfo->getEnergy();

	int mytype = gettype(pfo->getType());
	_h_pfo_id[massrange]->Fill( mytype );
	TVector3 thismom( pfo->getMomentum()[0], pfo->getMomentum()[1], pfo->getMomentum()[2] );
	_h_pfotauang[massrange]->Fill( std::min( thismom.Angle(mcTauMom[0]), thismom.Angle(mcTauMom[1]) ) , std::max ( thismom.Angle(mcTauMom[0]), thismom.Angle(mcTauMom[1]) ) );

	bool closestIsMinus = thismom.Angle( mcTauMom[0] ) < thismom.Angle( mcTauMom[1] );
	int itau = int(!closestIsMinus);


	if ( thismom.Angle( mcTauMom[itau] ) > 0.5 ) continue; // reject pfos at large angle to tau (e.g. isr...)


	nPfoByType[itau][mytype]++;

	TLorentzVector fourmom( pfo->getMomentum()[0], pfo->getMomentum()[1], pfo->getMomentum()[2], pfo->getEnergy() );

	tot4mom[itau]+=fourmom;

	if ( pfo->getType()==22 ) photon4moms[itau].push_back(fourmom);

	_h_ntrk   [massrange]->Fill(pfo->getTracks().size());
	if ( pfo->getCharge() != 0 && pfo->getTracks().size()>0 ) {
	  Track* trk = pfo->getTracks()[0];
	  const TrackState * atip = trk->getTrackState(TrackState::AtIP);

	  _h_trkD0   [massrange]->Fill(atip->getD0());
	  _h_trkD0err[massrange]->Fill(sqrt(atip->getCovMatrix()[0]));
	  _h_trkD0sig[massrange]->Fill(atip->getD0()/sqrt(atip->getCovMatrix()[0]));
	}	

      }

      _h_pfototen_mumu   ->Fill(pfototen);
      _h_pfototen_tautau ->Fill(pfototen);


      for (int i=0; i<2; i++) {
	_h_visMass[massrange][decMode[i]]->Fill( tot4mom[i].M() );

	if ( photon4moms[i].size()>1 ) {
	  for (size_t j=0; j<photon4moms[i].size()-1; j++ ) {
	    for (size_t k=j+1; k<photon4moms[i].size(); k++) {
	      _h_ggMass[massrange][decMode[i]]->Fill( (photon4moms[i][j] + photon4moms[i][k] ).M() );
	    }
	  }
	}

      }


      for (int i=0; i<2; i++) {
	int ntot(0);
	for (int j=0; j<NTYPE; j++) {
	  ntot+=nPfoByType[i][j];
	  _h_pfos_by_type[massrange]->Fill( decMode[i], nPfoByType[i][j] );
	  _h_pfos_by_dec_type[massrange][decMode[i]][j]->Fill(nPfoByType[i][j] );
	}
	_h_npfos[massrange][decMode[i]]->Fill(ntot);
      }

    } catch(DataNotAvailableException &e) {};

  }


  }


  return;
}


void validateDST_TauTauProcessor::check( LCEvent * ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void validateDST_TauTauProcessor::end(){
  std::cout << "validateDST_TauTauProcessor::end()  " << name()
            << std::endl ;

  _fout->Write(0);
  _fout->Close();
}

