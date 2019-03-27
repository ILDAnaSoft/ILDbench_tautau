#include "danielKeitaTauFinder.h"
#include <iostream>
#include <map>
#include <iomanip>

#include <marlin/Global.h>
#include "lcio.h"

#include "EVENT/LCCollection.h"
#include "EVENT/Track.h"
#include "EVENT/Vertex.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/ReconstructedParticleImpl.h"

#include "TVector3.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include "vertexInfo.h"

using namespace lcio ;
using namespace marlin ;
using std::cout;
using std::endl;
using std::setw;

danielKeitaTauFinderProcessor adanielKeitaTauFinderProcessor ;


danielKeitaTauFinderProcessor::danielKeitaTauFinderProcessor() : Processor("danielKeitaTauFinderProcessor") {
  
  // modify processor description
  _description = "danielKeitaTauFinderProcessor does whatever it does ..." ;

  registerProcessorParameter("outputFilename",
                             "name of output file",
                             _outfile,
                             std::string( "findTau") );

  _conesize = 0.1; 
  _highestPt = false; // how to choose seed tracks: false = momentum, true = pt
  _useDistilled = false;  // use the distilled pfo collection?


  _2tauEvtFitter = new tautau2fEventFitter();

  _relNavi  = 0;
  _relNavi2 = 0;

}


void danielKeitaTauFinderProcessor::init() { 
  cout << "hello from danielKeitaTauFinderProcessor::init" << endl;

  _nNotconvertedInTrk=0;
  _nconvertedInTrk=0;

  //  decLab[TDEC_E  ] = "E";
  //  decLab[TDEC_M  ] = "M";
  //  decLab[TDEC_PI ] = "PI";
  //  decLab[TDEC_RHO] = "RHO";
  //  decLab[TDEC_A3P] = "A3P";
  //  decLab[TDEC_A1P] = "A1P";
  //  decLab[TDEC_HAD] = "HAD";
  //  decLab[TDEC_K  ] = "K";  

  TString fnn(_outfile.c_str()); fnn+=".root";
  _fout = new TFile(fnn,"recreate");
  h_gamgam_mass = new TH1F( "ggmass", "ggmass", 200, 0, 0.4 );

  h_ttmass = new TH1F( "mcttmass", "mcttmass", 205, 0, 510 );
  hSEL_ttmass = new TH1F( "SEL_mcttmass", "SEL_mcttmass", 205, 0, 510 );

  h_conversionPos = new TH2F("convPos","convPos", 500, 0, 3000, 500, 0, 2000);
  h_v0Pos = new TH2F("v0Pos","v0Pos", 500, 0, 3000, 500, 0, 2000);

  h_conversionPos2 = new TH2F("convPos2","convPos2", 500, 0, 2000, 500, 0, 400);
  h_v0Pos2 = new TH2F("v0Pos2","v0Pos2", 500, 0, 2000, 500, 0, 400);

  h_neuPandoraCompound_mass   =  new TH1F("neuPandoraCompound_mass",   "neuPandoraCompound_mass", 100, 0, 3 );
  h_neuPandoraCompound_ntrk   =  new TH1F("neuPandoraCompound_ntrk",   "neuPandoraCompound_ntrk", 5, 0, 5 );
  h_neuPandoraCompound_vtxPos =  new TH2F("neuPandoraCompound_vtxPos", "neuPandoraCompound_vtxPos",  500, 0, 3000, 500, 0, 2000);
  h_neuPandoraCompound_vtxPos2 =  new TH2F("neuPandoraCompound_vtxPos2", "neuPandoraCompound_vtxPos2",  500, 0, 2000, 500, 0, 400);
  h_neuPandoraCompound_vtxChisq  =  new TH1F("neuPandoraCompound_vtxChisq",   "neuPandoraCompound_vtxChisq", 100, 0, 10 );


  h_mctautau_ecom_tauMinusCosth = new TH2F( "mctautau_ecom_tauMinusCosth",  "mctautau_ecom_tauMinusCosth", 100,0,500, 100, -1, 1 );


  h_gammaClus_nMajor           = new  TH1F("gammaClus_nMajor",        "gammaClus_nMajor",        10,-0.5, 9.5 );
  h_gammaClus_nSignf	       = new  TH1F("gammaClus_nSignf",	      "gammaClus_nSignf",        10,-0.5, 9.5 );
  h_gammaClus_nSignf_bMass     = new  TH2F("gammaClus_nSignf_bMass",  "gammaClus_nSignf_bMass",  10,-0.5, 9.5,100,0,1 );
  h_gammaClus_nSignf_bDist     = new  TH2F("gammaClus_nSignf_bDist",  "gammaClus_nSignf_bDist",  10,-0.5, 9.5,100,0,100 );
  h_neutronClus_nMajor	       = new  TH1F("neutronClus_nMajor",      "neutronClus_nMajor",      10,-0.5, 9.5 );
  h_neutronClus_nSignf	       = new  TH1F("neutronClus_nSignf",      "neutronClus_nSignf",      10,-0.5, 9.5 );
  h_neutronClus_nSignf_bMass   = new  TH2F("neutronClus_nSignf_bMass","neutronClus_nSignf_bMass",10,-0.5, 9.5,100,0,1 );
  h_neutronClus_nSignf_bDist   = new  TH2F("neutronClus_nSignf_bDist","neutronClus_nSignf_bDist",10,-0.5, 9.5,100,0,100 );

  h_neutronClus_eOnP_beforeAfter = new TH2F( "neutronClus_eOnP_beforeAfter",  "neutronClus_eOnP_beforeAfter", 100, 0, 10, 100, 0, 10 );

  h_neutralHadronPFO_mainMCcontrib_pdgEn = new TH2F("neutralHadronPFO_mainMCcontrib_pdgEn", "neutralHadronPFO_mainMCcontrib_pdgEn", 100, 0, 100, 10, -0.5, 9.5 );
  h_neutralHadronPFO_mainMCcontrib_pdgEn->GetYaxis()->SetBinLabel( 1, "undef" );
  h_neutralHadronPFO_mainMCcontrib_pdgEn->GetYaxis()->SetBinLabel( 2, "chHad" );
  h_neutralHadronPFO_mainMCcontrib_pdgEn->GetYaxis()->SetBinLabel( 3, "photon" );
  h_neutralHadronPFO_mainMCcontrib_pdgEn->GetYaxis()->SetBinLabel( 4, "electron" );
  h_neutralHadronPFO_mainMCcontrib_pdgEn->GetYaxis()->SetBinLabel( 5, "muon" );
  h_neutralHadronPFO_mainMCcontrib_pdgEn->GetYaxis()->SetBinLabel( 6, "neuHad" );

  h_rhoDecaySinglePhoClus_reason = new TH1F( "rhoDecaySinglePhoClus_reason", "rhoDecaySinglePhoClus_reason", 10, -0.5, 9.5);
  h_rhoDecaySinglePhoClus_reason->GetXaxis()->SetBinLabel( 1, "converted" );
  h_rhoDecaySinglePhoClus_reason->GetXaxis()->SetBinLabel( 2, "noPFO (lowen)" );
  h_rhoDecaySinglePhoClus_reason->GetXaxis()->SetBinLabel( 3, "noPFO (other)" );
  h_rhoDecaySinglePhoClus_reason->GetXaxis()->SetBinLabel( 4, "merged(phoClus)" );
  h_rhoDecaySinglePhoClus_reason->GetXaxis()->SetBinLabel( 5, "merged(nhadClus)" );
  h_rhoDecaySinglePhoClus_reason->GetXaxis()->SetBinLabel( 6, "merged(chgClus)" );
  h_rhoDecaySinglePhoClus_reason->GetXaxis()->SetBinLabel( 7, "attachedToChg" );
  h_rhoDecaySinglePhoClus_reason->GetXaxis()->SetBinLabel( 10, "other" );

  h_rhoDecaySinglePhoClus_mergedGammaClusterEvals = new TH2F( "rhoDecaySinglePhoClus_mergedGammaClusterEvals", "rhoDecaySinglePhoClus_mergedGammaClusterEvals", 100,0,20,100,0,20);
  h_rhoDecaySinglePhoClus_singleGammaClusterEvals = new TH2F( "rhoDecaySinglePhoClus_singleGammaClusterEvals", "rhoDecaySinglePhoClus_singleGammaClusterEvals", 100,0,20,100,0,20);


  for (int iss=0; iss<NCLASS; iss++) {
    TString samp="sample";
    samp+=iss; samp+="_";

    h_nIsoMuons[iss] = new TH1F(  samp+ "nIsoMuons",  samp+ "nIsoMuons", 10, -0.5, 9.5 );
    h_nIsoElectrons[iss] = new TH1F(  samp+ "nIsoElectrons",  samp+ "nIsoElectrons", 10, -0.5, 9.5 );

    h_pi_mcPolarExact[iss] = new TH1F( samp+ "pi_mcPolarExact",  samp+ "pi_mcPolarExact", 110, -1.1, 1.1 );
    h_rho_mcPolarExact[iss] = new TH1F( samp+ "rho_mcPolarExact",  samp+ "rho_mcPolarExact", 110, -1.1, 1.1);
    h_rho_mcPolarExactApprox[iss] = new TH2F( samp+ "rho_mcPolarExactApprox",  samp+ "rho_mcPolarExactApprox",110, -1.1, 1.1, 110, -1.1, 1.1);

//    h_rho_mcapprox_cospsi_helPos[iss]   = new TH1F( samp+ "rho_mcapprox_cospsi_helPos",  samp+ "rho_mcapprox_cospsi_helPos", 100, -1.5, 1.5 );
//    h_rho_mcapprox_costheta_helPos[iss] = new TH1F( samp+ "rho_mcapprox_costheta_helPos",  samp+ "rho_mcapprox_costheta_helPos", 100, -1.5, 1.5 );
//    h_rho_mcapprox_cosbeta_helPos[iss]  = new TH1F( samp+ "rho_mcapprox_cosbeta_helPos",  samp+ "rho_mcapprox_cosbeta_helPos", 100, -1.5, 1.5 );
//
//    h_rho_mcapprox_cospsi_helNeg[iss]   = new TH1F( samp+ "rho_mcapprox_cospsi_helNeg",  samp+ "rho_mcapprox_cospsi_helNeg", 100, -1.5, 1.5 );
//    h_rho_mcapprox_costheta_helNeg[iss] = new TH1F( samp+ "rho_mcapprox_costheta_helNeg",  samp+ "rho_mcapprox_costheta_helNeg", 100, -1.5, 1.5 );
//    h_rho_mcapprox_cosbeta_helNeg[iss]  = new TH1F( samp+ "rho_mcapprox_cosbeta_helNeg",  samp+ "rho_mcapprox_cosbeta_helNeg", 100, -1.5, 1.5 );

    h_pi_mcPolarExact_helNeg[iss]  = new TH1F( samp+  "pi_mcPolarExact_helNeg",  samp+  "pi_mcPolarExact_helNeg", 110, -1.1, 1.1 );
    h_pi_mcPolarExact_helPos[iss]  = new TH1F( samp+  "pi_mcPolarExact_helPos",  samp+  "pi_mcPolarExact_helPos", 110, -1.1, 1.1 );
    h_rho_mcPolarExact_helNeg[iss] = new TH1F( samp+ "rho_mcPolarExact_helNeg",  samp+ "rho_mcPolarExact_helNeg", 110, -1.1, 1.1);
    h_rho_mcPolarExact_helPos[iss] = new TH1F( samp+ "rho_mcPolarExact_helPos",  samp+ "rho_mcPolarExact_helPos", 110, -1.1, 1.1);


    hSEL_pi_mcPolarExact_helNeg[iss]  = new TH1F( samp+ "SEL_pi_mcPolarExact_helNeg",  samp+  "SEL_pi_mcPolarExact_helNeg", 110, -1.1, 1.1 );
    hSEL_pi_mcPolarExact_helPos[iss]  = new TH1F( samp+ "SEL_pi_mcPolarExact_helPos",  samp+  "SEL_pi_mcPolarExact_helPos", 110, -1.1, 1.1 );
    hSEL_rho_mcPolarExact_helNeg[iss] = new TH1F( samp+ "SEL_rho_mcPolarExact_helNeg",  samp+ "SEL_rho_mcPolarExact_helNeg", 110, -1.1, 1.1);
    hSEL_rho_mcPolarExact_helPos[iss] = new TH1F( samp+ "SEL_rho_mcPolarExact_helPos",  samp+ "SEL_rho_mcPolarExact_helPos", 110, -1.1, 1.1);



    h_pi_mcPolarApproxCheatEn[iss] = new TH1F( samp+ "pi_mcPolarApproxCheatEn",  samp+ "pi_mcPolarApproxCheatEn", 110, -1.1, 1.1 );
    h_rho_mcPolarApproxCheatEn[iss] = new TH1F( samp+ "rho_mcPolarApproxCheatEn",  samp+ "rho_mcPolarApproxCheatEn", 110, -1.1, 1.1);

    h_pi_mcPolarApproxCheatEn_helNeg[iss]  = new TH1F( samp+  "pi_mcPolarApproxCheatEn_helNeg",  samp+  "pi_mcPolarApproxCheatEn_helNeg", 110, -1.1, 1.1 );
    h_pi_mcPolarApproxCheatEn_helPos[iss]  = new TH1F( samp+  "pi_mcPolarApproxCheatEn_helPos",  samp+  "pi_mcPolarApproxCheatEn_helPos", 110, -1.1, 1.1 );
    h_rho_mcPolarApproxCheatEn_helNeg[iss] = new TH1F( samp+ "rho_mcPolarApproxCheatEn_helNeg",  samp+ "rho_mcPolarApproxCheatEn_helNeg", 110, -1.1, 1.1);
    h_rho_mcPolarApproxCheatEn_helPos[iss] = new TH1F( samp+ "rho_mcPolarApproxCheatEn_helPos",  samp+ "rho_mcPolarApproxCheatEn_helPos", 110, -1.1, 1.1);



    h_pi_mcPolarApprox[iss] = new TH1F( samp+ "pi_mcPolarApprox",  samp+ "pi_mcPolarApprox", 110, -1.1, 1.1 );
    h_rho_mcPolarApprox[iss] = new TH1F( samp+ "rho_mcPolarApprox",  samp+ "rho_mcPolarApprox", 110, -1.1, 1.1);

    h_pi_mcPolarApprox_helNeg[iss]  = new TH1F( samp+  "pi_mcPolarApprox_helNeg",  samp+  "pi_mcPolarApprox_helNeg", 110, -1.1, 1.1 );
    h_pi_mcPolarApprox_helPos[iss]  = new TH1F( samp+  "pi_mcPolarApprox_helPos",  samp+  "pi_mcPolarApprox_helPos", 110, -1.1, 1.1 );
    h_rho_mcPolarApprox_helNeg[iss] = new TH1F( samp+ "rho_mcPolarApprox_helNeg",  samp+ "rho_mcPolarApprox_helNeg", 110, -1.1, 1.1);
    h_rho_mcPolarApprox_helPos[iss] = new TH1F( samp+ "rho_mcPolarApprox_helPos",  samp+ "rho_mcPolarApprox_helPos", 110, -1.1, 1.1);


    hSEL_pi_mcPolarApprox_helNeg[iss]  = new TH1F( samp+ "SEL_pi_mcPolarApprox_helNeg",  samp+ "SEL_pi_mcPolarApprox_helNeg", 110, -1.1, 1.1 );
    hSEL_pi_mcPolarApprox_helPos[iss]  = new TH1F( samp+ "SEL_pi_mcPolarApprox_helPos",  samp+ "SEL_pi_mcPolarApprox_helPos", 110, -1.1, 1.1 );
    hSEL_rho_mcPolarApprox_helNeg[iss] = new TH1F( samp+ "SEL_rho_mcPolarApprox_helNeg", samp+ "SEL_rho_mcPolarApprox_helNeg", 110, -1.1, 1.1);
    hSEL_rho_mcPolarApprox_helPos[iss] = new TH1F( samp+ "SEL_rho_mcPolarApprox_helPos", samp+ "SEL_rho_mcPolarApprox_helPos", 110, -1.1, 1.1);




    h_pirho_tauMinusCosth_mcPolarExact[iss] = new TH2F( samp+ "pirho_tauMinusCosth_mcPolarExact", samp+ "pirho_tauMinusCosth_mcPolarExact", 22, -1.1, 1.1, 22, -1.10, 1.10 );

    h_pirho_mcPolarExactPlusMinus[iss] = new TH2F( samp+ "pirho_mcPolarExactPlusMinus", samp+ "pirho_mcPolarExactPlusMinus", 22, -1.1, 1.1, 22, -1.10, 1.10 );
    h_pirho_mcPolarApproxPlusMinus[iss] = new TH2F( samp+ "pirho_mcPolarApproxPlusMinus", samp+ "pirho_mcPolarApproxPlusMinus", 22, -1.1, 1.1, 22, -1.10, 1.10 );

    h_pi_mcPolarComp[iss] = new TH2F( samp+ "pi_mcPolarComp",  samp+ "pi_mcPolarComp", 110, -1.10, 1.10 , 110, -1.10, 1.10 );
    h_rho_mcPolarComp[iss] = new TH2F( samp+ "rho_mcPolarComp",  samp+ "rho_mcPolarComp", 110, -1.10, 1.10 , 110, -1.10, 1.10);

    h_pi_mcPolar[iss] = new TH1F( samp+ "pi_mcPolar",  samp+ "pi_mcPolar", 110, 0., 1.10 );

    h_rho_mcPolar[iss] = new TH2F( samp+ "rho_mcPolar",  samp+ "rho_mcPolar", 110, 0., 1.1,  110, 0, 1.1);
    h_rho_mcPolar2[iss] = new TH2F( samp+ "rho_mcPolar2",  samp+ "rho_mcPolar2", 110, 0., 1.1, 110, 0., 1.1 );

    h_mc_tauMinus_costh[iss] = new TH1F( samp+ "mc_tauMinus_costh",  samp+ "mc_tauMinus_costh", 110, -1.1, 1.1 );
    hSEL_mc_tauMinus_costh[iss] = new TH1F( samp+ "SELmc_tauMinus_costh",  samp+ "SELmc_tauMinus_costh", 110, -1.1, 1.1 );

    h_mc_tauSpin[iss] = new TH2F( samp+ "mc_tauSpin", samp+ "mc_tauSpin", 3, -1.5, 1.5, 3,  -1.5, 1.5 );

    h_npfo_chg_nneu[iss] = new TH2F( samp+ "npfo_chg_nneu", samp+ "npfo_chg_nneu", 20, -0.5, 19.5, 20, -0.5, 19.5);

    h_rawmass_trimmass[iss] = new TH2F( samp+ "rawmass_trimmass", samp+  "rawmass_trimmass", 100, 0., 3., 100, 0., 3. );

    h_visjetAcolinearity[iss]    = new TH1F( samp+ "visjetAcolinearity",   samp+"visjetAcolinearity"  ,  100, 0, 0.25);
    h_visjetAcoplanarity[iss]    = new TH1F( samp+ "visjetAcoplanarity",   samp+"visjetAcoplanarity"  ,  100, 0, 0.25);

    h_prongAngle[iss]    = new TH1F( samp+ "prongAngle", samp+"prongAngle", 300,0,3.2);
    h_outsideEnergy[iss] = new TH1F( samp+ "outsideEnergy", samp+"outsideEnergy",100,0,100);
    h_insideEnergy[iss]  = new TH1F( samp+ "insideEnergy", samp+"insideEnergy",200,0,600);
    h_outsidePt[iss]     = new TH1F( samp+ "outsidePt", samp+"outsidePt",100,0,100);

    h_tchg_tmass[iss] = new TH2F( samp+ "tchg_tmass", samp+"tchg_tmass", 5,-2.5,2.5,100,0,5 );
    h_ngam_tmass[iss] = new TH2F( samp+ "ngam_tmass", samp+"ngam_tmass", 10,-0.5,9.5,100,0,5 );
    h_nchg_tmass[iss] = new TH2F( samp+ "nchg_tmass", samp+"nchg_tmass", 10,-0.5,9.5,100,0,5 );
    h_nnhad_tmass[iss] = new TH2F( samp+ "nnhad_tmass", samp+"nnhad_tmass", 10,-0.5,9.5,100,0,5 );
    h_ngam_nchg[iss] = new TH2F( samp+ "ngam_nchg", samp+"ngam_nchg", 10,-0.5,9.5,10,-0.5,9.5 );

    h_mcTau_costh[iss]  = new TH2F(  samp+ "MCtau_costh", samp+ "MCtau_costh", 100, 0, 1, 100, 0, 1 );

    hSEL_tchg_tmass[iss]           = new TH2F( samp+ "SEL_tchg_tmass", samp+"SEL_tchg_tmass", 5,-2.5,2.5,100,0,5 );

    hSEL_visjetAcolinearity[iss]    = new TH1F( samp+ "SEL_visjetAcolinearity",   samp+"SEL_visjetAcolinearity"  ,  100, 0, 0.25);
    hSEL_visjetAcoplanarity[iss]    = new TH1F( samp+ "SEL_visjetAcoplanarity",   samp+"SEL_visjetAcoplanarity"  ,  100, 0, 0.25);

    hSEL_prongAngle[iss]    = new TH1F( samp+ "SEL_prongAngle",   samp+"SEL_prongAngle"  , 300,0,3.2);
    hSEL_outsideEnergy[iss] = new TH1F( samp+ "SEL_outsideEnergy",samp+"SEL_outsideEnergy",100,0,100);
    hSEL_insideEnergy[iss]  = new TH1F( samp+ "SEL_insideEnergy", samp+"SEL_insideEnergy" ,600,0,600);
    hSEL_outsidePt[iss]     = new TH1F( samp+ "SEL_outsidePt",    samp+"SEL_outsidePt"    ,100,0,100);

    hSEL_ngam_tmass[iss]           = new TH2F( samp+ "SEL_ngam_tmass",          samp+"SEL_ngam_tmass", 10,-0.5,9.5,100,0,5 );
    hSEL_nchg_tmass[iss]           = new TH2F( samp+ "SEL_nchg_tmass",          samp+"SEL_nchg_tmass", 10,-0.5,9.5,100,0,5 );
    hSEL_nnhad_tmass[iss]          = new TH2F( samp+ "SEL_nnhad_tmass",         samp+"SEL_nnhad_tmass", 10,-0.5,9.5,100,0,5 );
    hSEL_ngam_nchg[iss]            = new TH2F( samp+ "SEL_ngam_nchg",           samp+"SEL_ngam_nchg", 10,-0.5,9.5,10,-0.5,9.5 );

    hSEL_mcTau_costh[iss]  = new TH2F(  samp+ "SEL_MCtau_costh", samp+ "SEL_MCtau_costh", 100, 0, 1, 100, 0, 1 );

    h_dec_tmass[iss]= new TH2F( samp+"h_dec_tmass", samp+"h_dec_tmass",tauUtils::NDECAYS,-.5,tauUtils::NDECAYS-0.5,100,0,5);
    h_dec_ngam[iss] = new TH2F( samp+"h_dec_ngam" , samp+"h_dec_ngam" ,tauUtils::NDECAYS,-.5,tauUtils::NDECAYS-0.5,10,-0.5,9.5);
    h_dec_nchg[iss] = new TH2F( samp+"h_dec_nchg" , samp+"h_dec_nchg" ,tauUtils::NDECAYS,-.5,tauUtils::NDECAYS-0.5,10,-0.5,9.5);
    h_dec[iss]      = new TH1F( samp+"h_dec"      , samp+"h_dec"      ,tauUtils::NDECAYS,-.5,tauUtils::NDECAYS-0.5);
    hSEL_dec[iss]   = new TH1F( samp+"hSEL_dec"   , samp+"hSEL_dec"   ,tauUtils::NDECAYS,-.5,tauUtils::NDECAYS-0.5);

    
    for (int i=0; i<tauUtils::NDECAYS; i++) {
      h_dec_tmass[iss]->GetXaxis()->SetBinLabel( i+1, tauUtils::getTauDecLab(i) );
      h_dec_ngam[iss]->GetXaxis()->SetBinLabel( i+1, tauUtils::getTauDecLab(i) );
      h_dec_nchg[iss]->GetXaxis()->SetBinLabel( i+1, tauUtils::getTauDecLab(i) );
      h_dec[iss]->GetXaxis()->SetBinLabel( i+1, tauUtils::getTauDecLab(i) );
      hSEL_dec[iss]->GetXaxis()->SetBinLabel( i+1, tauUtils::getTauDecLab(i) );
    }

    //h_a3p_cone_chgMC_pfo[iss] = new TH2F( samp+"a3p_cone_chgMC_pfo", samp+"a3p_cone_chgMC_pfo", 8,-0.5,7.5,8,-0.5,7.5 );
    //h_a3p_cone_chgMC_trk[iss] = new TH2F( samp+"a3p_cone_chgMC_trk", samp+"a3p_cone_chgMC_trk", 8,-0.5,7.5,8,-0.5,7.5 );
    //h_a3p_cone_trk_pfo  [iss] = new TH2F( samp+"a3p_cone_trk_pfo",   samp+"a3p_cone_trk_pfo",   8,-0.5,7.5,8,-0.5,7.5 );

    for (int idec=0; idec<NDEC; idec++) {

      TString lab=tauUtils::getTauDecLab(idec);
      if      (idec==NDEC-2) lab="OTHER";
      else if (idec==NDEC-1)   lab="ALL";

      h_dec_seed_clusterWidth1[iss][idec] = new TH1F( samp+"_"+lab+"_seed_logclusterWidth1", samp+"_"+lab+"_seed_logclusterWidth1", 100,-1,7);
      h_dec_seed_clusterWidth2[iss][idec] = new TH1F( samp+"_"+lab+"_seed_logclusterWidth2", samp+"_"+lab+"_seed_logclusterWidth2", 100,-1,7);
      h_dec_seed_clusterLength[iss][idec] = new TH1F( samp+"_"+lab+"_seed_logclusterLength", samp+"_"+lab+"_seed_logclusterLength", 100,-1,7);

      h_dec_seed_energy[iss][idec] = new TH1F( samp+"_"+lab+"_seed_energy", samp+"_"+lab+"_seed_energy", 100,0,300);


      h_dec_cone_ntracks    [iss][idec] = new TH1F( samp+"_"+lab+"_cone_ntracks",      samp+"_"+lab+"_cone_ntracks", 8,-0.5,7.5 );
      h_dec_cone_trackNTPC  [iss][idec] = new TH1F( samp+"_"+lab+"_cone_trackNTPC",    samp+"_"+lab+"_cone_trackNTPC", 50,0,250 );
      h_dec_cone_nchpfo     [iss][idec] = new TH1F( samp+"_"+lab+"_cone_nchpfo",      samp+"_"+lab+"_cone_nchpfo", 8,-0.5,7.5 );
      h_dec_cone_nnhadpfo   [iss][idec] = new TH1F( samp+"_"+lab+"_cone_nnhadpfo",    samp+"_"+lab+"_cone_nnhadpfo", 8,-0.5,7.5 );
      h_dec_cone_ngammapfo  [iss][idec] = new TH1F( samp+"_"+lab+"_cone_ngammapfo",   samp+"_"+lab+"_cone_ngammapfo", 8,-0.5,7.5 );
      h_dec_cone_npi0pfo  [iss][idec] = new TH1F( samp+"_"+lab+"_cone_npi0pfo",   samp+"_"+lab+"_cone_npi0pfo", 8,-0.5,7.5 );
      h_dec_cone_ncompoundpfo[iss][idec] = new TH1F( samp+"_"+lab+"_cone_ncompoundpfo",   samp+"_"+lab+"_cone_ncompoundpfo", 8,-0.5,7.5 );
      h_dec_cone_nhaden     [iss][idec] = new TH1F( samp+"_"+lab+"_cone_nhaden",      samp+"_"+lab+"_cone_nhaden",  100,0.,100 );
      h_dec_cone_gammaen    [iss][idec] = new TH1F( samp+"_"+lab+"_cone_gammaen",     samp+"_"+lab+"_cone_gammaen", 100,0.,300 );
      h_dec_cone_nhadenFrac [iss][idec] = new TH1F( samp+"_"+lab+"_cone_nhadenFrac",  samp+"_"+lab+"_cone_nhadenFrac", 100,0.,1. );
      h_dec_cone_gammaenFrac[iss][idec] = new TH1F( samp+"_"+lab+"_cone_gammaenFrac", samp+"_"+lab+"_cone_gammaenFrac", 100,0.,1. );



      h_dec_cone_visMass    [iss][idec] = new TH1F( samp+"_"+lab+"_cone_visMass",     samp+"_"+lab+"_cone_visMass", 100, 0, 3 );
      h_dec_cone_neutralvisMass[iss][idec] = new TH1F( samp+"_"+lab+"_cone_neutralvisMass", samp+"_"+lab+"_cone_neutralvisMass", 100, 0, 3 );
      h_dec_cone_visMassDiff    [iss][idec] = new TH1F( samp+"_"+lab+"_cone_visMassDiff",     samp+"_"+lab+"_cone_visMassDiff", 100, -1,1 );
      h_dec_cone_neutralvisMassDiff[iss][idec] = new TH1F( samp+"_"+lab+"_cone_neutralvisMassDiff", samp+"_"+lab+"_cone_neutralvisMassDiff", 100, -1,1 );

      h_dec_cone_vis_neutral_Mass[iss][idec] = new TH2F( samp+"_"+lab+"_cone_vis_neutral_Mass",   samp+"_"+lab+"_cone_vis_neutral_Mass", 100, 0, 3, 100, 0, 3 );
      h_dec_cone_vis_neutral_Mass_0g[iss][idec] = new TH2F( samp+"_"+lab+"_cone_vis_neutral_Mass_0g",   samp+"_"+lab+"_cone_vis_neutral_Mass_0g", 100, 0, 3, 100, 0, 3 );
      h_dec_cone_vis_neutral_Mass_1g[iss][idec] = new TH2F( samp+"_"+lab+"_cone_vis_neutral_Mass_1g",   samp+"_"+lab+"_cone_vis_neutral_Mass_1g", 100, 0, 3, 100, 0, 3 );
      h_dec_cone_vis_neutral_Mass_2g[iss][idec] = new TH2F( samp+"_"+lab+"_cone_vis_neutral_Mass_2g",   samp+"_"+lab+"_cone_vis_neutral_Mass_2g", 100, 0, 3, 100, 0, 3 );
      h_dec_cone_vis_neutral_Mass_3g[iss][idec] = new TH2F( samp+"_"+lab+"_cone_vis_neutral_Mass_3g",   samp+"_"+lab+"_cone_vis_neutral_Mass_3g", 100, 0, 3, 100, 0, 3 );
      h_dec_cone_vis_neutral_Mass_4g[iss][idec] = new TH2F( samp+"_"+lab+"_cone_vis_neutral_Mass_4g",   samp+"_"+lab+"_cone_vis_neutral_Mass_4g", 100, 0, 3, 100, 0, 3 );

      h_dec_cone_visEnergyDiff    [iss][idec] = new TH1F( samp+"_"+lab+"_cone_visEnergyDiff",     samp+"_"+lab+"_cone_visEnergyDiff", 100, -50,50 );
      h_dec_cone_neutralvisEnergyDiff[iss][idec] = new TH1F( samp+"_"+lab+"_cone_neutralvisEnergyDiff", samp+"_"+lab+"_cone_neutralvisEnergyDiff", 100, -50,50 );

      // h_dec_coneTRIM_visMass    [iss][idec] = new TH1F( samp+"_"+lab+"_coneTRIM_visMass",     samp+"_"+lab+"_coneTRIM_visMass", 100, 0, 2 );
      // h_dec_coneTRIM_neutralvisMass[iss][idec] = new TH1F( samp+"_"+lab+"_coneTRIM_neutralvisMass", samp+"_"+lab+"_coneTRIM_neutralvisMass", 100, 0, 2 );
      // h_dec_coneTRIM_vis_neutral_Mass[iss][idec] = new TH2F( samp+"_"+lab+"_coneTRIM_vis_neutral_Mass",   samp+"_"+lab+"_coneTRIM_vis_neutral_Mass", 100, 0, 2, 100, 0, 2 );
      // h_dec_coneTRIM_vis_neutral_Mass_0g[iss][idec] = new TH2F( samp+"_"+lab+"_coneTRIM_vis_neutral_Mass_0g",   samp+"_"+lab+"_coneTRIM_vis_neutral_Mass_0g", 100, 0, 2, 100, 0, 2 );
      // h_dec_coneTRIM_vis_neutral_Mass_1g[iss][idec] = new TH2F( samp+"_"+lab+"_coneTRIM_vis_neutral_Mass_1g",   samp+"_"+lab+"_coneTRIM_vis_neutral_Mass_1g", 100, 0, 2, 100, 0, 2 );
      // h_dec_coneTRIM_vis_neutral_Mass_2g[iss][idec] = new TH2F( samp+"_"+lab+"_coneTRIM_vis_neutral_Mass_2g",   samp+"_"+lab+"_coneTRIM_vis_neutral_Mass_2g", 100, 0, 2, 100, 0, 2 );
      // h_dec_coneTRIM_vis_neutral_Mass_3g[iss][idec] = new TH2F( samp+"_"+lab+"_coneTRIM_vis_neutral_Mass_3g",   samp+"_"+lab+"_coneTRIM_vis_neutral_Mass_3g", 100, 0, 2, 100, 0, 2 );
      // h_dec_coneTRIM_vis_neutral_Mass_4g[iss][idec] = new TH2F( samp+"_"+lab+"_coneTRIM_vis_neutral_Mass_4g",   samp+"_"+lab+"_coneTRIM_vis_neutral_Mass_4g", 100, 0, 2, 100, 0, 2 );
      // h_dec_coneTRIM_visMassDiff    [iss][idec] = new TH1F( samp+"_"+lab+"_coneTRIM_visMassDiff",     samp+"_"+lab+"_coneTRIM_visMassDiff", 100, -1,1 );
      // h_dec_coneTRIM_neutralvisMassDiff[iss][idec] = new TH1F( samp+"_"+lab+"_coneTRIM_neutralvisMassDiff", samp+"_"+lab+"_coneTRIM_neutralvisMassDiff", 100, -1,1 );
      // h_dec_coneTRIM_visEnergyDiff    [iss][idec] = new TH1F( samp+"_"+lab+"_coneTRIM_visEnergyDiff",     samp+"_"+lab+"_coneTRIM_visEnergyDiff", 100, -50,50 );
      // h_dec_coneTRIM_neutralvisEnergyDiff[iss][idec] = new TH1F( samp+"_"+lab+"_coneTRIM_neutralvisEnergyDiff", samp+"_"+lab+"_coneTRIM_neutralvisEnergyDiff", 100, -50,50 );

      h_dec_mc_ngamma               [iss][idec] = new TH1F( samp+"_"+lab+"_mc_ngamma"             , samp+"_"+lab+"_mc_ngamma"             , 10, -0.5, 9.5 ); 
      h_dec_mc_vismass		    [iss][idec] = new TH1F( samp+"_"+lab+"_mc_vismass"            , samp+"_"+lab+"_mc_vismass"            , 100, 0, 3 );
      h_dec_mc_visneutralmass	    [iss][idec] = new TH1F( samp+"_"+lab+"_mc_visneutralmass"     , samp+"_"+lab+"_mc_visneutralmass"     , 100, 0, 3 );
      h_dec_mc_cone_ngamma	    [iss][idec] = new TH1F( samp+"_"+lab+"_mc_cone_ngamma"        , samp+"_"+lab+"_mc_cone_ngamma"        , 10, -0.5, 9.5 );
      h_dec_mc_cone_vismass	    [iss][idec] = new TH1F( samp+"_"+lab+"_mc_cone_vismass"       , samp+"_"+lab+"_mc_cone_vismass"       , 100, 0, 3 );
      h_dec_mc_cone_visneutralmass  [iss][idec] = new TH1F( samp+"_"+lab+"_mc_cone_visneutralmass", samp+"_"+lab+"_mc_cone_visneutralmass", 100, 0, 3 );

      h_dec_mcall_cone_ngamma	    [iss][idec] = new TH1F( samp+"_"+lab+"_mcall_cone_ngamma"        , samp+"_"+lab+"_mcall_cone_ngamma"        , 10, -0.5, 9.5 );
      h_dec_mcall_cone_vismass	    [iss][idec] = new TH1F( samp+"_"+lab+"_mcall_cone_vismass"       , samp+"_"+lab+"_mcall_cone_vismass"       , 100, 0, 3 );
      h_dec_mcall_cone_visneutralmass  [iss][idec] = new TH1F( samp+"_"+lab+"_mcall_cone_visneutralmass", samp+"_"+lab+"_mcall_cone_visneutralmass", 100, 0, 3 );


      hSEL_dec_seed_clusterWidth1[iss][idec]  = new TH1F( samp+"_"+lab+"_SELseed_logclusterWidth1",  samp+"_"+lab+"_SELseed_logclusterWidth1", 100,-1,7);
      hSEL_dec_seed_clusterWidth2[iss][idec]  = new TH1F( samp+"_"+lab+"_SELseed_logclusterWidth2",  samp+"_"+lab+"_SELseed_logclusterWidth2", 100,-1,7);
      hSEL_dec_seed_clusterLength[iss][idec] = new TH1F( samp+"_"+lab+"_SELseed_logclusterLength", samp+"_"+lab+"_SELseed_logclusterLength", 100,-1,7);
      hSEL_dec_seed_energy[iss][idec]        = new TH1F( samp+"_"+lab+"_SELseed_energy",           samp+"_"+lab+"_SELseed_energy", 100,0,300);


      hSEL_dec_cone_nchpfo     [iss][idec] = new TH1F( samp+"_"+lab+"_SELcone_nchpfo",      samp+"_"+lab+"_SELcone_nchpfo", 8,-0.5,7.5 );
      hSEL_dec_cone_ntracks    [iss][idec] = new TH1F( samp+"_"+lab+"_SELcone_ntracks",     samp+"_"+lab+"_SELcone_ntracks", 8,-0.5,7.5 );
      hSEL_dec_cone_trackNTPC  [iss][idec] = new TH1F( samp+"_"+lab+"_SELcone_trackNTPC",   samp+"_"+lab+"_SELcone_trackNTPC", 50,0,250 );
      hSEL_dec_cone_nnhadpfo   [iss][idec] = new TH1F( samp+"_"+lab+"_SELcone_nnhadpfo",    samp+"_"+lab+"_SELcone_nnhadpfo", 8,-0.5,7.5 );
      hSEL_dec_cone_ngammapfo  [iss][idec] = new TH1F( samp+"_"+lab+"_SELcone_ngammapfo",   samp+"_"+lab+"_SELcone_ngammapfo", 8,-0.5,7.5 );
      hSEL_dec_cone_ncompoundpfo  [iss][idec] = new TH1F( samp+"_"+lab+"_SELcone_ncompoundpfo",   samp+"_"+lab+"_SELcone_ncompoundpfo", 8,-0.5,7.5 );
      hSEL_dec_cone_nhaden     [iss][idec] = new TH1F( samp+"_"+lab+"_SELcone_nhaden",      samp+"_"+lab+"_SELcone_nhaden",  100,0.,100 );
      hSEL_dec_cone_gammaen    [iss][idec] = new TH1F( samp+"_"+lab+"_SELcone_gammaen",     samp+"_"+lab+"_SELcone_gamma", 100,0.,300 );
      hSEL_dec_cone_nhadenFrac [iss][idec] = new TH1F( samp+"_"+lab+"_SELcone_nhadenFrac",  samp+"_"+lab+"_SELcone_nhadenFrac", 100,0.,1. );
      hSEL_dec_cone_gammaenFrac[iss][idec] = new TH1F( samp+"_"+lab+"_SELcone_gammaenFrac", samp+"_"+lab+"_SELcone_gammaenFrac", 100,0.,1. );

      hSEL_dec_cone_visMass    [iss][idec]    = new TH1F( samp+"_"+lab+"_SELcone_visMass",        samp+"_"+lab+"_SELcone_visMass", 100, 0, 3 );
      hSEL_dec_cone_neutralvisMass[iss][idec] = new TH1F( samp+"_"+lab+"_SELcone_neutralvisMass", samp+"_"+lab+"_SELcone_neutralvisMass", 100, 0, 3 );

      hSEL_dec_cone_vis_neutral_Mass[iss][idec] = new TH2F( samp+"_"+lab+"_SELcone_vis_neutral_Mass",   samp+"_"+lab+"_SELcone_vis_neutral_Mass", 100, 0, 3, 100, 0, 3 );
      hSEL_dec_cone_vis_neutral_Mass_0g[iss][idec] = new TH2F( samp+"_"+lab+"_SELcone_vis_neutral_Mass_0g",   samp+"_"+lab+"_SELcone_vis_neutral_Mass_0g", 100, 0, 3, 100, 0, 3 );
      hSEL_dec_cone_vis_neutral_Mass_1g[iss][idec] = new TH2F( samp+"_"+lab+"_SELcone_vis_neutral_Mass_1g",   samp+"_"+lab+"_SELcone_vis_neutral_Mass_1g", 100, 0, 3, 100, 0, 3 );
      hSEL_dec_cone_vis_neutral_Mass_2g[iss][idec] = new TH2F( samp+"_"+lab+"_SELcone_vis_neutral_Mass_2g",   samp+"_"+lab+"_SELcone_vis_neutral_Mass_2g", 100, 0, 3, 100, 0, 3 );
      hSEL_dec_cone_vis_neutral_Mass_3g[iss][idec] = new TH2F( samp+"_"+lab+"_SELcone_vis_neutral_Mass_3g",   samp+"_"+lab+"_SELcone_vis_neutral_Mass_3g", 100, 0, 3, 100, 0, 3 );
      hSEL_dec_cone_vis_neutral_Mass_4g[iss][idec] = new TH2F( samp+"_"+lab+"_SELcone_vis_neutral_Mass_4g",   samp+"_"+lab+"_SELcone_vis_neutral_Mass_4g", 100, 0, 3, 100, 0, 3 );

      hSEL_dec_cone_visMassDiff    [iss][idec]    = new TH1F( samp+"_"+lab+"_SELcone_visMassDiff",        samp+"_"+lab+"_SELcone_visMassDiff", 100, -1,1 );
      hSEL_dec_cone_neutralvisMassDiff[iss][idec] = new TH1F( samp+"_"+lab+"_SELcone_neutralvisMassDiff", samp+"_"+lab+"_SELcone_neutralvisMassDiff", 100, -1,1 );

      hSEL_dec_cone_visEnergyDiff    [iss][idec]    = new TH1F( samp+"_"+lab+"_SELcone_visEnergyDiff",        samp+"_"+lab+"_SELcone_visEnergyDiff", 100, -50,50 );
      hSEL_dec_cone_neutralvisEnergyDiff[iss][idec] = new TH1F( samp+"_"+lab+"_SELcone_neutralvisEnergyDiff", samp+"_"+lab+"_SELcone_neutralvisEnergyDiff", 100, -50,50 );

      hSEL_dec_coneTRIM_visMass    [iss][idec]    = new TH1F( samp+"_"+lab+"_SELconeTRIM_visMass",        samp+"_"+lab+"_SELconeTRIM_visMass", 100, 0, 2 );
      hSEL_dec_coneTRIM_neutralvisMass[iss][idec] = new TH1F( samp+"_"+lab+"_SELconeTRIM_neutralvisMass", samp+"_"+lab+"_SELconeTRIM_neutralvisMass", 100, 0, 2 );


      hSEL_dec_coneTRIM_vis_neutral_Mass[iss][idec] = new TH2F( samp+"_"+lab+"_SELconeTRIM_vis_neutral_Mass",   samp+"_"+lab+"_SELconeTRIM_vis_neutral_Mass", 100, 0, 2, 100, 0, 2 );
      hSEL_dec_coneTRIM_vis_neutral_Mass_0g[iss][idec] = new TH2F( samp+"_"+lab+"_SELconeTRIM_vis_neutral_Mass_0g",   samp+"_"+lab+"_SELconeTRIM_vis_neutral_Mass_0g", 100, 0, 2, 100, 0, 2 );
      hSEL_dec_coneTRIM_vis_neutral_Mass_1g[iss][idec] = new TH2F( samp+"_"+lab+"_SELconeTRIM_vis_neutral_Mass_1g",   samp+"_"+lab+"_SELconeTRIM_vis_neutral_Mass_1g", 100, 0, 2, 100, 0, 2 );
      hSEL_dec_coneTRIM_vis_neutral_Mass_2g[iss][idec] = new TH2F( samp+"_"+lab+"_SELconeTRIM_vis_neutral_Mass_2g",   samp+"_"+lab+"_SELconeTRIM_vis_neutral_Mass_2g", 100, 0, 2, 100, 0, 2 );
      hSEL_dec_coneTRIM_vis_neutral_Mass_3g[iss][idec] = new TH2F( samp+"_"+lab+"_SELconeTRIM_vis_neutral_Mass_3g",   samp+"_"+lab+"_SELconeTRIM_vis_neutral_Mass_3g", 100, 0, 2, 100, 0, 2 );
      hSEL_dec_coneTRIM_vis_neutral_Mass_4g[iss][idec] = new TH2F( samp+"_"+lab+"_SELconeTRIM_vis_neutral_Mass_4g",   samp+"_"+lab+"_SELconeTRIM_vis_neutral_Mass_4g", 100, 0, 2, 100, 0, 2 );

      hSEL_dec_coneTRIM_visMassDiff    [iss][idec]    = new TH1F( samp+"_"+lab+"_SELconeTRIM_visMassDiff",        samp+"_"+lab+"_SELconeTRIM_visMassDiff", 100, -1,1 );
      hSEL_dec_coneTRIM_neutralvisMassDiff[iss][idec] = new TH1F( samp+"_"+lab+"_SELconeTRIM_neutralvisMassDiff", samp+"_"+lab+"_SELconeTRIM_neutralvisMassDiff", 100, -1,1 );

      hSEL_dec_coneTRIM_visEnergyDiff    [iss][idec]    = new TH1F( samp+"_"+lab+"_SELconeTRIM_visEnergyDiff",        samp+"_"+lab+"_SELconeTRIM_visEnergyDiff", 100, -50,50 );
      hSEL_dec_coneTRIM_neutralvisEnergyDiff[iss][idec] = new TH1F( samp+"_"+lab+"_SELconeTRIM_neutralvisEnergyDiff", samp+"_"+lab+"_SELconeTRIM_neutralvisEnergyDiff", 100, -50,50 );

    }


    // selection observables

    h_seed_dZ0[iss] = new TH1F( samp+"seed_dZ0", samp+"seed_dZ0", 100, -1,1);
    h_seed_dZ0ave[iss] = new TH1F( samp+"seed_avedZ0", samp+"seed_avedZ0", 100, -1,1);

    h_maxSeedEn[iss] = new TH1F( samp+"maxSeedEn", samp+"maxSeedEn", 100, 0, 300 );
    h_minSeedEn[iss] = new TH1F( samp+"minSeedEn", samp+"minSeedEn", 100, 0, 300 );
    h_maxSeedEn_eonp[iss] = new TH1F( samp+"maxSeedEn_eonp", samp+"maxSeedEn_eonp", 100, 0, 2. );
    h_minSeedEn_eonp[iss] = new TH1F( samp+"minSeedEn_eonp", samp+"minSeedEn_eonp", 100, 0, 2. );
    h_maxSeedEn_caloen[iss] = new TH1F( samp+"maxSeedEn_caloen", samp+"maxSeedEn_caloen", 100, 0, 50. );
    h_minSeedEn_caloen[iss] = new TH1F( samp+"minSeedEn_caloen", samp+"minSeedEn_caloen", 100, 0, 50. );

    h_seedClusterWidth1[iss] = new TH1F( samp+"logseedClusterWidth1",  samp+"logseedClusterWidth1", 100, -1, 7 );
    h_seedClusterWidth2[iss] = new TH1F( samp+"logseedClusterWidth2",  samp+"logseedClusterWidth2", 100, -1, 7 );
    h_seedClusterLength[iss] = new TH1F( samp+"logseedClusterLength",  samp+"logseedClusterLength", 100, -1, 7 );

    h_seedEcalEn[iss] = new TH1F( samp+"seedEcalEn", samp+"seedEcalEn", 100, 0, 300 );
    h_seedEcalEonP[iss] = new TH1F( samp+"seedEcalEonP", samp+"seedEcalEonP", 100, 0, 3 );    
    h_ooconeMaxGammaEn[iss] = new TH1F( samp+"ooconeMaxGammaEn", samp+"ooconeMaxGammaEn", 100, 0, 150 );
    h_ooconeMaxPFOEn[iss] = new TH1F( samp+"ooconeMaxPFOEn", samp+"ooconeMaxPFOEn", 100, 0, 150 );
    h_coneMass     [iss] = new TH1F( samp+"coneMass", samp+"coneMass", 100, 0, 2 );
    h_seedSeedAngle[iss] = new TH1F( samp+"seedSeedAngle", samp+"seedSeedAngle", 100, 0, 0.25 );
    h_seedSeedDphi [iss] = new TH1F( samp+"seedSeedDphi", samp+"seedSeedDphi", 100, 0, 0.25 );
    h_maxSeedCosth[iss] = new TH1F( samp+"maxSeedCosth", samp+"maxSeedCosth", 100, 0, 1. );
    h_minSeedCosth[iss] = new TH1F( samp+"minSeedCosth", samp+"minSeedCosth", 100, 0, 1. );

    hSEL_maxSeedEn[iss]        = new TH1F( samp+"SELmaxSeedEn",        samp+"SELmaxSeedEn", 100, 0, 300 );
    hSEL_minSeedEn[iss]        = new TH1F( samp+"SELminSeedEn",        samp+"SELminSeedEn", 100, 0, 300 );
    hSEL_maxSeedEn_eonp[iss]        = new TH1F( samp+"SELmaxSeedEn_eonp",        samp+"SELmaxSeedEn_eonp", 100, 0, 2. );
    hSEL_minSeedEn_eonp[iss]        = new TH1F( samp+"SELminSeedEn_eonp",        samp+"SELminSeedEn_eonp", 100, 0, 2. );
    hSEL_maxSeedEn_caloen[iss]        = new TH1F( samp+"SELmaxSeedEn_caloen",        samp+"SELmaxSeedEn_caloen", 100, 0, 50. );
    hSEL_minSeedEn_caloen[iss]        = new TH1F( samp+"SELminSeedEn_caloen",        samp+"SELminSeedEn_caloen", 100, 0, 50. );
    hSEL_seedEcalEn[iss]       = new TH1F( samp+"SELseedEcalEn",       samp+"SELseedEcalEn", 100, 0, 300 );
    hSEL_seedEcalEonP[iss]     = new TH1F( samp+"SELseedEcalEonP",     samp+"SELseedEcalEonP", 100, 0, 3 );    
    hSEL_ooconeMaxGammaEn[iss] = new TH1F( samp+"SELooconeMaxGammaEn", samp+"SELooconeMaxGammaEn", 100, 0, 150 );
    hSEL_ooconeMaxPFOEn[iss] = new TH1F( samp+"SELooconeMaxPFOEn", samp+"SELooconeMaxPFOEn", 100, 0, 150 );
    hSEL_coneMass     [iss] = new TH1F( samp+"SEL_coneMass", samp+"SEL_coneMass", 100, 0, 2 );
    hSEL_seedSeedAngle[iss] = new TH1F( samp+"SEL_seedSeedAngle", samp+"SEL_seedSeedAngle", 100, 0., 0.25 );
    hSEL_seedSeedDphi [iss] = new TH1F( samp+"SEL_seedSeedDphi", samp+"SEL_seedSeedDphi", 100, 0., 0.25 );
    hSEL_maxSeedCosth[iss]        = new TH1F( samp+"SELmaxSeedCosth",        samp+"SELmaxSeedCosth", 100, 0, 1. );
    hSEL_minSeedCosth[iss]        = new TH1F( samp+"SELminSeedCosth",        samp+"SELminSeedCosth", 100, 0, 1. );

    hSEL_tauMinusCosth[iss] = new TH1F( samp+"SELtauMinusCosth",  samp+"SELtauMinusCosth", 100, -1., 1. );

    hSEL_mcdec[iss] = new TH1F( samp+"SEL_mcdec", samp+"SEL_mcdec", tauUtils::NDECAYS+1,-1.5,tauUtils::NDECAYS-0.5 );

    hSEL_rec_rho_mcdec[iss] = new TH1F( samp+"SEL_rec_rho_mcdec", samp+"SEL_rec_rho_mcdec", tauUtils::NDECAYS+1,-1.5,tauUtils::NDECAYS-0.5 );
    hSEL_rec_rho_pol  [iss] = new TH1F( samp+"SEL_rec_rho_pol"  , samp+"SEL_rec_rho_pol"  , 150, -1.5, 1.5 );
    hSEL_rec_rho_pol_MCpos  [iss] = new TH1F( samp+"SEL_rec_rho_pol_MCpos"  , samp+"SEL_rec_rho_pol_MCpos"  , 150, -1.5, 1.5 );
    hSEL_rec_rho_pol_MCneg  [iss] = new TH1F( samp+"SEL_rec_rho_pol_MCneg"  , samp+"SEL_rec_rho_pol_MCneg"  , 150, -1.5, 1.5 );
    hSEL_rec_rho_pol_MCoth  [iss] = new TH1F( samp+"SEL_rec_rho_pol_MCoth"  , samp+"SEL_rec_rho_pol_MCoth"  , 150, -1.5, 1.5 );
    hSEL_rec_rho_coneMass [iss] = new TH1F( samp+"SEL_rec_rho_coneMass"  ,samp+"SEL_rec_rho_coneMass"  , 100, 0, 2 );

    hSEL_rec_a1p_mcdec[iss] = new TH1F( samp+"SEL_rec_a1p_mcdec", samp+"SEL_rec_a1p_mcdec", tauUtils::NDECAYS+1,-1.5,tauUtils::NDECAYS-0.5 );
    hSEL_rec_a1p_pol  [iss] = new TH1F( samp+"SEL_rec_a1p_pol"  , samp+"SEL_rec_a1p_pol"  , 100, -1., 1. );
    hSEL_rec_a1p_coneMass [iss] = new TH1F( samp+"SEL_rec_a1p_coneMass"  ,samp+"SEL_rec_a1p_coneMass"  , 100, 0, 2 );

    hSEL_rec_pi_mcdec [iss] = new TH1F( samp+"SEL_rec_pi_mcdec" , samp+"SEL_rec_pi_mcdec" , tauUtils::NDECAYS+1,-1.5,tauUtils::NDECAYS-0.5 );
    hSEL_rec_pi_pol   [iss] = new TH1F( samp+"SEL_rec_pi_pol"   , samp+"SEL_rec_pi_pol"   , 150, -1.5, 1.5 );
    hSEL_rec_pi_pol_MCpos   [iss] = new TH1F( samp+"SEL_rec_pi_pol_MCpos"   , samp+"SEL_rec_pi_pol_MCpos"   , 150, -1.5, 1.5 );
    hSEL_rec_pi_pol_MCneg   [iss] = new TH1F( samp+"SEL_rec_pi_pol_MCneg"   , samp+"SEL_rec_pi_pol_MCneg"   , 150, -1.5, 1.5 );
    hSEL_rec_pi_pol_MCoth   [iss] = new TH1F( samp+"SEL_rec_pi_pol_MCoth"   , samp+"SEL_rec_pi_pol_MCoth"   , 150, -1.5, 1.5 );
    hSEL_rec_pi_coneMass [iss] = new TH1F( samp+"SEL_rec_pi_coneMass"  ,samp+"SEL_rec_pi_coneMass"  , 100, 0, 2 );

    hSEL_recpi_mcall_pol [iss] = new TH2F(  samp+"SEL_recpi_mcall_pol",  samp+"SEL_recpi_mcall_pol",  150, -1.5, 1.5, 150, -1.5, 1.5 );
    hSEL_recpi_mcpi_pol  [iss] = new TH2F(  samp+"SEL_recpi_mcpi_pol",   samp+"SEL_recpi_mcpi_pol",   150, -1.5, 1.5, 150, -1.5, 1.5 );
    hSEL_recrho_mcall_pol[iss] = new TH2F(  samp+"SEL_recrho_mcall_pol", samp+"SEL_recrho_mcall_pol", 150, -1.5, 1.5, 150, -1.5, 1.5 );
    hSEL_recrho_mcrho_pol[iss] = new TH2F(  samp+"SEL_recrho_mcrho_pol", samp+"SEL_recrho_mcrho_pol", 150, -1.5, 1.5, 150, -1.5, 1.5 );
    
    hSEL_recpi_mcall_dpol [iss] = new TH1F(  samp+"SEL_recpi_mcall_dpol",  samp+"SEL_recpi_mcall_dpol",  500, -1.5, 1.5 );
    hSEL_recpi_mcpi_dpol  [iss] = new TH1F(  samp+"SEL_recpi_mcpi_dpol",   samp+"SEL_recpi_mcpi_dpol",   500, -1.5, 1.5 );
    hSEL_recrho_mcall_dpol[iss] = new TH1F(  samp+"SEL_recrho_mcall_dpol", samp+"SEL_recrho_mcall_dpol", 500, -1.5, 1.5 );
    hSEL_recrho_mcrho_dpol[iss] = new TH1F(  samp+"SEL_recrho_mcrho_dpol", samp+"SEL_recrho_mcrho_dpol", 500, -1.5, 1.5 );
    
    //    TString sdec[3]={"pipi","pirho","rhorho"};
    // for (int idecdec=0; idecdec<3; idecdec++) {
    //   //_hFit_nSolutions_cone_ip[iss][idecdec] = new TH2F( samp+sdec[idecdec]+"_nSolutions_cone_ip",  samp+sdec[idecdec]+"_nSolutions_cone_ip", 5,-0.5,4.5, 5,-0.5,4.5);
    //   //
    //   //_hFit_nSolutions_cone[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_nSolutions_cone",  samp+sdec[idecdec]+"_nSolutions_cone", 5,-0.5,4.5);
    //   //_hFit_nSolutions_ip[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_nSolutions_ip",  samp+sdec[idecdec]+"_nSolutions_ip", 5,-0.5,4.5);
    //   //
    //   //_hFit_coneFit_pol12[iss][idecdec] = new TH2F( samp+sdec[idecdec]+"_coneFit_pol12",  samp+sdec[idecdec]+"_coneFit_pol12", 20,-1,1,20,-1,1);
    //   //_hFit_coneFit_pol[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_coneFit_pol",  samp+sdec[idecdec]+"_coneFit_pol", 20,-1,1);
    //   //_hFit_coneFit_fitFlag[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_coneFit_fitFlag",  samp+sdec[idecdec]+"_coneFit_fitFlag", 10,-0.5, 9.5);
    // 
    //   // _hFit_ipFit_pol12[iss][idecdec] = new TH2F( samp+sdec[idecdec]+"_ipFit_pol12",  samp+sdec[idecdec]+"_ipFit_pol12", 20,-1,1,20,-1,1);
    //   // _hFit_ipFit_pol[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_ipFit_pol",  samp+sdec[idecdec]+"_ipFit_pol", 20,-1,1);
    //   // 
    //   // _hFit_ipFit_allSols_totE[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_ipFit_allSols_totE", samp+sdec[idecdec]+"_ipFit_allSols_totE",100,0,1000);
    //   // _hFit_ipFit_allSols_totPt[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_ipFit_allSols_totPt", samp+sdec[idecdec]+"_ipFit_allSols_totPt",100,0,500);
    //   // _hFit_ipFit_allSols_totPz[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_ipFit_allSols_totPz", samp+sdec[idecdec]+"_ipFit_allSols_totPz",100,-200,200);
    //   // _hFit_ipFit_allSols_decayLength[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_ipFit_allSols_decL", samp+sdec[idecdec]+"_ipFit_allSols_decL",100,-250,250);
    //   // 
    //   // _hFit_ipFit_reasonableSols_totE[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_ipFit_reasonableSols_totE", samp+sdec[idecdec]+"_ipFit_reasonableSols_totE",100,0,1000);
    //   // _hFit_ipFit_reasonableSols_totPt[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_ipFit_reasonableSols_totPt", samp+sdec[idecdec]+"_ipFit_reasonableSols_totPt",100,0,500);
    //   // _hFit_ipFit_reasonableSols_totPz[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_ipFit_reasonableSols_totPz", samp+sdec[idecdec]+"_ipFit_reasonableSols_totPz",100,-200,200);
    //   // _hFit_ipFit_reasonableSols_decayLength[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_ipFit_reasonableSols_decL", samp+sdec[idecdec]+"_ipFit_reasonableSols_decL",100,-250,250);
    //   // 
    //   // _hFit_ipFit_bestSol_totE[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_ipFit_bestSol_totE", samp+sdec[idecdec]+"_ipFit_bestSol_totE",100,0,1000);
    //   // _hFit_ipFit_bestSol_totPt[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_ipFit_bestSol_totPt", samp+sdec[idecdec]+"_ipFit_bestSol_totPt",100,0,500);
    //   // _hFit_ipFit_bestSol_totPz[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_ipFit_bestSol_totPz", samp+sdec[idecdec]+"_ipFit_bestSol_totPz",100,-200,200);
    //   // _hFit_ipFit_bestSol_decayLength[iss][idecdec] = new TH1F( samp+sdec[idecdec]+"_ipFit_bestSol_decL", samp+sdec[idecdec]+"_ipFit_bestSol_decL",100,-250,250);
    // }

    // TString tdec[2]={"pi","rho"};
    // for (int idecdec=0; idecdec<2; idecdec++) {
    //   // _hFit_coneFit_polRecMc[iss][idecdec] = new TH2F( samp+tdec[idecdec]+"_coneFit_polRecMc", samp+tdec[idecdec]+"_coneFit_polRecMc",  20,-1,1,20,-1,1);
    //   // _hFit_coneFit_polRec[iss][idecdec] = new TH1F( samp+tdec[idecdec]+"_coneFit_polRec", samp+tdec[idecdec]+"_coneFit_polRec",  20,-1,1);
    //   // _hFit_coneFit_polRecMcDiff[iss][idecdec] = new TH1F( samp+tdec[idecdec]+"_coneFit_polRecMcDiff", samp+tdec[idecdec]+"_coneFit_polRecMcDiff",  100,-2,2);
    //   // 
    //   // _hFit_ipFit_polRecMc[iss][idecdec] = new TH2F( samp+tdec[idecdec]+"_ipFit_polRecMc", samp+tdec[idecdec]+"_ipFit_polRecMc",  20,-1,1,20,-1,1);
    //   // _hFit_ipFit_polRec[iss][idecdec] = new TH1F( samp+tdec[idecdec]+"_ipFit_polRec", samp+tdec[idecdec]+"_ipFit_polRec",  20,-1,1);
    //   // _hFit_ipFit_polRecMcDiff[iss][idecdec] = new TH1F( samp+tdec[idecdec]+"_ipFit_polRecMcDiff", samp+tdec[idecdec]+"_ipFit_polRecMcDiff",  100,-2,2);
    // }

  }




  for (int i=0; i<NCLASS; i++) {
    _nOrig[i]=0;
    _nPresel[i]=0;
    _nTwoSeeds[i]=0;
    _nTwoWellMatchedSeeds[i]=0;
    _nGoodSeedDir[i]=0;
    _nSel[i]=0;
    _nSel_secondseenen[i]=0;
    _nSel_tmass[i]=0;
    _nSel_outofcone[i]=0;
    _nSel_acoLin[i]=0;
    _nSel_acoPlan[i]=0;
    _nSel_chg[i]=0;
    _nSel_lepton[i]=0;
    //    _nSel_seedcaloen[i]=0;
    //    _nSel_singleElec[i]=0;
    _nSel_seedCluster[i]=0;
    _nSel_isr[i]=0;

    _nCumulSel[i]=0;
    _nCumulSel_secondseenen[i]=0;
    _nCumulSel_tmass[i]=0;
    _nCumulSel_outofcone[i]=0;
    _nCumulSel_acoLin[i]=0;
    _nCumulSel_acoPlan[i]=0;
    _nCumulSel_chg[i]=0;
    _nCumulSel_lepton[i]=0;
    //    _nCumulSel_seedcaloen[i]=0;
    //    _nCumulSel_singleElec[i]=0;
    _nCumulSel_seedCluster[i]=0;
    _nCumulSel_isr[i]=0;
  }

  return;
}

void danielKeitaTauFinderProcessor::processRunHeader( LCRunHeader* run) { 
  cout << "hello from danielKeitaTauFinderProcessor::processRunHeader : " ;

  for (int i=0; i<NCLASS; i++) {
    cout << _nOrig[i] << " ";
  }
  cout << endl;

}

void danielKeitaTauFinderProcessor::processEvent( LCEvent * evt ) { 

  //  cout << "processEvent " << evt->getEventNumber() << endl;

  int isample(0); // MC sample - 0:high mass tt, 1:med mass tt, 2: lowmass tt, 3: mm


  float tautauInvMass(0);
  float taucosth[2]={-1,-1};
  std::vector <MCParticle*> stableMCtaudaughters[2];
  std::vector <MCParticle*> stableMCtaudaughtersInCone[2];
  std::vector <MCParticle*> stableMCallInCone[2];
  std::vector <MCParticle*> stableMCallNeutralInCone[2];
  std::vector <MCParticle*> stableMCneutraltaudaughtersInCone[2];
  std::vector <MCParticle*> allMCtaudaughters[2];
  std::vector <MCParticle*> finalmctaus;
  int MCdecayMode[2]={-1,-1};
  
  TLorentzVector mcTau4mom[2];
  TLorentzVector mcTauVis4mom[2];
  TLorentzVector mcTauNeutralVis4mom[2];

  float mctauhelicity[2]={0,0};
  float mcexactPolarimeter[2]={-99,-99};
  float mcapproxPolarimeter[2]={-99,-99};
  float mcapproxCheatEnPolarimeter[2]={-99,-99};


  MCParticle* mctauminus(0);

  TVector3 pvposMC(-99,-99,-99);

  float mctauMinusCosTh(-2);

  classLabels[0]="2T,HM,h";
  classLabels[1]="2T,HM,l";
  classLabels[2]="2T,MM";
  classLabels[3]="2T,LM";
  classLabels[4]="0T";
  classLabels[5]="1T";
  classLabels[6]="gt2T";


  try {
    LCCollection* mccol = evt->getCollection( "MCParticle" );

    finalmctaus = tauUtils::findFinalTaus( mccol );

    if ( finalmctaus.size()==0 ) isample=4;
    else if ( finalmctaus.size()==1 ) isample=5;
    else if ( finalmctaus.size()!=2 ) isample=6;
    else {

      for (size_t i=0; i<finalmctaus.size(); i++) {

	if ( pvposMC.X()<-98 ) {
	  pvposMC.SetXYZ( finalmctaus[i]->getVertex()[0] , finalmctaus[i]->getVertex()[1] , finalmctaus[i]->getVertex()[2] );
	}

	stableMCtaudaughters[i] = tauUtils::getstablemctauDaughters( finalmctaus[i] );
	allMCtaudaughters[i] = tauUtils::getmctauDaughters( finalmctaus[i] );
	MCdecayMode[i] = tauUtils::getMCdecayMode( stableMCtaudaughters[i] );

	mcTau4mom[i] = tauUtils::getTLV( finalmctaus[i] );
	mcTauVis4mom[i] = tauUtils::getVisibleTau4mom( finalmctaus[i] );
	mcTauNeutralVis4mom[i] = tauUtils::getVisibleNeutralTau4mom( finalmctaus[i] );

	stableMCtaudaughtersInCone[i].clear();
	stableMCallInCone[i].clear();
	stableMCallNeutralInCone[i].clear();
	// find seed MC particle
	MCParticle* mcseed(0);
	for ( auto m : stableMCtaudaughters[i] ) {
	  if ( isCharged ( m ) && ( !mcseed || m->getEnergy()>mcseed->getEnergy() ) ) {
	    mcseed=m;
	  }
	}
	stableMCtaudaughtersInCone[i].push_back( mcseed );
	stableMCallInCone[i].push_back( mcseed );

	TVector3 seedDir( mcseed->getMomentum() );
	for ( auto m : stableMCtaudaughters[i] ) {
	  if (  m != mcseed && 
		! tauUtils::isNeutrino( m ) && 
		TVector3( m->getMomentum() ).Angle( seedDir ) < _conesize ) {
	    stableMCtaudaughtersInCone[i].push_back( m );
	    if (isNeutral(m)) {
	      stableMCneutraltaudaughtersInCone[i].push_back( m );
	    }
	  }
	}

	for (int j=0; j<mccol->getNumberOfElements(); j++) {
	  MCParticle* mcp =  dynamic_cast<MCParticle*> (mccol->getElementAt(j));

	  //	  if ( j<20 ) 
	  //	    cout << j << " " << mcp->getPDG() << " " << mcp->getGeneratorStatus() << " "<<  TVector3(mcp->getMomentum() ).Angle( seedDir );

	  if (  mcp != mcseed && 
		mcp->getGeneratorStatus()==1 &&
		! tauUtils::isNeutrino( mcp ) && 
		TVector3( mcp->getMomentum() ).Angle( seedDir ) < _conesize ) {
	    stableMCallInCone[i].push_back( mcp );
	    if (isNeutral(mcp)) {
	      stableMCallNeutralInCone[i].push_back( mcp );
	    }
	    //	    if ( j<20 ) 
	    //	      cout << " ADDING..." ;

	  }

	  //	  if ( j<20 ) 
	  //	    cout << endl;

	}



      }

      tautauInvMass = ( mcTau4mom[0] + mcTau4mom[1] ).M();
      if ( tautauInvMass > 480 ) {
	if ( MCdecayMode[0] != tauUtils::decayEl && MCdecayMode[0] != tauUtils::decayMu &&
	     MCdecayMode[1] != tauUtils::decayEl && MCdecayMode[1] != tauUtils::decayMu ) {
	  isample=0;
	} else {
	  isample=1;
	}
      } else if ( tautauInvMass > 75 ) isample=2;
      else isample=3;

      taucosth[0] = fabs( mcTau4mom[0].Vect().CosTheta() );
      taucosth[1] = fabs( mcTau4mom[1].Vect().CosTheta() );
      if ( taucosth[0] > taucosth[1] ) {
	float temp = taucosth[0];
	taucosth[0] = taucosth[1];
	taucosth[1] = temp;
      }

      
      if      ( finalmctaus[0]->getCharge()<0 ) mctauminus=finalmctaus[0];
      else if ( finalmctaus[1]->getCharge()<0 ) mctauminus=finalmctaus[1];
      if ( mctauminus ) {
	mctauMinusCosTh = tauUtils::getTLV(mctauminus).Vect().CosTheta();
	h_mc_tauMinus_costh[isample]->Fill( mctauMinusCosTh );
	h_mctautau_ecom_tauMinusCosth->Fill(tautauInvMass, mctauMinusCosTh);
      }

      h_mc_tauSpin[isample]->Fill( finalmctaus[0]->getSpin()[2],  finalmctaus[1]->getSpin()[2] );


    }

    _nOrig[isample]++;

    //    static TVector3 getPolarimeter_rho ( TLorentzVector charged, TLorentzVector neutral, TLorentzVector neutrino );
    //    static TVector3 getPolarimeter_pi  ( TLorentzVector charged, TLorentzVector neutrino );

    if ( finalmctaus.size()==2 ) {


      for (size_t i=0; i<finalmctaus.size(); i++) {

	mctauhelicity[i] = finalmctaus[i]->getCharge()*finalmctaus[i]->getSpin()[2] > 0 ? -1 : +1 ;

	if (  MCdecayMode[i] == tauUtils::decayChPi ) {

	  TLorentzVector tlv_charged;
	  TLorentzVector tlv_neutrino;
	  float mc_pi_en(-1);

	  for (int j=0; j< stableMCtaudaughters[i].size(); j++ ) {
	    if ( abs( stableMCtaudaughters[i][j]->getPDG()) == 211 ) {
	      tlv_charged.SetXYZT( stableMCtaudaughters[i][j]->getMomentum()[0],
				   stableMCtaudaughters[i][j]->getMomentum()[1],
				   stableMCtaudaughters[i][j]->getMomentum()[2],
				   stableMCtaudaughters[i][j]->getEnergy() );
	      mc_pi_en = stableMCtaudaughters[i][j]->getEnergy();
	    } else if ( abs( stableMCtaudaughters[i][j]->getPDG()) == 16 ) {
	      tlv_neutrino.SetXYZT( stableMCtaudaughters[i][j]->getMomentum()[0],
				    stableMCtaudaughters[i][j]->getMomentum()[1],
				    stableMCtaudaughters[i][j]->getMomentum()[2],
				    stableMCtaudaughters[i][j]->getEnergy() );
	    }
	  }

	  // switch sign for exact
	  mcexactPolarimeter[i] = -1.*cos( tauUtils::getPolarimeter_pi( tlv_charged, tlv_neutrino ).Angle( TVector3( finalmctaus[i]->getMomentum() ) ) );
	  mcapproxPolarimeter[i] = tauUtils::getLongPolarimeterFromEn_pi( tlv_charged, 250. );
	  mcapproxCheatEnPolarimeter[i] = tauUtils::getLongPolarimeterFromEn_pi( tlv_charged, mcTau4mom[i].E() );

	  h_pi_mcPolarExact[isample]->Fill(mcexactPolarimeter[i]);
	  h_pi_mcPolarApprox[isample]->Fill(mcapproxPolarimeter[i]);
	  h_pi_mcPolarApproxCheatEn[isample]->Fill(mcapproxCheatEnPolarimeter[i]);

	  if ( mctauhelicity[i]<0 ) {
	    h_pi_mcPolarExact_helNeg[isample]->Fill(mcexactPolarimeter[i]);
	    h_pi_mcPolarApprox_helNeg[isample]->Fill(mcapproxPolarimeter[i]);
	    h_pi_mcPolarApproxCheatEn_helNeg[isample]->Fill(mcapproxCheatEnPolarimeter[i]);
	  } else {
	    h_pi_mcPolarExact_helPos[isample]->Fill(mcexactPolarimeter[i]);
	    h_pi_mcPolarApprox_helPos[isample]->Fill(mcapproxPolarimeter[i]);
	    h_pi_mcPolarApproxCheatEn_helPos[isample]->Fill(mcapproxCheatEnPolarimeter[i]);
	  }

	  h_pirho_tauMinusCosth_mcPolarExact[isample]->Fill(mctauMinusCosTh, mcexactPolarimeter[i] );

	  //	  float polar = 2*( mc_pi_en/250. -0.5 ); // normalise to -1 -> 1

	  h_pi_mcPolar[isample]->Fill( mc_pi_en/250. );

	  h_pi_mcPolarComp[isample]->Fill( mcexactPolarimeter[i], mc_pi_en/250. );

	} else if (  MCdecayMode[i] == tauUtils::decayRho ) {
	  float polar = -1;
	  float pi0en(0);
	  float piChgen(0);

	  TLorentzVector tlv_charged;
	  TLorentzVector tlv_neutrino;
	  TLorentzVector tlv_neutral;


	  for (int j=0; j<  allMCtaudaughters[i] .size(); j++ ) {
	    if ( abs( allMCtaudaughters[i][j]->getPDG()) == 211 ) {
	      piChgen= allMCtaudaughters[i][j]->getEnergy();
	      tlv_charged.SetXYZT( allMCtaudaughters[i][j]->getMomentum()[0],
				   allMCtaudaughters[i][j]->getMomentum()[1],
				   allMCtaudaughters[i][j]->getMomentum()[2],
				   allMCtaudaughters[i][j]->getEnergy() );
	    } else if ( abs( allMCtaudaughters[i][j]->getPDG()) == 16 ) {
	      tlv_neutrino.SetXYZT( allMCtaudaughters[i][j]->getMomentum()[0],
				    allMCtaudaughters[i][j]->getMomentum()[1],
				    allMCtaudaughters[i][j]->getMomentum()[2],
				    allMCtaudaughters[i][j]->getEnergy() );
	    } else if ( allMCtaudaughters[i][j]->getPDG() == 111 ) {
	      tlv_neutral.SetXYZT( allMCtaudaughters[i][j]->getMomentum()[0],
				   allMCtaudaughters[i][j]->getMomentum()[1],
				   allMCtaudaughters[i][j]->getMomentum()[2],
				   allMCtaudaughters[i][j]->getEnergy() );
	      pi0en= allMCtaudaughters[i][j]->getEnergy();
	    }
	  }
	  polar = (piChgen-pi0en)/250.;

	  mcexactPolarimeter[i] = -1.*cos( tauUtils::getPolarimeter_rho(tlv_charged, tlv_neutral, tlv_neutrino ).Angle( TVector3( finalmctaus[i]->getMomentum() ) ) );

//	  float ebeam=250.;
//
//	  float etau=mcTau4mom[i].E();
//
//	  float mtau=1.777;
//	  float mtausq=pow(mtau,2);
//
//	  //	  float sqrts = 2.*etau;
//	  float sqrts = 2.*ebeam;
//	  float s = pow(sqrts,2);
//
//	  float eh = ( tlv_charged+tlv_neutral ).E();
//	  
//	  float x = 2.*eh/sqrts;
//	  float Qsq = (tlv_charged+tlv_neutral).M2();
//
//
//	  float cosPsi = ( x*( mtausq + Qsq ) - 2.*Qsq ) / ( (mtausq-Qsq)*sqrt( x*x - 4*Qsq/s )  );
//	  
//	  //float cosTheta = (2.*eh/etau) - 1.;
//
//	  //// a. do it explicitly using the tlvs
//	  //TLorentzVector h = tlv_charged+tlv_neutral;
//	  //h.Boost( - mcTau4mom[i].BoostVector() );
//	  //float costh_ceck = cos ( h.Angle( mcTau4mom[i].Vect() ) );
//
//
//	  // b. get from the energies in lab using true tau energy
//	  float gam = etau/mtau;
//	  float beta = sqrt( gam*gam - 1 )/gam;
//	  float ph_rf = (mtausq - Qsq)/(2*mtau);
//	  float eh_rf = sqrt( ph_rf*ph_rf + Qsq );
//	  float costh_exact = ( eh - gam*eh_rf ) / ( gam*beta*ph_rf );
//
//	  // I checked that a and b are consistent if we use true tau energy
//
//	  // c. get from the energies in lab using ebeam
//	  gam = ebeam/mtau;
//	  beta = sqrt( gam*gam - 1 )/gam;
//	  ph_rf = (mtausq - Qsq)/(2*mtau);
//	  eh_rf = sqrt( ph_rf*ph_rf + Qsq );
//	  float costh_ebeam = ( eh - gam*eh_rf ) / ( gam*beta*ph_rf );
//
//
//
//	  //	  cout << "BLAH  from tlv " << costh_ceck << " : using tauEn " << costh_exact << " using beamEn " << costh_ebeam << endl;
//
//	  TLorentzVector h=tlv_charged+tlv_neutral;
//	  TLorentzVector n=tlv_charged;
//	  n.Boost( -h.BoostVector() );
//	  float cosBeta = cos ( n.Vect().Angle( h.BoostVector() ) );
//	  
//	  float sin2Psi = sin( 2.*acos(cosPsi ) );
//	  float sinTheta = sqrt( 1. - pow( costh_exact,2 ) );
//
//	  float threeCosPsiterm = (3.*cosPsi-1)/2.;
//	  float threeCos2Betaterm = (3.*pow(cosBeta,2)-1)/2.;
//	  float mtausqOnQsq = mtausq/Qsq;

	  mcapproxCheatEnPolarimeter[i] = tauUtils::getLongPolarimeterFromEn_rho( tlv_charged, tlv_neutral, mcTau4mom[i].E() );
// (
// 					   ( -2.+mtausqOnQsq+2.*(1.+mtausqOnQsq)*threeCosPsiterm*threeCos2Betaterm )*costh_exact + 
// 					   3*sqrt(mtausqOnQsq)*threeCos2Betaterm*sin2Psi*sinTheta  
// 					   ) / ( 
// 						2. + mtausqOnQsq - 2.*(1.-mtausqOnQsq)*threeCosPsiterm*threeCos2Betaterm 
// 						 ) ;
	  
	  mcapproxPolarimeter[i] = tauUtils::getLongPolarimeterFromEn_rho( tlv_charged, tlv_neutral, 250. );
//(
//				    ( -2.+mtausqOnQsq+2.*(1.+mtausqOnQsq)*threeCosPsiterm*threeCos2Betaterm )*costh_ebeam + 
//				    3*sqrt(mtausqOnQsq)*threeCos2Betaterm*sin2Psi*sinTheta  
//				    ) / ( 
//					 2. + mtausqOnQsq - 2.*(1.-mtausqOnQsq)*threeCosPsiterm*threeCos2Betaterm 
//					  ) ;
	  


	  //	  h_rho_mcPolar[isample]->Fill(piChgen/250., pi0en/250.);
	  h_rho_mcPolar[isample]->Fill( (piChgen+pi0en)/250., fabs(piChgen-pi0en)/250.);
	  h_rho_mcPolar2[isample]->Fill( pi0en/(piChgen+pi0en), (piChgen+pi0en)/250. );

	  h_rho_mcPolarExact[isample]->Fill(mcexactPolarimeter[i]);

	  h_rho_mcPolarApprox[isample]->Fill(mcapproxPolarimeter[i]);
	  h_rho_mcPolarApproxCheatEn[isample]->Fill(mcapproxCheatEnPolarimeter[i]);

	  h_rho_mcPolarExactApprox[isample]->Fill(mcexactPolarimeter[i], mcapproxPolarimeter[i]);

          if ( mctauhelicity[i]<0 ) {

	    //h_rho_mcapprox_cospsi_helNeg[isample] -> Fill(cosPsi);
	    //h_rho_mcapprox_costheta_helNeg[isample] -> Fill(costh_exact);
	    //h_rho_mcapprox_cosbeta_helNeg[isample] -> Fill(cosBeta);
	  


	    h_rho_mcPolarExact_helNeg[isample]->Fill(mcexactPolarimeter[i]);
	    h_rho_mcPolarApprox_helNeg[isample]->Fill(mcapproxPolarimeter[i]);
	    h_rho_mcPolarApproxCheatEn_helNeg[isample]->Fill(mcapproxCheatEnPolarimeter[i]);
	  } else {

	    //h_rho_mcapprox_cospsi_helPos[isample] -> Fill(cosPsi);
	    //h_rho_mcapprox_costheta_helPos[isample] -> Fill(costh_exact);
	    //h_rho_mcapprox_cosbeta_helPos[isample] -> Fill(cosBeta);

	    h_rho_mcPolarExact_helPos[isample]->Fill(mcexactPolarimeter[i]);
	    h_rho_mcPolarApproxCheatEn_helPos[isample]->Fill(mcapproxCheatEnPolarimeter[i]);
	    h_rho_mcPolarApprox_helPos[isample]->Fill(mcapproxPolarimeter[i]);
	  }

	  h_pirho_tauMinusCosth_mcPolarExact[isample]->Fill(mctauMinusCosTh, mcexactPolarimeter[i] );

	}
      }

      if ( mcexactPolarimeter[0]>-1.99 && mcexactPolarimeter[1]>-1.99 ) {
	h_pirho_mcPolarExactPlusMinus[isample]->Fill(mcexactPolarimeter[0], mcexactPolarimeter[1]);
	h_pirho_mcPolarApproxPlusMinus[isample]->Fill(mcexactPolarimeter[0], mcexactPolarimeter[1]);
      }

    }

    // look for converted photons from pi0
    for (int j=0; j<mccol->getNumberOfElements(); j++) {
      MCParticle* mcp = dynamic_cast<MCParticle*> (mccol->getElementAt(j));

      if ( mcp->getPDG()==22 && mcp->getParents().size()>0 && mcp->getParents()[0]->getPDG()==111 ) {

	float radius = sqrt( pow(mcp->getEndpoint()[0],2) + pow (mcp->getEndpoint()[1], 2 ) );

	//	cout << "MC photon endpoint r, z " << radius << " " << mcp->getEndpoint()[2] << endl;

	h_conversionPos->Fill(fabs( mcp->getEndpoint()[2] ), radius);
	h_conversionPos2->Fill(fabs( mcp->getEndpoint()[2] ), radius);
	
	if ( fabs( mcp->getEndpoint()[2] ) < 2400 && radius < 1600 ) {
	  _nconvertedInTrk++;
	} else {
	  _nNotconvertedInTrk++;
	}


      }
    }


  } catch(DataNotAvailableException &e) {};


  if ( tautauInvMass>0 ) {
    h_ttmass->Fill(tautauInvMass);
    h_mcTau_costh[isample]->Fill(taucosth[0], taucosth[1]);
  }

  //  cout << tautauInvMass << " " << isample << endl;


  int nIsoElectrons(0);
  int nIsoMuons(0);
  try {
    LCCollection* isolepCol =  evt->getCollection( "ISOLeptons");
    if ( isolepCol->getNumberOfElements()>0 ) {
      IntVec ISOtypes;
      isolepCol->getParameters().getIntVals( "ISOLepType", ISOtypes );
      for (int i=0; i<isolepCol->getNumberOfElements(); i++) {
	if ( abs(ISOtypes[i]) == 11 ) nIsoElectrons++;
	else if ( abs(ISOtypes[i]) == 13 ) nIsoMuons++;
      }
    }
  } catch(DataNotAvailableException &e) {};
  h_nIsoMuons[isample]->Fill(nIsoMuons);
  h_nIsoElectrons[isample]->Fill(nIsoElectrons);




  bool verboseFF=false;

 if (verboseFF) cout << "hello 1" << endl;

  //-----------------------------------
  // apply simple preselection for tau pair events
  //-----------------------------------
  int nchpfo(0);
  int nneupfo(0);
  try {
    LCCollection* pfocol = _useDistilled ? evt->getCollection( "DistilledPFOs" ) : evt->getCollection( "PandoraPFOs") ;
    for (int j=0; j<pfocol->getNumberOfElements(); j++) {
      ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (pfocol->getElementAt(j));
      //      if ( fabs(pfo->getCharge())>0.1 ) {
      if ( isCharged(pfo) ) {
	nchpfo++;
      } else {
	nneupfo++;
      }
    }
  } catch(DataNotAvailableException &e) {};

  h_npfo_chg_nneu[isample]->Fill( nchpfo, nneupfo );


  if ( nchpfo < 2 || nchpfo>12 ) return;

  _nPresel[isample]++;



  //-----------------------------

  if ( _relNavi )  delete _relNavi;
  if ( _relNavi2 ) delete _relNavi2;

  try {
    LCCollection* linkcol = evt->getCollection( "RecoMCTruthLink" );
    _relNavi = new LCRelationNavigator(linkcol);
  } catch(DataNotAvailableException &e) {};

  try {
    LCCollection* linkcol = evt->getCollection( "MCTruthRecoLink" );
    _relNavi2 = new LCRelationNavigator(linkcol);
  } catch(DataNotAvailableException &e) {};


  //-----------------------------

  // look for electron and muon PFOs (using PID info)
  std::vector < ReconstructedParticle* > rp_electrons;
  std::vector < ReconstructedParticle* > rp_muons;
  PIDHandler* pidHandler(0);
  int myID ;
  try {
    LCCollection* pfocol = evt->getCollection( "PandoraPFOs");
    pidHandler=new PIDHandler(pfocol);
    int iLikePid =  pidHandler->getAlgorithmID( "LikelihoodPID" );
    for (int j=0; j<pfocol->getNumberOfElements(); j++) {
      ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (pfocol->getElementAt(j));
      ParticleIDVec likePids = pidHandler->getParticleIDs ( pfo, iLikePid );
      int likePdg = abs( pidHandler->getParticleID( pfo,  iLikePid ).getPDG() );
      if ( likePdg==11 ) {
	rp_electrons.push_back(pfo);
      } else if ( likePdg==13 ) {
	rp_muons.push_back(pfo);
      }
    }
  } catch(DataNotAvailableException &e) {};

  if ( pidHandler ) delete pidHandler;

  //-----------------------------

 if (verboseFF) cout << "hello 2" << endl;

  // try looking at the recosntructed v0 collection
  // possibly some conversions?
  try {
    LCCollection* v0_pfocol = evt->getCollection( "BuildUpVertex_V0_RP" );
    LCCollection* v0_vtxcol = evt->getCollection( "BuildUpVertex_V0" );
    if ( v0_pfocol->getNumberOfElements()>0 ) {
      //std::cout << "event " << evt->getEventNumber() << " number of built up vertices: " << v0_pfocol->getNumberOfElements() << std::endl;
      for (int j=0; j<v0_pfocol->getNumberOfElements(); j++) {
	ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (v0_pfocol->getElementAt(j));
	Vertex* vtx = dynamic_cast<Vertex*> (v0_vtxcol->getElementAt(j));

	h_v0Pos->Fill( fabs(vtx->getPosition()[2]) , sqrt( pow( vtx->getPosition()[0], 2 ) +  pow( vtx->getPosition()[1], 2 ) ) ); 
	h_v0Pos2->Fill( fabs(vtx->getPosition()[2]) , sqrt( pow( vtx->getPosition()[0], 2 ) +  pow( vtx->getPosition()[1], 2 ) ) ); 

	//cout << "PFO " << j << " type " << pfo->getType() << " mass " << pfo->getMass() << " en " << pfo->getEnergy() << " npart " << pfo->getParticles().size() << " " << endl;
	//cout << " vertex info " << vtx->getPosition()[0] << " " <<  vtx->getPosition()[1] << " " <<  vtx->getPosition()[2] << endl;
	for ( size_t jj=0; jj<pfo->getParticles().size(); jj++) {
	  ReconstructedParticle* pfo2 = pfo->getParticles()[jj];
          MCParticle* bestmatch = tauUtils::getBestTrackMatch( pfo2, _relNavi );
	  //cout << "   " << jj << " " << pfo2->getType() << " " << pfo2->getMass() << " " << pfo2->getEnergy() << " " << bestmatch;
	  //if (bestmatch) {
	    //cout << " MC match: " << bestmatch->getPDG() << " " << bestmatch->getEnergy();
	    //if ( bestmatch->getParents().size()>0 ) cout << " MC parent: " <<  bestmatch->getParents()[0]->getPDG();
	  //}

	  //  if ( find( rp_muons.begin(),  rp_muons.end(), pfo2 ) != rp_muons.end() ) cout << " LIKE_MUON " ;
	  //	  if ( find( rp_electrons.begin(),  rp_electrons.end(), pfo2 ) != rp_electrons.end() ) cout << " LIKE_ELEC " ;


	  //	  cout << endl;
	}
      }
    }
  } catch(DataNotAvailableException &e) {};

  // look for v0 pfos found by Pandora
  try {
    LCCollection* pfocol = evt->getCollection( "PandoraPFOs" );
    for (int j=0; j<pfocol->getNumberOfElements(); j++) {
      ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (pfocol->getElementAt(j));
      //      if ( fabs(pfo->getCharge())<0.1 ) {
      if ( isNeutral(pfo) ) {
	if ( pfo->getTracks().size()>0 || pfo->getParticles().size()>0 ) {
	  //cout << "neutral Pandora PFO with ntrks, nparts = " << 
	  //  pfo->getTracks().size() << " " << pfo->getParticles().size() << 
	  //  " type, mass = " << pfo->getType() << " " << pfo->getMass() << 
	  //  " vertex? " << pfo->getStartVertex() << endl;
	  h_neuPandoraCompound_mass->Fill(pfo->getMass());
	  h_neuPandoraCompound_ntrk->Fill(pfo->getTracks().size());
	  // causes crash on exit ? // 
	  vertexInfo* vtxInfo = new vertexInfo();
	  for ( size_t itrk=0; itrk< pfo->getTracks().size(); itrk++) {
	    vtxInfo->addTrack( pfo->getTracks()[itrk] );
	  }
	  const float *refPt = pfo->getTracks()[0]->getTrackState ( EVENT::TrackState::AtFirstHit ) -> getReferencePoint();
	  vtxInfo->setSeedPos( TVector3( refPt[0], refPt[1], refPt[2] ) );
	  TVector3 vtxpos = vtxInfo->getVertexPosition();
	  h_neuPandoraCompound_vtxPos->Fill(fabs(vtxpos[2]), sqrt( pow(vtxpos[0],2) + pow(vtxpos[1],2) ) );
	  h_neuPandoraCompound_vtxPos2->Fill(fabs(vtxpos[2]), sqrt( pow(vtxpos[0],2) + pow(vtxpos[1],2) ) );
	  h_neuPandoraCompound_vtxChisq->Fill( vtxInfo->getVertexChisq() );
	  
	  delete vtxInfo;
	}
      }
    }
  } catch(DataNotAvailableException &e) {};


 if (verboseFF) cout << "hello 3" << endl;

    
  try {
    LCCollection* pfocol = _useDistilled ? evt->getCollection( "DistilledPFOs" ) : evt->getCollection( "PandoraPFOs") ; 

    //identify which really come from tau
    std::vector < ReconstructedParticle* > PFOfromtau[2];
    std::vector < ReconstructedParticle* > noMatchPFOs;
    std::vector < ReconstructedParticle* > underlyingPFOs;

    for (int j=0; j<pfocol->getNumberOfElements(); j++) {
      ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (pfocol->getElementAt(j));

      std::vector < ReconstructedParticle* > subpfos;
      subpfos.push_back( pfo );
      if ( pfo->getParticles().size()>0 ) {
	for ( size_t kk=0; kk<pfo->getParticles().size(); kk++ ) {
	  subpfos.push_back( pfo->getParticles()[kk] );
	}
      }
      //      cout << "matching " << pfo->getType() << " " << pfo->getEnergy() << " nsub " << subpfos.size() << endl;

      std::vector < MCParticle* > mcpmatch;
      for ( size_t kk=0; kk<subpfos.size(); kk++ ) {
	//	cout << " subpfo " << kk << " charge, en " <<  subpfos[kk]->getCharge() << " " << subpfos[kk]->getEnergy() << endl;
	if ( isCharged( subpfos[kk] ) ) {
	  MCParticle* mcp = tauUtils::getBestTrackMatch ( subpfos[kk], _relNavi );
	  if ( mcp ) mcpmatch.push_back(mcp);
	} else {
	  MCParticle* mcp = tauUtils::getBestCaloMatch ( subpfos[kk], _relNavi );
	  if ( mcp ) mcpmatch.push_back(mcp);
	}
      }


//      // let's look at some PFOs in more detail
//
//      if ( pfo->getCharge()==0 ) {
//
//	cout << "neutral PFO of type " << pfo->getType() << " energy " << pfo->getEnergy() << endl;
//
//	cout << "reco->mc related MCparticles: ";
//
//	std::map < MCParticle*, std::pair < float, float > > conts;
//
//	for (size_t k=0; k<_relNavi->getRelatedToObjects( pfo ).size(); k++) {
//	  MCParticle* mcp = dynamic_cast <MCParticle*> (_relNavi->getRelatedToObjects(pfo)[k]);
//	  float wgt = _relNavi->getRelatedToWeights(pfo)[k];
//	  float clusterwgt = (int(wgt)/10000)/1000. ;
//
//	  int parPdg = mcp->getParents().size()==0 ? 0 : mcp->getParents()[0]->getPDG();
//
//	  cout << mcp->getPDG() << " " << mcp->getEnergy() << " " << parPdg << " " << clusterwgt << " , ";
//
//	  conts[mcp]=std::pair < float, float > ( clusterwgt, 0. );
//
//	}
//	cout << endl;
//
//	cout << "mc->reco related MCparticles: ";
//	for (size_t k=0; k<_relNavi2->getRelatedFromObjects( pfo ).size(); k++) {
//	  MCParticle* mcp = dynamic_cast <MCParticle*> (_relNavi2->getRelatedFromObjects(pfo)[k]);
//	  float wgt = _relNavi2->getRelatedFromWeights(pfo)[k];
//	  float clusterwgt = (int(wgt)/10000)/1000. ;
//	  cout << mcp->getPDG() << " " << mcp->getEnergy() << " " << clusterwgt << " , ";
//
//	  conts[mcp].second=clusterwgt;
//
//	}
//	cout << endl;
//
//	int majorContributors(0);
//	int significantContributors(0);
//	int significantPhotonContributors(0);
//	for (auto& kv : conts) {
//
//	  std::cout << kv.first->getPDG() << " " << kv.first->getEnergy() << " has value " << kv.second.first << " " << kv.second.second << std::endl;
//
//	  if (  kv.second.first > 0.1 &&  // this MC contributes to at least 10% of cluster energy
//		kv.second.second > 0.5    // at least 50% of this MC's calo energy deposit is in this cluster
//		) {
//	    majorContributors++;
//	  }
//
//	  if ( kv.second.second > 0.5 ) {   // at least 50% of this MC's calo energy deposit is in this cluster
//	    significantContributors++;
//
//	    if (  kv.first->getPDG()==22 ) {
//	      significantPhotonContributors++;
//	    }
//
//	  }
//	}
//
//
//	cout << "major/significant contributors " << majorContributors << " " << significantContributors << endl;
//
//	// try to get Bulos mass
//	//	cout << "nClusters = " << pfo->getClusters().size() << endl;
//	for (size_t k=0; k< pfo->getClusters().size(); k++) {
//	  Cluster* cl = pfo->getClusters()[k];
////	  cout << "cluster " << k << " shape : " ;
////	  for (size_t kk=0; kk< cl->getShape().size(); kk++) {
////	    cout << cl->getShape()[kk] << " ";
////	  }
////	  cout << endl;
////
////	  cout << "position error: ";
////	  for (int kk=0; kk<6; kk++) {
////	    cout << cl->getPositionError()[kk] << " ";
////	  }
////	  cout << endl;
//
//	  double temp[3][3];
//	  temp[0][0]=cl->getPositionError()[0];
//	  temp[0][1]=cl->getPositionError()[1];
//	  temp[0][2]=cl->getPositionError()[2];
//	  temp[1][1]=cl->getPositionError()[3];
//	  temp[1][2]=cl->getPositionError()[4];
//	  temp[2][2]=cl->getPositionError()[5];
//	  temp[1][0]=temp[0][1];
//	  temp[2][0]=temp[0][2];
//	  temp[2][1]=temp[1][2];
//
//	  TMatrixDSym errmat(3, &(temp[0][0]) );
//
//	  //cout << errmat[0][0] << " " << errmat[0][1] << " " << errmat[0][2] << endl;
//	  //	  cout << errmat[1][0] << " " << errmat[1][1] << " " << errmat[1][2] << endl;
//	  //	  cout << errmat[2][0] << " " << errmat[2][1] << " " << errmat[2][2] << endl;
//
//	  TMatrixD eigenVectors(3,3);
//	  TVectorD eigenValues(3);
//	  eigenVectors = errmat.EigenVectors( eigenValues );
//
//	  std::vector <float> evals;
//	  for (int kk=0; kk<3; kk++) {
//	    //	    cout << "eigen vec #" << kk << " : " << eigenVectors[kk][0] << " " << eigenVectors[kk][1] << " " << eigenVectors[kk][2] << " : " << eigenValues[kk] << endl;
//	    evals.push_back( fabs(eigenValues[kk]) );
//	  }
//	  std::sort(evals.begin(), evals.end());
//	  //	  cout << evals[0] << " " << evals[1] << " " << evals[2] << endl;
//
//	  float clusterDist = sqrt( pow( cl->getPosition()[0], 2 ) + 
//				    pow( cl->getPosition()[1], 2 ) + 
//				    pow( cl->getPosition()[2], 2 ) );
//
//	  float bulosDist = sqrt( pow(evals[1],2) - pow(evals[0],2) );
//	  float bulosAngle = bulosDist/clusterDist;
//
//	  float bulosMass = cl->getEnergy()*bulosAngle/2.;
//
//	  //	  cout << "bulos mass: " << bulosMass << endl;
//
//
//	  if (pfo->getType()==22) {
//	    h_gammaClus_nMajor->Fill(majorContributors);
//	    h_gammaClus_nSignf->Fill(significantContributors);
//	    h_gammaClus_nSignf_bMass->Fill( significantPhotonContributors, bulosMass );
//	    h_gammaClus_nSignf_bDist->Fill( significantPhotonContributors, bulosDist );
//	  } else if  (pfo->getType()==2112) {
//	    h_neutronClus_nMajor->Fill(majorContributors);
//	    h_neutronClus_nSignf->Fill(significantContributors);
//	    h_neutronClus_nSignf_bMass->Fill( significantPhotonContributors, bulosMass );
//	    h_neutronClus_nSignf_bDist->Fill( significantPhotonContributors, bulosDist );
//	  }
//
//
//	}
//
//
//      }


      //      for ( size_t k=0 ; k<mcpmatch.size(); k++) {
      //	cout << k << " " << mcpmatch[k]->getPDG() << " " <<  mcpmatch[k]->getEnergy() << endl;
      //      }
	  

      bool taumatched(false);
      if ( mcpmatch.size()==0 ) {
	noMatchPFOs.push_back( pfo );
      } else {
	for (int itau=0; itau<2; itau++) {
	  for (int jj=0; jj<mcpmatch.size(); jj++) {
	    if ( std::find( allMCtaudaughters[itau].begin(), allMCtaudaughters[itau].end(), mcpmatch[jj] ) != allMCtaudaughters[itau].end() ) {
	      PFOfromtau[itau].push_back( pfo );
	      taumatched=true;
	      break;
	    }
	  }
	  if ( taumatched ) break;
	}
	if (!taumatched) underlyingPFOs.push_back( pfo );
      }
    }

    // first find the best seed prongs
    std::pair< float , ReconstructedParticle* > highestPtChargedPFO[2];
    highestPtChargedPFO[0]=std::pair< float , ReconstructedParticle* >(-1, NULL);
    highestPtChargedPFO[1]=std::pair< float , ReconstructedParticle* >(-1, NULL);

    // first highest pt track
    for (int j=0; j<pfocol->getNumberOfElements(); j++) {
      ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (pfocol->getElementAt(j));
      //      if ( fabs(pfo->getCharge())>0.1 ) {
      if ( isCharged(pfo) ) {

	float pt = _highestPt ? 
	  sqrt( pow( pfo->getMomentum()[0], 2 ) + pow( pfo->getMomentum()[1], 2 ) ) : // pt
	  sqrt( pow( pfo->getMomentum()[0], 2 ) + pow( pfo->getMomentum()[1], 2 ) +  pow( pfo->getMomentum()[2], 2 ) ); //  momentum


	if ( pt > highestPtChargedPFO[0].first ) {
	  highestPtChargedPFO[0] = std::pair< float , ReconstructedParticle* > ( pt, pfo );
	} 
      }
    }


 if (verboseFF) cout << "hello 4" << endl;


    //if ( ! highestPtChargedPFO[0].second ) {
    //  cout << "did not find first seed track!" << endl;
    //} else {
    if ( highestPtChargedPFO[0].second ) {
      // find second cone seed
      TVector3 coneSeedDir[2];
      coneSeedDir[0] = tauUtils::getTLV(highestPtChargedPFO[0].second).Vect();
      // then highest pt at least pi/2 away from first in deltaPhi
      for (int j=0; j<pfocol->getNumberOfElements(); j++) {
	ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (pfocol->getElementAt(j));
	TVector3 mom( pfo->getMomentum()[0],pfo->getMomentum()[1],pfo->getMomentum()[2]);
	if ( isCharged(pfo) && fabs(mom.DeltaPhi(coneSeedDir[0]))>TMath::Pi()/2. ) {
	  float pt = _highestPt ? mom.Pt() : mom.Mag();
	  if ( pt > highestPtChargedPFO[1].first ) {
	    highestPtChargedPFO[1] = std::pair< float , ReconstructedParticle* > ( pt, pfo );
	  } 
	}
      }

      // if ( ! highestPtChargedPFO[1].second ) {
      // 	// find second cone seed
      // 	// then highest pt at least pi/2 away from first in deltaPhi
      // 	for (int j=0; j<pfocol->getNumberOfElements(); j++) {
      // 	  ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (pfocol->getElementAt(j));
      // 	  TVector3 mom( pfo->getMomentum()[0],pfo->getMomentum()[1],pfo->getMomentum()[2]);
      // 	}
      // }

      if ( highestPtChargedPFO[1].second ) {
	_nTwoSeeds[isample]++;

	float seedecalen[2]={0,0};
	float seedcaloen[2]={0,0};
	float seedClWidth1[2]={0,0};
	float seedClWidth2[2]={0,0};
	float seedClLength[2]={0,0};
	float seedecaleonp[2]={0,0};
	for (int ij=0; ij<2; ij++) {
	  float ecalen(0);
	  float caloen(0);
	  ClusterVec clv = highestPtChargedPFO[ij].second->getClusters();
	  Cluster* maxCl(0);
	  for (int ic=0; ic<clv.size(); ic++) {
	    Cluster* cl = clv[ic];
	    caloen+=cl->getEnergy();
	    ecalen+=cl->getSubdetectorEnergies()[0];
	    ecalen+=cl->getSubdetectorEnergies()[3]; // lumical
	    ecalen+=cl->getSubdetectorEnergies()[5]; // beamcal

	    if ( !maxCl || cl->getEnergy() > maxCl->getEnergy() ) {
	      maxCl=cl;
	    }

	  }
	  seedcaloen[ij]=caloen;
	  seedecalen[ij]=ecalen;


	  if ( maxCl ) {
	    std::vector < float > clevals = tauUtils::getClusterEigenvalues( maxCl );
	    seedClWidth1[ij]=fabs( clevals[0] );
	    seedClWidth2[ij]=fabs( clevals[1] );
	    seedClLength[ij]=fabs( clevals[2] );
	  } else {
	    seedClWidth1[ij]=0;
	    seedClWidth2[ij]=0;
	    seedClLength[ij]=0;
	  }
	  

	  float mom(0);
	  for (int i=0; i<3; i++) 
	    mom+=pow(highestPtChargedPFO[ij].second->getMomentum()[i], 2 );
	  mom=sqrt(mom);
	  seedecaleonp[ij]=ecalen / mom;
	}

	h_maxSeedEn[isample]->Fill( highestPtChargedPFO[0].first );
	h_minSeedEn[isample]->Fill( highestPtChargedPFO[1].first );

	h_maxSeedEn_eonp[isample]->Fill( seedecaleonp[0] );
	h_minSeedEn_eonp[isample]->Fill( seedecaleonp[1] );

	h_maxSeedEn_caloen[isample]->Fill( seedcaloen[0] );
	h_minSeedEn_caloen[isample]->Fill( seedcaloen[1] );


	// cout << "a " << highestPtChargedPFO[0].second->getTracks().size() << endl;
	// cout << "b " << highestPtChargedPFO[0].second->getTracks()[0]->getZ0() << endl;
	// cout << "c " << highestPtChargedPFO[1].second->getTracks().size() << endl;
	// cout << "d " << highestPtChargedPFO[1].second->getTracks()[0]->getZ0() << endl;

	h_seed_dZ0[isample]->Fill( highestPtChargedPFO[0].second->getTracks()[0]->getZ0() - pvposMC.Z() );
	h_seed_dZ0[isample]->Fill( highestPtChargedPFO[1].second->getTracks()[0]->getZ0() - pvposMC.Z() );
	h_seed_dZ0ave[isample]->Fill( ( highestPtChargedPFO[0].second->getTracks()[0]->getZ0() + highestPtChargedPFO[1].second->getTracks()[0]->getZ0() )/2. - pvposMC.Z()  );

	//	cout << "blah" << endl;

	float costh0= fabs( TVector3(  highestPtChargedPFO[0].second->getMomentum() ).CosTheta() );
	float costh1= fabs( TVector3(  highestPtChargedPFO[1].second->getMomentum() ).CosTheta() );

	float costhMax= costh0 > costh1 ? costh0 : costh1;
	float costhMin= costh0 > costh1 ? costh1 : costh0;

	h_maxSeedCosth[isample]->Fill( costhMax );
        h_minSeedCosth[isample]->Fill( costhMin );


	h_seedEcalEn[isample]->Fill( seedecalen[0] );
	h_seedEcalEn[isample]->Fill( seedecalen[1] );

	h_seedClusterWidth1[isample]->Fill( log10(seedClWidth1[0]) );
	h_seedClusterWidth1[isample]->Fill( log10(seedClWidth1[1]) );
	h_seedClusterWidth2[isample]->Fill( log10(seedClWidth2[0]) );
	h_seedClusterWidth2[isample]->Fill( log10(seedClWidth2[1]) );

	

	h_seedClusterLength[isample]->Fill( log10(seedClLength[0]) );
	h_seedClusterLength[isample]->Fill( log10(seedClLength[1]) );

	h_seedEcalEonP[isample]->Fill( seedecaleonp[0] );
	h_seedEcalEonP[isample]->Fill( seedecaleonp[1] );
	  

	coneSeedDir[1] = tauUtils::getTLV(highestPtChargedPFO[1].second).Vect();
	float prongangle =  TMath::Pi() - tauUtils::getTLV( highestPtChargedPFO[0].second ).Angle( tauUtils::getTLV( highestPtChargedPFO[1].second ).Vect() );
	// try to associate seeds with MC tau directions
	MCParticle* matchedTauByDir[2]={0,0};
	int matchedTauByDirDecay[2]={-1,-1};
	int matchedTauByDirIndex[2]={-1,-1};
	if ( finalmctaus.size()==2 ) {
	  for (int jseed=0; jseed<2; jseed++) { // the two seeds
	    for (int itau=0; itau<2; itau++) {
	      float tauSeedAngle =  tauUtils::getTLV( finalmctaus[itau] ).Angle( coneSeedDir[jseed] );
	      if ( tauSeedAngle < _conesize ) {
		matchedTauByDir[jseed]=finalmctaus[itau];
		matchedTauByDirDecay[jseed]=tauUtils::getMCdecayMode(matchedTauByDir[jseed]);
		matchedTauByDirIndex[jseed]=itau;
	      }
	    }
	  }
	  assert( !matchedTauByDir[0] || matchedTauByDir[0] != matchedTauByDir[1] );
	}

	bool goodSeedDir = matchedTauByDir[0] && matchedTauByDir[1];
	if (goodSeedDir) {
	  _nGoodSeedDir[isample]++;
	}

	// check if the chosen seed tracks are really from tau
	int ttre[2]={-1,-1};
	for (int jseed=0; jseed<2; jseed++) { // the two seeds
	  ReconstructedParticle* pfo = highestPtChargedPFO[jseed].second;
	  MCParticle* bestmatch = tauUtils::getBestTrackMatch( pfo, _relNavi );
	  if ( !bestmatch ) {
	    cout << "could not find MCP match..." << pfo->getType() << " " << pfo->getEnergy() << " " << _relNavi->getRelatedToWeights(pfo).size() << endl;
	  } else {
	    // which tau does this MC particle come from?
	    for (int jtau=0; jtau<2; jtau++) {
	      if ( find( allMCtaudaughters[jtau].begin(), allMCtaudaughters[jtau].end(), bestmatch ) != allMCtaudaughters[jtau].end() ) {
		ttre[jseed]=jtau;
	      }
	    }
	  }
	}
	bool goodMCMatchSeeds = ttre[0]>=0 && ttre[1]>=0 && ttre[0]!=ttre[1];
	if (goodMCMatchSeeds) {
	  _nTwoWellMatchedSeeds[isample]++;
	}
	// then gather particles within cone around seeds
	std::vector < ReconstructedParticle* > coneparticles[2];
	std::vector < ReconstructedParticle* > outsideconeparticles;
	std::vector < ReconstructedParticle* > outsideconephotons;
	for (int j=0; j<pfocol->getNumberOfElements(); j++) {
	  ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (pfocol->getElementAt(j));
	  TLorentzVector ttr = tauUtils::getTLV(pfo);
	  float angle0=ttr.Vect().Angle( coneSeedDir[0] );
	  float angle1=ttr.Vect().Angle( coneSeedDir[1] );
	  if ( angle0<_conesize ) {
	    coneparticles[0].push_back( pfo );
	  } else if ( angle1<_conesize ) {
	    coneparticles[1].push_back( pfo );
	  } else {
	    // cout << "outdside cone: " << pfo->getType() << " " << pfo->getEnergy() << " " << ttr.CosTheta() << " " << angle0 << " " << angle1 << endl;
	    outsideconeparticles.push_back( pfo );
	    if ( pfo->getType()==22 ) 
	      outsideconephotons.push_back( pfo );
	  }
	}

 if (verboseFF) cout << "hello 6" << endl;


	TLorentzVector jet4mom[2];
	TLorentzVector jetNeutral4mom[2];
	int cone_chg[2]={0,0};
	int cone_ngam[2]={0,0};
	int cone_npi0[2]={0,0};
	int cone_nchg[2]={0,0};
	int cone_nnhad[2]={0,0};

	int nTrkInCone[2]={0,0};

	std::vector < Track* > tracksInCone[2];

	for (int j=0; j<2; j++) {
	  //	  if ( isample==0 ) { // extra checks for multi-prongs
	  // int decay=-1;
	  // if ( matchedTauByDir[j] ) {
	  //   decay =  tauUtils::getMCdecayMode(  matchedTauByDir[j]  );
	  // }
	  int decay = matchedTauByDirDecay[j];


	  try {
	    LCCollection* trkcol = evt->getCollection( "MarlinTrkTracks" );
	    for (int jj=0; jj<trkcol->getNumberOfElements(); jj++) {
	      Track* trk = dynamic_cast<Track*> (trkcol->getElementAt(jj));
	      TLorentzVector trktlv = tauUtils::getFourMomentum( trk->getTrackState( TrackState::AtIP ), 0.139 );
	      float angle=trktlv.Angle( coneSeedDir[j] );
	      if ( angle < _conesize ) {
		nTrkInCone[j]++;
		tracksInCone[j].push_back(trk);
	      }
	    }
	  } catch(DataNotAvailableException &e) {};




	  //if ( decay == tauUtils::decayA1_3p ) {
	  //  TVector3 tauDir =  tauUtils::getTLV( matchedTauByDir[j] ).Vect();
	  //  int nmcInCone(0);
	  //  std::vector < EVENT::MCParticle* > stmc = tauUtils::getstablemctauDaughters( matchedTauByDir[j] );
	  //  for (size_t k=0; k<stmc.size(); k++) {
	  //    float angle=tauUtils::getTLV(  stmc[k] ).Angle( tauDir );
	  //    if ( isCharged(stmc[k]) && angle<_conesize ) nmcInCone++;
	  //  }
	  //  int npfoInCone(0);
	  //  for (int jj=0; jj<pfocol->getNumberOfElements(); jj++) {
	  //    ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (pfocol->getElementAt(jj));
	  //    if ( isCharged(pfo)) {
	  //	float angle = tauUtils::getTLV( pfo ).Angle( tauDir );
	  //	if ( angle<_conesize ) npfoInCone++;
	  //    }
	  //  }

	  //  h_a3p_cone_chgMC_pfo[isample]->Fill(nmcInCone,npfoInCone);
	  //  h_a3p_cone_chgMC_trk[isample]->Fill(nmcInCone,nTrkInCone);
	  //  h_a3p_cone_trk_pfo  [isample]->Fill(nTrkInCone,npfoInCone);
	  //}

	  int ngam(0);
	  float engam(0);
	  int nhad(0);
	  float enhad(0);
	  int nch(0);
	  float ench(0);
	  int ncompound(0);
	  int npi0(0);

	  for ( size_t i=0; i<coneparticles[j].size(); i++) {
	    if ( isCharged ( coneparticles[j][i] ) ) {
	      nch++;
	      ench+=coneparticles[j][i]->getEnergy();
	    } else {
	      if ( coneparticles[j][i]->getType()==22 ) {
		ngam++;
		engam+=coneparticles[j][i]->getEnergy();
	      } else if ( coneparticles[j][i]->getType()==111 ) {
		npi0++;
		ngam+=2;
		engam+=coneparticles[j][i]->getEnergy();
	      } else if  ( coneparticles[j][i]->getType()==2112 ) {
		nhad++;
		enhad+=coneparticles[j][i]->getEnergy();
	      } else if (  coneparticles[j][i]->getType()==221 || coneparticles[j][i]->getType()==310 || abs(coneparticles[j][i]->getType())==3122 || coneparticles[j][i]->getType()==331 ) {
		ncompound++;
	      } else {
		cout << "unknown PFO type !!!" <<  coneparticles[j][i]->getType() << endl;
		assert(0);
	      }
	    }
	  }

	  jet4mom[j].SetXYZT(0,0,0,0);
	  jetNeutral4mom[j].SetXYZT(0,0,0,0);
	  for ( size_t i=0; i<coneparticles[j].size(); i++) {
	    TLorentzVector thisTLV = tauUtils::getTLV( coneparticles[j][i] );
	    jet4mom[j] += thisTLV;

	    if ( isNeutral( coneparticles[j][i] ) ) {
	      jetNeutral4mom[j] += thisTLV;
	    }

	    cone_chg[j]+=coneparticles[j][i]->getCharge();
	    if ( coneparticles[j][i]->getType()==22 )       cone_ngam[j]++;
	    else if ( coneparticles[j][i]->getType()==111 ) cone_ngam[j]+=2;
	    else if ( isCharged( coneparticles[j][i] ) ) cone_nchg[j]++;
	    if  ( coneparticles[j][i]->getType()==2112 ) cone_nnhad[j]++;
	  }

	  int nggMCtaudau(0);
	  int nggMCtaudauCone=0;
	  int nggMCallCone=0;
	  if ( matchedTauByDirIndex[j]>=0 ) {
	    for ( auto m : stableMCtaudaughters[matchedTauByDirIndex[j]] ) {
	      if ( m->getPDG()==22 ) nggMCtaudau++;
	    }
	    for ( auto m : stableMCtaudaughtersInCone[matchedTauByDirIndex[j]] ) {
              if ( m->getPDG()==22 ) nggMCtaudauCone++;
            }
	    for ( auto m : stableMCallInCone[matchedTauByDirIndex[j]] ) {
              if ( m->getPDG()==22 ) nggMCallCone++;
            }
	  }

	  for (int jj=0; jj<2; jj++) {

	    int idec(-1);

	    if (jj==0) {
	      idec=NDEC-1; // all events
	    } else if ( decay>=0 ) {
	      idec= std::min( decay, NDEC-2 );
	      // decay<6 ? decay : NDEC-2;
	    }

	    if ( idec>=0 ) {

	      h_dec_seed_clusterWidth1[isample][idec]->Fill( log10(seedClWidth1[j]) );
	      h_dec_seed_clusterWidth2[isample][idec]->Fill( log10(seedClWidth2[j]) );
	      h_dec_seed_clusterLength[isample][idec]->Fill( log10(seedClLength[j]) );
	      h_dec_seed_energy[isample][idec]->Fill( highestPtChargedPFO[j].second->getEnergy() );

	      h_dec_cone_nchpfo[isample][idec]->Fill(nch);
	      h_dec_cone_ntracks[isample][idec]->Fill(nTrkInCone[j]);

	      for ( auto tt : tracksInCone[j] ) {
		h_dec_cone_trackNTPC[isample][idec]->Fill(tt->getSubdetectorHitNumbers()[6]);
	      }


	      h_dec_cone_nnhadpfo[isample][idec]->Fill(nhad);
	      h_dec_cone_ngammapfo[isample][idec]->Fill(ngam);
	      h_dec_cone_npi0pfo[isample][idec]->Fill(npi0);
	      h_dec_cone_ncompoundpfo[isample][idec]->Fill(ncompound);
	      h_dec_cone_nhaden[isample][idec]->Fill(enhad);
	      h_dec_cone_gammaen[isample][idec]->Fill(engam);
	      h_dec_cone_nhadenFrac[isample][idec]->Fill(enhad/(enhad+engam+ench));
	      h_dec_cone_gammaenFrac[isample][idec]->Fill(engam/(enhad+engam+ench));
	      
	      h_dec_cone_visMass[isample][idec]->Fill( jet4mom[j].M() );
	      h_dec_cone_neutralvisMass[isample][idec]->Fill( jetNeutral4mom[j].M() );

	      h_dec_cone_vis_neutral_Mass[isample][idec]->Fill( jet4mom[j].M(), jetNeutral4mom[j].M() );
	      if ( ngam==0) h_dec_cone_vis_neutral_Mass_0g[isample][idec]->Fill( jet4mom[j].M(), jetNeutral4mom[j].M() );
	      else if ( ngam==1) h_dec_cone_vis_neutral_Mass_1g[isample][idec]->Fill( jet4mom[j].M(), jetNeutral4mom[j].M() );
	      else if ( ngam==2) h_dec_cone_vis_neutral_Mass_2g[isample][idec]->Fill( jet4mom[j].M(), jetNeutral4mom[j].M() );
	      else if ( ngam==3) h_dec_cone_vis_neutral_Mass_3g[isample][idec]->Fill( jet4mom[j].M(), jetNeutral4mom[j].M() );
	      else if ( ngam>=4) h_dec_cone_vis_neutral_Mass_4g[isample][idec]->Fill( jet4mom[j].M(), jetNeutral4mom[j].M() );
	      



	      h_dec_cone_visMassDiff[isample][idec]->Fill( jet4mom[j].M() - mcTauVis4mom[matchedTauByDirIndex[j]].M() );
	      h_dec_cone_neutralvisMassDiff[isample][idec]->Fill( jetNeutral4mom[j].M()  - mcTauNeutralVis4mom[matchedTauByDirIndex[j]].M() );
	      
	      h_dec_cone_visEnergyDiff[isample][idec]->Fill( jet4mom[j].E() - mcTauVis4mom[matchedTauByDirIndex[j]].E() );
	      h_dec_cone_neutralvisEnergyDiff[isample][idec]->Fill( jetNeutral4mom[j].E()  - mcTauNeutralVis4mom[matchedTauByDirIndex[j]].E() );


	      h_dec_mc_ngamma               [isample][idec]->Fill( nggMCtaudau );
	      h_dec_mc_vismass              [isample][idec]->Fill( mcTauVis4mom[matchedTauByDirIndex[j]].M() );
	      h_dec_mc_visneutralmass       [isample][idec]->Fill( mcTauNeutralVis4mom[matchedTauByDirIndex[j]].M() );

	      h_dec_mc_cone_ngamma          [isample][idec]->Fill( nggMCtaudauCone );
	      h_dec_mc_cone_vismass         [isample][idec]->Fill( tauUtils::getTLV( stableMCtaudaughtersInCone[matchedTauByDirIndex[j]] ).M() );
	      h_dec_mc_cone_visneutralmass  [isample][idec]->Fill( tauUtils::getTLV( stableMCneutraltaudaughtersInCone[matchedTauByDirIndex[j]] ).M() );

	      h_dec_mcall_cone_ngamma          [isample][idec]->Fill( nggMCallCone );
	      h_dec_mcall_cone_vismass         [isample][idec]->Fill( tauUtils::getTLV( stableMCallInCone[matchedTauByDirIndex[j]] ).M() );
	      h_dec_mcall_cone_visneutralmass  [isample][idec]->Fill( tauUtils::getTLV( stableMCallNeutralInCone[matchedTauByDirIndex[j]] ).M() );
	    }


	  }


	} // 2 cone loop

	// should look for singleton photons outside cone, which make pi0 mass with single photon inside cone and dont take jet mass above mtau
	for ( size_t i=0; i<outsideconephotons.size(); i++) {
	  TLorentzVector outsidePhotonTlv = tauUtils::getTLV( outsideconephotons[i] );
	  for ( int j=0; j<2; j++) {
	    for ( size_t k=0; k<coneparticles[j].size(); k++) {
	      if ( coneparticles[j][k]->getType()==22 ) {
		TLorentzVector insidePhotonTlv = tauUtils::getTLV( coneparticles[j][k] );
		h_gamgam_mass->Fill( (outsidePhotonTlv+insidePhotonTlv).M() );
	      }
	    }
	  }
	} // photons outside cone
	// this seems not so useful (most pi0 already reconstructed in distilled collection?)

 if (verboseFF) cout << "hello 6.4" << endl;

	//------------------------------------------------
	// do a first selection
	//------------------------------------------------

	// look to remove di-muon di-electron events
	// if either cone has a single electron or muon with the beam energy, tag as di-lepton event
	// bool mumu=false;
	// bool elel=false;
	// 
	// int nch[2]={0,0};
	// ReconstructedParticle* rpch[2]={0,0};
	// 
	// for (int icone=0; icone<2; icone++) {
	//   for ( size_t ip=0; ip<coneparticles[icone].size(); ip++) {
	//     ReconstructedParticle* rp = coneparticles[icone][ip];
	//     if ( fabs( rp->getCharge() )>0.1 ) {
	//       rpch[icone]=rp;
	//       nch[icone]++;
	//     }
	//   }
	// }
	// 
	// if ( nch[0]==1 && nch[1]==1 ) {
	//   // same flavour lepton in each cone
	//   //    nb could also ask them to be prompt -> keep events with 2 leptonic tau decays
	//   if ( find( rp_muons.begin(), rp_muons.end(), rpch[0] )  != rp_muons.end() &&
	//        find( rp_muons.begin(), rp_muons.end(), rpch[1] )  != rp_muons.end() ) {
	//     mumu=true;
	//   }
	//   if ( find( rp_electrons.begin(), rp_electrons.end(), rpch[0] ) != rp_electrons.end() &&
	//        find( rp_electrons.begin(), rp_electrons.end(), rpch[1] ) != rp_electrons.end() ) {
	//     elel=true;
	//   }
	//   // or both charged have high energy, and at least one is ID as lepton   (high energy events)
	//   if ( rpch[0]->getEnergy()>240 && rpch[1]->getEnergy()>240 ) {
	//     if (  find( rp_muons.begin(), rp_muons.end(), rpch[0] )  != rp_muons.end() ||
	// 	  find( rp_muons.begin(), rp_muons.end(), rpch[1] )  != rp_muons.end() ) {
	//       mumu=true;
	//     }
	//     if ( find( rp_electrons.begin(), rp_electrons.end(), rpch[0] ) != rp_electrons.end() ||
	// 	 find( rp_electrons.begin(), rp_electrons.end(), rpch[1] ) != rp_electrons.end() ) {
	//       elel=true;
	//     }
	//   }
	// }
	// 
	// // check that neither seed looks like high energy electron
	// //for ( int ij=0; ij<2; ij++) {
	// //  if ( highestPtChargedPFO[ij].first>200 && seedecaleonp[ij]>0.85 ) elel=true;
	// //	}


	// check energy outside cone
	float outsideConeEnergy(0);
	float outsideConePt(0);
	float maxoutsidepfoen(0);
	for (size_t i=0; i<outsideconeparticles.size(); i++) {
	  TLorentzVector outTlv = tauUtils::getTLV( outsideconeparticles[i] );
	  outsideConeEnergy+=outTlv.E();
	  outsideConePt+=outTlv.Pt();
	  if ( outsideconeparticles[i]->getEnergy() > maxoutsidepfoen ) {
	    maxoutsidepfoen=outsideconeparticles[i]->getEnergy();
	  }
	}

	// remove events with high energy visible photon outside the cones
	bool largeISR=false;
	float maxisr(0);
	for (size_t ii=0 ; ii<outsideconephotons.size(); ii++) {
	  if ( outsideconephotons[ii]->getEnergy() > maxisr ) {
	    maxisr=outsideconephotons[ii]->getEnergy();
	  }

	  if ( outsideconephotons[ii]->getEnergy()>10 ) {
	    largeISR=true;
	    break;
	  }
	}
	float insideConeEnergy = jet4mom[0].E() + jet4mom[1].E() ;

	float visjetacolinearity = TMath::Pi() - jet4mom[0].Vect().Angle(jet4mom[1].Vect());
	float visjetacoplanarity = TMath::Pi() - fabs(jet4mom[0].Vect().DeltaPhi(jet4mom[1].Vect()));

	h_ooconeMaxGammaEn[isample]->Fill( maxisr );
	h_ooconeMaxPFOEn[isample]->Fill( maxoutsidepfoen );

	h_visjetAcolinearity[isample]  ->Fill( visjetacolinearity );
	h_visjetAcoplanarity [isample] ->Fill( visjetacoplanarity );

	h_prongAngle[isample]   ->Fill( prongangle );


	h_outsideEnergy[isample]->Fill( outsideConeEnergy );
	h_insideEnergy[isample] ->Fill( insideConeEnergy );
	h_outsidePt[isample]    ->Fill( outsideConePt );


	for (int j=0; j<2; j++) {
	  h_coneMass     [isample]->Fill( jet4mom[j].M() );
	}
	h_seedSeedAngle[isample]->Fill(  TMath::Pi() - coneSeedDir[0].Angle(coneSeedDir[1])  );
	h_seedSeedDphi [isample]->Fill(  TMath::Pi() - fabs(coneSeedDir[0].DeltaPhi(coneSeedDir[1]))  );



	//bool prongAngleCut = prongangle < 0.07;
	//bool deltaPhiCut =  fabs(coneSeedDir[0].DeltaPhi(coneSeedDir[1])) > 3.1;

	bool acoPlanCut = visjetacolinearity < 0.05;
	bool acoLinCut = visjetacoplanarity < 0.075;

	bool outsideEnergyCut = outsideConeEnergy < 40. && outsideConePt < 20.;
	//	bool dilepCut = !mumu && !elel;
	bool leptonCut =  nIsoElectrons==0 && nIsoMuons == 0;
	bool isrCut = !largeISR;
	bool secondSeedEnergyCut = highestPtChargedPFO[1].first < 200.;

	// bool seedCaloEnCut =  
	//   seedecaleonp[0]<0.85 && seedecaleonp[1]<0.85 && 
	//   ( highestPtChargedPFO[0].first < 25. || seedcaloen[0]>15 ) && 
	//   ( highestPtChargedPFO[1].first < 25. || seedcaloen[1]>15 ) ;

	bool seedClusterCut = 
	  seedcaloen[1] > 5. && 
	  log10(seedClWidth2[0])>1. && log10(seedClWidth2[0])<3.6 &&
	  log10(seedClWidth2[1])>1. && log10(seedClWidth2[1])<3.6 &&
	  log10(seedClLength[0])>1.6 && log10(seedClLength[0])<4.0 && 
	  log10(seedClLength[1])>1.6 && log10(seedClLength[1])<4.0 ;


	//	bool singleElectronSeedCut = seedecaleonp[0]<0.85 || seedecaleonp[1]<0.85;
	//	bool singleElectronSeedCut = true;

	if (secondSeedEnergyCut) _nSel_secondseenen[isample]++;
	if (outsideEnergyCut)  _nSel_outofcone[isample]++;
	//	if (deltaPhiCut)       _nSel_deltaPhi[isample]++;
	if (acoPlanCut)        _nSel_acoPlan[isample]++;
	if (acoLinCut)         _nSel_acoLin[isample]++;
	if (leptonCut)          _nSel_lepton[isample]++;
	//	if (seedCaloEnCut)      _nSel_seedcaloen[isample]++;
	if (seedClusterCut) _nSel_seedCluster[isample]++;
	//	if (singleElectronSeedCut) _nSel_singleElec[isample]++;
	if (isrCut)            _nSel_isr[isample]++;

	bool cumulSel=true;


	if (cumulSel && secondSeedEnergyCut ) _nCumulSel_secondseenen[isample]++;
	else cumulSel=false;

	if (cumulSel && outsideEnergyCut)  _nCumulSel_outofcone[isample]++;
	else cumulSel=false;

	if (cumulSel && acoLinCut)     _nCumulSel_acoLin[isample]++;
	else cumulSel=false;

	//	if (cumulSel && deltaPhiCut)       _nCumulSel_deltaPhi[isample]++;
	if (cumulSel && acoPlanCut)       _nCumulSel_acoPlan[isample]++;
	else cumulSel=false;

	if (cumulSel && isrCut)            _nCumulSel_isr[isample]++;
	else cumulSel=false;

	//if (cumulSel && singleElectronSeedCut)  _nCumulSel_singleElec[isample]++;
	//else cumulSel=false;

	if (cumulSel && leptonCut)          _nCumulSel_lepton[isample]++;
	else cumulSel=false;

	//	if (cumulSel && seedCaloEnCut)          _nCumulSel_seedcaloen[isample]++;
	if (cumulSel && seedClusterCut)   _nCumulSel_seedCluster[isample]++;
        else cumulSel=false;

	//	bool select = secondSeedEnergyCut && acoLinCut && outsideEnergyCut && leptonCut && isrCut;

	if ( !cumulSel ) return;


 if (verboseFF) cout << "hello 7" << endl;

	//------------------------------------
	//  end of first selection
	//---------------------------------

	//	cout << "Event " << evt->getEventNumber() << "passed 1st selection (sample = " << isample << endl;

	// -----------------------------------
	// try to trim the cones somewhat
	// -----------------------------------


	TLorentzVector trimmed_jet4mom[2];
	int trimmed_jetChg[2];

	// -----------------------------------
	// if neutral hadrons are inside the cone, they are likely fragments of the charged hadron...
	// -----------------------------------

	std::vector < ReconstructedParticle* > coneparticlesTRIM[2];
        for (int j=0; j<2; j++) {

	  std::vector  < ReconstructedParticle* > tempChHad;
	  std::vector  < ReconstructedParticle* > tempNeuHad;

	  for ( size_t k=0; k<coneparticles[j].size(); k++) {
	    if ( coneparticles[j][k]->getType()!=2112 ) {
	      coneparticlesTRIM[j].push_back(coneparticles[j][k]);

	      if ( coneparticles[j][k]->getType()==211 ) tempChHad.push_back(coneparticles[j][k]);

	    } else {
	      // cout << "neutral hadron in cone" << endl;

	      tempNeuHad.push_back(coneparticles[j][k]);


	      ReconstructedParticle* pfo = coneparticles[j][k];

	      //cout << "reco->mc related MCparticles: ";
	      //for (size_t kk=0; kk<_relNavi->getRelatedToObjects( pfo ).size(); kk++) {
	      //	MCParticle* mcp = dynamic_cast <MCParticle*> (_relNavi->getRelatedToObjects( pfo )[kk]);
	      //	float wgt = _relNavi->getRelatedToWeights(pfo)[kk];
	      //	float clusterwgt = (int(wgt)/10000)/1000. ;
	      //	int parPdg = mcp->getParents().size()==0 ? 0 : mcp->getParents()[0]->getPDG();
	      //	cout << mcp->getPDG() << " " << mcp->getEnergy() << " " << parPdg << " " << clusterwgt << " , ";
	      //}
	      //cout << endl;
	      //
	      //cout << "mc->reco related MCparticles: ";
	      //for (size_t kk=0; kk<_relNavi2->getRelatedFromObjects( pfo ).size(); kk++) {
	      //	MCParticle* mcp = dynamic_cast <MCParticle*> (_relNavi2->getRelatedFromObjects( pfo )[kk]);
	      //	float wgt = _relNavi2->getRelatedFromWeights( pfo )[kk];
	      //	float clusterwgt = (int(wgt)/10000)/1000. ;
	      //	cout << mcp->getPDG() << " " << mcp->getEnergy() << " " << clusterwgt << " , ";
	      //}
	      //cout << endl;


	      MCParticle* mcpp = dynamic_cast <MCParticle*> (_relNavi->getRelatedToObjects( pfo )[0]);
	      //cout << "main contrib: "<< mcpp->getPDG() << " chg " << mcpp->getCharge() << " prodVtx " << sqrt( pow(mcpp->getVertex()[0],2)+pow(mcpp->getVertex()[1],2)+pow(mcpp->getVertex()[2],2) ) ;
	      //if (mcpp->getParents().size()>0 ) {
	      //	cout << " parent " << mcpp->getParents()[0]->getPDG() ;
	      // }
	      //cout <<endl;

	      int icode(0);
	      switch ( abs( mcpp->getPDG() ) ) {
	      case 211:
	      case 2212:        // proton
	      case 321:         // chg K
	      case 3312:        // sigma-
		icode=1; break;
	      case 22:
		icode=2; break;
	      case 11:
		icode=3; break;
	      case 13:          // mu
		icode=4; break;
	      case 2112:        // neutron
	      case 130:         // Klong
	      case 310:         // Kshrt
	      case 3122:        // lambda
		icode=5; break;
	      default:
		cout << "define code for "<< abs( mcpp->getPDG() ) << endl;
	      }

	      h_neutralHadronPFO_mainMCcontrib_pdgEn ->Fill( pfo->getEnergy(), icode );

	      //	      cout << mcpp->getPDG() << endl;

	    }
	  } // loop over cone particles

	  if ( tempChHad.size()==1 && tempNeuHad.size()==1 ) {
	    float chgMom=0;
	    for (int ll=0; ll<3; ll++) {
	      chgMom += pow( tempChHad[0]->getMomentum()[ll], 2 );
	    }
	    chgMom=sqrt(chgMom);

	    float chgClusEn = 0;
	    for (size_t ll=0; ll<tempChHad[0]->getClusters().size(); ll++) {
	      chgClusEn+=tempChHad[0]->getClusters()[ll]->getEnergy();
	    }
	    float neuClusEn=0;
	    for (size_t ll=0; ll<tempNeuHad[0]->getClusters().size(); ll++) {
	      neuClusEn+=tempNeuHad[0]->getClusters()[ll]->getEnergy();
	    }


	    //	    cout << "chg/neu " << chgMom << " " << chgClusEn << " " << neuClusEn << endl;

	    float eOnP_orig = chgClusEn/chgMom;
	    float eOnP_merge = (chgClusEn+neuClusEn)/chgMom;

	    h_neutronClus_eOnP_beforeAfter->Fill(eOnP_orig, eOnP_merge);

	  }


	  TLorentzVector totalfmomTRIM;
	  for ( size_t k=0; k<coneparticlesTRIM[j].size(); k++) {
	    // cout << coneparticles[j][k]->getType() << " " << coneparticles[j][k]->getEnergy() << endl;
	    TLorentzVector fmom =  tauUtils::getTLV( coneparticlesTRIM[j][k] ) ;
	    totalfmomTRIM+=fmom;
	  }

	  // now look at the charge
	  // cout << "charge = " << cone_chg[j] << endl;
	  //if ( isample==0 && fabs( cone_chg[j] ) != 1 ) cout << "BAD CHARGE!!" << endl;
	  int ifurthest(-1);
	  float maxangle(0);
	  if ( fabs( cone_chg[j] ) <0.5 ) {
	    // look for charged pfo farthest from cone axis
	    for ( size_t k=0; k<coneparticlesTRIM[j].size(); k++) {
	      if ( isCharged ( coneparticlesTRIM[j][k] )  ) {
		TLorentzVector fmom =  tauUtils::getTLV( coneparticlesTRIM[j][k] ) ;
		float angle = fmom.Vect().Angle( totalfmomTRIM.Vect() );
		// cout << angle << " " << coneparticles[j][k]->getCharge() << endl;
		if ( angle>maxangle ) {
		  maxangle=angle;
		  ifurthest=k;
		}
	      }
	    }
	    if ( ifurthest>=0 ) {
	      //cout << "erasing furthest charged particle : energy " << coneparticles[j][ifurthest]->getEnergy() << " angle " << maxangle << endl;
	      coneparticlesTRIM[j].erase( coneparticlesTRIM[j].begin()+ifurthest );
	    } else {
	      cout << "ERROR! cannot find furthest charged cone member!" << endl;
	      assert(0);
	    }
	  }

	  // now look at the invariant mass
	  totalfmomTRIM.SetXYZT(0,0,0,0);
	  int trimChg(0);
	  for ( size_t k=0; k<coneparticlesTRIM[j].size(); k++) {
	    TLorentzVector fmom =  tauUtils::getTLV( coneparticlesTRIM[j][k] ) ;
	    totalfmomTRIM+=fmom;
	    trimChg+=coneparticlesTRIM[j][k]->getCharge();
	  }
	  // cout << "invariant mass after first trim: " << totalfmom.M() << endl;

	  // if ( totalfmom.M()> 1.77 && isample==0 ) cout << "TOO MASSIVE" << endl;

	  trimmed_jet4mom[j]=totalfmomTRIM;
	  trimmed_jetChg[j]=trimChg;

	  h_rawmass_trimmass[isample]->Fill(jet4mom[j].M() , trimmed_jet4mom[j].M() );

	}

	// now lets see how well we did by comparing to MC
	float taumatchedEnergy(0);
	float nontaumatchedEnergy(0);
	float unmatchedEnergy(0);
	for ( int j=0; j<2; j++) { // the 2 cones
	  for ( size_t k=0; k<coneparticles[j].size(); k++) {
	    ReconstructedParticle* rp = coneparticles[j][k];
	    int mtau(-1);
	    for (int itau=0; itau<2; itau++) {
	      if ( find( PFOfromtau[itau].begin(), PFOfromtau[itau].end(), rp )!=PFOfromtau[itau].end() ) {
		mtau=itau;
		break;
	      }
	    }
	    if ( mtau>=0 ) {
	      taumatchedEnergy+=rp->getEnergy();
	    } else {
	      if ( find ( noMatchPFOs.begin(), noMatchPFOs.end(), rp ) != noMatchPFOs.end() ) {
		unmatchedEnergy+=rp->getEnergy(); // not matched to a MC particle
	      } else {
		nontaumatchedEnergy+=rp->getEnergy(); // matched to non-tau MC particle
	      }
	    }
	  }
	}
	//	cout << "MC tau-matched / non-tau / unmatched energy = " << taumatchedEnergy << " / " << nontaumatchedEnergy << " / " << unmatchedEnergy << endl;


 if (verboseFF) cout << "hello 8" << endl;

	//---------second selection on tau jets

	float maxMass = 1.77;
	// bool tauMassSel = jet4mom[0].M()<maxMass && jet4mom[1].M()<maxMass;
	bool tauMassSel = trimmed_jet4mom[0].M()<maxMass && trimmed_jet4mom[1].M()<maxMass;


	// bool chargeCut = cone_chg[0]*cone_chg[1]==-1;
	bool chargeCut = trimmed_jetChg[0]*trimmed_jetChg[1]==-1;

	//	select = select && tauMassSel && chargeCut;

	if (chargeCut)         _nSel_chg[isample]++;
	if (tauMassSel)        _nSel_tmass[isample]++;
	//	if (select)            _nSel[isample]++;

	if (cumulSel && tauMassSel)        _nCumulSel_tmass[isample]++;
	else cumulSel=false;

	if (cumulSel && chargeCut)         _nCumulSel_chg[isample]++;
	else cumulSel=false;

	//	if (cumulSel && select)            _nCumulSel[isample]++;

	//---------------------
	
	if ( tautauInvMass>0 && cumulSel ) {
	  hSEL_ttmass->Fill(tautauInvMass);
	  hSEL_mcTau_costh[isample]->Fill(taucosth[0], taucosth[1]);
	}

	for (int j=0; j<2; j++) {	    // j = cone index
	  if ( finalmctaus.size()==2 && matchedTauByDir[0] && matchedTauByDir[1] ) {
	    int mctauIndex(-1);
	    if ( matchedTauByDir[j] == finalmctaus[0] ) mctauIndex=0;
	    else if ( matchedTauByDir[j] == finalmctaus[1] ) mctauIndex=1;
	    else {
	      cout << "unmatched tau by dir! " << endl;
	      assert(0);
	    }
	    h_dec_tmass [isample]-> Fill ( MCdecayMode[mctauIndex], trimmed_jet4mom[j].M() );
	    h_dec_ngam  [isample]-> Fill ( MCdecayMode[mctauIndex], cone_ngam[j] );
	    h_dec_nchg  [isample]-> Fill ( MCdecayMode[mctauIndex], cone_nchg[j] );
	    h_dec       [isample]-> Fill ( MCdecayMode[mctauIndex] );
	  }
	  h_tchg_tmass[isample]->Fill( trimmed_jetChg[j], trimmed_jet4mom[j].M() );
	  h_ngam_tmass[isample]->Fill( cone_ngam[j], trimmed_jet4mom[j].M() );
	  h_nchg_tmass[isample]->Fill( cone_nchg[j], trimmed_jet4mom[j].M() );
	  h_nnhad_tmass[isample]->Fill( cone_nnhad[j], trimmed_jet4mom[j].M() );
	  h_ngam_nchg[isample]->Fill ( cone_ngam[j], cone_nchg[j] );
	}

	if ( cumulSel ) {

	  if ( mctauminus ) hSEL_mc_tauMinus_costh[isample]->Fill( tauUtils::getTLV(mctauminus).Vect().CosTheta() );


	  for (size_t i=0; i<finalmctaus.size(); i++) {
	    if (  MCdecayMode[i] == tauUtils::decayChPi ) {
	      if ( mctauhelicity[i]<0 ) {
		hSEL_pi_mcPolarExact_helNeg[isample]->Fill(mcexactPolarimeter[i]);
		hSEL_pi_mcPolarApprox_helNeg[isample]->Fill(mcapproxPolarimeter[i]);
	      } else {
		hSEL_pi_mcPolarExact_helPos[isample]->Fill(mcexactPolarimeter[i]);
		hSEL_pi_mcPolarApprox_helPos[isample]->Fill(mcapproxPolarimeter[i]);
	      }
	    } else if (  MCdecayMode[i] == tauUtils::decayRho ) {
	      if ( mctauhelicity[i]<0 ) {
		hSEL_rho_mcPolarExact_helNeg[isample]->Fill(mcexactPolarimeter[i]);
		hSEL_rho_mcPolarApprox_helNeg[isample]->Fill(mcapproxPolarimeter[i]);
	      } else {
		hSEL_rho_mcPolarExact_helPos[isample]->Fill(mcexactPolarimeter[i]);
		hSEL_rho_mcPolarApprox_helPos[isample]->Fill(mcapproxPolarimeter[i]);
	      }
	    }
	  }







	  hSEL_visjetAcolinearity[isample]   ->Fill( visjetacolinearity );
	  hSEL_visjetAcoplanarity [isample]->Fill( visjetacoplanarity );

	  hSEL_prongAngle[isample]   ->Fill( prongangle );
	  hSEL_outsideEnergy[isample]->Fill( outsideConeEnergy );
	  hSEL_insideEnergy[isample] ->Fill( trimmed_jet4mom[0].E() + trimmed_jet4mom[1].E() );
	  hSEL_outsidePt[isample]    ->Fill( outsideConePt );
	  for (int j=0; j<2; j++) {
	    hSEL_dec[isample]    -> Fill ( MCdecayMode[j] );
	    hSEL_tchg_tmass[isample]->Fill( trimmed_jetChg[j], trimmed_jet4mom[j].M() );
	    hSEL_ngam_tmass[isample]->Fill( cone_ngam[j], trimmed_jet4mom[j].M() );
	    hSEL_nchg_tmass[isample]->Fill( cone_nchg[j], trimmed_jet4mom[j].M() );
	    hSEL_nnhad_tmass[isample]->Fill( cone_nnhad[j], trimmed_jet4mom[j].M() );
	    hSEL_ngam_nchg[isample]->Fill(  cone_ngam[j], cone_nchg[j] );
	  }
	  hSEL_maxSeedEn[isample]->Fill( highestPtChargedPFO[0].first );
	  hSEL_minSeedEn[isample]->Fill( highestPtChargedPFO[1].first );

	  hSEL_maxSeedEn_eonp[isample]->Fill( seedecaleonp[0] );
	  hSEL_minSeedEn_eonp[isample]->Fill( seedecaleonp[1] );

	  hSEL_maxSeedEn_caloen[isample]->Fill( seedcaloen[0] );
	  hSEL_minSeedEn_caloen[isample]->Fill( seedcaloen[1] );

	  hSEL_maxSeedCosth[isample]->Fill( costhMax );
	  hSEL_minSeedCosth[isample]->Fill( costhMin );

	  hSEL_seedEcalEn[isample]->Fill( seedecalen[0] );
	  hSEL_seedEcalEn[isample]->Fill( seedecalen[1] );
	  hSEL_seedEcalEonP[isample]->Fill( seedecaleonp[0] );
	  hSEL_seedEcalEonP[isample]->Fill( seedecaleonp[1] );
	  hSEL_ooconeMaxGammaEn[isample]->Fill( maxisr );
	  hSEL_ooconeMaxPFOEn[isample]->Fill( maxoutsidepfoen );

	  for (int j=0; j<2; j++) {
	    hSEL_coneMass     [isample]->Fill( trimmed_jet4mom[j].M() );
	  }
	  hSEL_seedSeedAngle[isample]->Fill( TMath::Pi() - coneSeedDir[0].Angle(coneSeedDir[1])  );
	  hSEL_seedSeedDphi [isample]->Fill( TMath::Pi() - fabs(coneSeedDir[0].DeltaPhi(coneSeedDir[1]))  );

	  int iminus(0);
	  if ( trimmed_jetChg[0]>0 ) {
	    iminus=1;
	  }
	  hSEL_tauMinusCosth [isample]->Fill( trimmed_jet4mom[iminus].Vect().CosTheta() );

	  int recoDecay[2]={-1,-1};
	  TLorentzVector trimneutraljet4mom[2];

	  for (int j=0; j<2; j++) {

	    trimneutraljet4mom[j].SetXYZT(0,0,0,0);

	    int ngam(0);
	    float engam(0);
	    int nhad(0);
	    float enhad(0);
	    int nch(0);
	    float ench(0);
	    int ncompound(0);
	    for ( size_t i=0; i<coneparticlesTRIM[j].size(); i++) {
	      if ( isCharged ( coneparticlesTRIM[j][i] ) ) {
		nch++;
		ench+=coneparticlesTRIM[j][i]->getEnergy();
	      } else {
		if ( coneparticlesTRIM[j][i]->getType()==22 ) {
		  ngam++;
		  engam+=coneparticlesTRIM[j][i]->getEnergy();

		  trimneutraljet4mom[j] += TLorentzVector( coneparticlesTRIM[j][i]->getMomentum()[0], 
						       coneparticlesTRIM[j][i]->getMomentum()[1], 
						       coneparticlesTRIM[j][i]->getMomentum()[2], 
						       coneparticlesTRIM[j][i]->getEnergy() );

		} else if ( coneparticlesTRIM[j][i]->getType()==111 ) {
		  ngam+=2;
		  engam+=coneparticlesTRIM[j][i]->getEnergy();

		  trimneutraljet4mom[j] += TLorentzVector( coneparticlesTRIM[j][i]->getMomentum()[0], 
						       coneparticlesTRIM[j][i]->getMomentum()[1], 
						       coneparticlesTRIM[j][i]->getMomentum()[2], 
						       coneparticlesTRIM[j][i]->getEnergy() );

		} else if  ( coneparticlesTRIM[j][i]->getType()==2112 ) {
		  nhad++;
		  enhad+=coneparticlesTRIM[j][i]->getEnergy();
		} else if (  coneparticlesTRIM[j][i]->getType()==221 || coneparticlesTRIM[j][i]->getType()==310 || abs(coneparticlesTRIM[j][i]->getType())==3122 ||  coneparticlesTRIM[j][i]->getType()==331 ) {
		  ncompound++;
		} else {
		  cout << "unknown PFO type !!!" <<  coneparticlesTRIM[j][i]->getType() << endl;
		  assert(0);
		}
	      }
	    }
	    

	    int decay= matchedTauByDirDecay[j];


	    // look at merged photon events (?)
	    if ( decay == tauUtils::decayRho ) {

	      if (ngam==2) { // expected # of photon-like PFOs

		for ( auto cp : coneparticlesTRIM[j] ) {
		  if ( cp->getType()==22 && cp->getClusters().size()==1 ) {
		    std::vector <float> evals = tauUtils::getClusterEigenvalues( cp->getClusters()[0] );
		    h_rhoDecaySinglePhoClus_singleGammaClusterEvals->Fill(evals[0], evals[1]);
		  }
		}
		
	      } else if ( ngam==1 ) {
		  

		ReconstructedParticle* rpgam(0);
		  
		for ( auto cp : coneparticlesTRIM[j] ) {
		  if ( cp->getType()==22 ) {
		    rpgam = cp;
		    break;
		  }
		}

		int nclus =  rpgam->getClusters().size();
		std::vector <float> evals;
		if ( nclus==1 ) {
		  Cluster* cl = rpgam->getClusters()[0];
		  evals = tauUtils::getClusterEigenvalues(cl);
		  // cout << "single gam cluster: eigenvalues " << evals[0] << " " << evals[1] << " " << evals[2] << endl;
		}


		std::vector <MCParticle*> asd = tauUtils::getstablemctauDaughters( matchedTauByDir[j] );
		std::vector <MCParticle*> mcgamma;
		for ( auto t : asd ) {
		  if ( t->getPDG()==22 ) mcgamma.push_back(t);
		}
		//		cout << "single phoClus rho decay: nmc photons = " << mcgamma.size() << endl;

		std::map<MCParticle*, ReconstructedParticle*> majorPFO;

		for ( auto mcj : mcgamma ) {
		  //		  cout << "mc photon energy:" << mcj->getEnergy() << " decayed in tracker? " << mcj->isDecayedInTracker() << endl;

		  // //		  cout << "mc->reco related PFOs: ";
		  // for (size_t kk=0; kk<_relNavi->getRelatedFromObjects( mcj ).size(); kk++) {
		  //   ReconstructedParticle* pfo = dynamic_cast <ReconstructedParticle*> (_relNavi->getRelatedFromObjects( mcj )[kk]);
		  //   float wgt = _relNavi->getRelatedFromWeights( mcj )[kk];
		  //   float clusterwgt = (int(wgt)/10000)/1000. ;
		  //   //		    cout << "en " << pfo->getEnergy() << " type " << pfo->getType() << " clwgt " << clusterwgt << " ; ";
		  // }
		  // //		  cout << endl;
		
		  float biggestclwt(0);
		  
		  //		  cout << "reco->mc related PFOs: ";
		  for (size_t kk=0; kk<_relNavi2->getRelatedToObjects( mcj ).size(); kk++) {
		    ReconstructedParticle* pfo = dynamic_cast <ReconstructedParticle*> (_relNavi2->getRelatedToObjects( mcj )[kk]);
		    float wgt = _relNavi2->getRelatedToWeights( mcj )[kk];
		    float clusterwgt = (int(wgt)/10000)/1000. ;
		    //		    cout << "en " << pfo->getEnergy() << " type " << pfo->getType() << " clwgt " << clusterwgt << " ; ";
		  
		    if ( clusterwgt>biggestclwt ) {
		      biggestclwt=clusterwgt;
		      majorPFO[mcj]=pfo;
		    }

		  }
		  //		  cout << endl;
		
		  //		cout << "mc->reco related MCparticles: ";
		  //		for (size_t kk=0; kk<_relNavi2->getRelatedFromObjects( pfo ).size(); kk++) {
		  //			MCParticle* mcp = dynamic_cast <MCParticle*> (_relNavi2->getRelatedFromObjects( pfo )[kk]);
		  //			float wgt = _relNavi2->getRelatedFromWeights( pfo )[kk];
		  //			float clusterwgt = (int(wgt)/10000)/1000. ;
		  //			cout << mcp->getPDG() << " " << mcp->getEnergy() << " " << clusterwgt << " , ";
		  //		}
		  //		cout << endl;
		  

		}

		////		cout << mcgamma.size() << " " << majorPFO.size() << endl;
		//for ( auto mcpf : majorPFO ) {
		//  cout << " mcgammen " << mcpf.first->getEnergy() << " pfoen " << mcpf.second->getEnergy() << " ; ";
		//}
		//cout << endl;

		if (  mcgamma.size()==2 ) {

		  float angle = TVector3( mcgamma[0]->getMomentum() ).Angle( TVector3( mcgamma[1]->getMomentum() ) );
		  //		  cout << "angle " << angle;


		  bool oneIsDecayed = mcgamma[0]->isDecayedInTracker() || mcgamma[1]->isDecayedInTracker();

		  if ( oneIsDecayed ) {
		    h_rhoDecaySinglePhoClus_reason->Fill( 0 );
		  } else {
		    if ( majorPFO.find( mcgamma[0] ) == majorPFO.end() || majorPFO.find( mcgamma[1] ) == majorPFO.end() ) {
		      
		      if (  mcgamma[0]->getEnergy()<0.3 || mcgamma[1]->getEnergy()<0.3 ) {
			h_rhoDecaySinglePhoClus_reason->Fill( 1 );
		      } else {
			h_rhoDecaySinglePhoClus_reason->Fill( 2 );
		      }
		      
		    } else {
		      
		      bool merged = majorPFO[mcgamma[0]]==majorPFO[mcgamma[1]];
		      //		      cout << " merged into pho clus ? " << merged;

		    


		      if ( merged ) {
			if ( majorPFO[mcgamma[0]]->getType()==22 ) {
			  h_rhoDecaySinglePhoClus_reason->Fill( 3 );
			  
			  if ( evals.size()>0 ) h_rhoDecaySinglePhoClus_mergedGammaClusterEvals->Fill(evals[0], evals[1]);


			  
			} else if (  majorPFO[mcgamma[0]]->getCharge()==0 ) {
			  h_rhoDecaySinglePhoClus_reason->Fill( 4 );
			} else {
			  h_rhoDecaySinglePhoClus_reason->Fill( 5 );
			}
		      } else {
			if ( majorPFO[mcgamma[0]]->getType()==211 || majorPFO[mcgamma[1]]->getType()==211 ) {
			  h_rhoDecaySinglePhoClus_reason->Fill( 6 );
			} else {
			  //			  cout << "hihihi" << endl;
			  h_rhoDecaySinglePhoClus_reason->Fill( 9 );
			}
		      }
		    }
		  }
		}
	      
		//		cout << endl;
	      }	      
	    }
	      


	    for (int jj=0; jj<2; jj++) {

	      int idec(-1);

	      if (jj==0) {
		idec=NDEC-1; // all events
	      } else if ( decay>=0 ) {
		idec= std::min( decay, NDEC-2 );
		//idec= decay<6 ? decay : NDEC-2;
	      }

	      if ( idec>=0 ) {

		hSEL_dec_seed_clusterWidth1[isample][idec]->Fill( log10(seedClWidth1[j]) );
		hSEL_dec_seed_clusterWidth2[isample][idec]->Fill( log10(seedClWidth2[j]) );
		hSEL_dec_seed_clusterLength[isample][idec]->Fill( log10(seedClLength[j]) );
		hSEL_dec_seed_energy[isample][idec]->Fill( highestPtChargedPFO[j].second->getEnergy() );

		hSEL_dec_cone_nchpfo[isample][idec]->Fill(nch);
		hSEL_dec_cone_ntracks[isample][idec]->Fill(nTrkInCone[j]);

		for ( auto tt : tracksInCone[j] ) {
		  hSEL_dec_cone_trackNTPC[isample][idec]->Fill(tt->getSubdetectorHitNumbers()[6]);
		}


		hSEL_dec_cone_nnhadpfo[isample][idec]->Fill(nhad);
		hSEL_dec_cone_ngammapfo[isample][idec]->Fill(ngam);
		hSEL_dec_cone_ncompoundpfo[isample][idec]->Fill(ncompound);
		hSEL_dec_cone_nhaden[isample][idec]->Fill(enhad);
		hSEL_dec_cone_gammaen[isample][idec]->Fill(engam);
		hSEL_dec_cone_nhadenFrac[isample][idec]->Fill(enhad/(enhad+engam+ench));
		hSEL_dec_cone_gammaenFrac[isample][idec]->Fill(engam/(enhad+engam+ench));

		hSEL_dec_cone_visMass[isample][idec]->Fill( jet4mom[j].M() );
		hSEL_dec_cone_neutralvisMass[isample][idec]->Fill( jetNeutral4mom[j].M() );
		hSEL_dec_cone_visMassDiff[isample][idec]->Fill( jet4mom[j].M() - mcTauVis4mom[matchedTauByDirIndex[j]].M() );
		hSEL_dec_cone_neutralvisMassDiff[isample][idec]->Fill( jetNeutral4mom[j].M()  - mcTauNeutralVis4mom[matchedTauByDirIndex[j]].M() );

		hSEL_dec_cone_vis_neutral_Mass[isample][idec]->Fill( jet4mom[j].M(),  jetNeutral4mom[j].M() );
		if ( ngam==0) hSEL_dec_cone_vis_neutral_Mass_0g[isample][idec]->Fill( jet4mom[j].M(),  jetNeutral4mom[j].M() );
		else if ( ngam==1) hSEL_dec_cone_vis_neutral_Mass_1g[isample][idec]->Fill( jet4mom[j].M(),  jetNeutral4mom[j].M() );
		else if ( ngam==2) hSEL_dec_cone_vis_neutral_Mass_2g[isample][idec]->Fill( jet4mom[j].M(),  jetNeutral4mom[j].M() );
		else if ( ngam==3) hSEL_dec_cone_vis_neutral_Mass_3g[isample][idec]->Fill( jet4mom[j].M(),  jetNeutral4mom[j].M() );
		else if ( ngam>=4) hSEL_dec_cone_vis_neutral_Mass_4g[isample][idec]->Fill( jet4mom[j].M(),  jetNeutral4mom[j].M() );

		hSEL_dec_cone_visEnergyDiff[isample][idec]->Fill( jet4mom[j].E() - mcTauVis4mom[matchedTauByDirIndex[j]].E() );
		hSEL_dec_cone_neutralvisEnergyDiff[isample][idec]->Fill( jetNeutral4mom[j].E()  - mcTauNeutralVis4mom[matchedTauByDirIndex[j]].E() );
		
		hSEL_dec_coneTRIM_visMass[isample][idec]->Fill( trimmed_jet4mom[j].M() );
		hSEL_dec_coneTRIM_neutralvisMass[isample][idec]->Fill( trimneutraljet4mom[j].M() );

		hSEL_dec_coneTRIM_vis_neutral_Mass[isample][idec]->Fill( trimmed_jet4mom[j].M(), trimneutraljet4mom[j].M() );
		if ( ngam==0) hSEL_dec_coneTRIM_vis_neutral_Mass_0g[isample][idec]->Fill( trimmed_jet4mom[j].M(), trimneutraljet4mom[j].M() );
		else if ( ngam==1) hSEL_dec_coneTRIM_vis_neutral_Mass_1g[isample][idec]->Fill( trimmed_jet4mom[j].M(), trimneutraljet4mom[j].M() );
		else if ( ngam==2) hSEL_dec_coneTRIM_vis_neutral_Mass_2g[isample][idec]->Fill( trimmed_jet4mom[j].M(), trimneutraljet4mom[j].M() );
		else if ( ngam==3) hSEL_dec_coneTRIM_vis_neutral_Mass_3g[isample][idec]->Fill( trimmed_jet4mom[j].M(), trimneutraljet4mom[j].M() );
		else if ( ngam>=4) hSEL_dec_coneTRIM_vis_neutral_Mass_4g[isample][idec]->Fill( trimmed_jet4mom[j].M(), trimneutraljet4mom[j].M() );

		hSEL_dec_coneTRIM_visMassDiff[isample][idec]->Fill( trimmed_jet4mom[j].M() - mcTauVis4mom[matchedTauByDirIndex[j]].M() );
		hSEL_dec_coneTRIM_neutralvisMassDiff[isample][idec]->Fill( trimneutraljet4mom[j].M()  - mcTauNeutralVis4mom[matchedTauByDirIndex[j]].M() );

		hSEL_dec_coneTRIM_visEnergyDiff[isample][idec]->Fill( trimmed_jet4mom[j].E() - mcTauVis4mom[matchedTauByDirIndex[j]].E() );
		hSEL_dec_coneTRIM_neutralvisEnergyDiff[isample][idec]->Fill( trimneutraljet4mom[j].E()  - mcTauNeutralVis4mom[matchedTauByDirIndex[j]].E() );
	      }

	    }

	    hSEL_mcdec[isample]->Fill(decay);


	    float jmass =  trimmed_jet4mom[j].M() ;
	    float jmassNeutral = trimneutraljet4mom[j].M();

	    if ( nch==1 ) {
	      if ( ngam==0 ) {
		if ( jmass > 0.11 && jmass < 0.2 )      recoDecay[j]=tauUtils::decayChPi;
	      } else if ( ngam==1 ) {
		if      ( engam < 5   && jmass < 0.4 )  recoDecay[j]=tauUtils::decayChPi;
		else if ( jmass > 0.4 && jmass < 0.9 )  recoDecay[j]=tauUtils::decayRho;
		else if ( jmass > 0.9 )                 recoDecay[j]=tauUtils::decayA1_1p;
	      } else if ( ngam==2 ) {
		if      ( jmass > 0.4 && jmass < 1.1 && jmassNeutral < 0.4 ) recoDecay[j]=tauUtils::decayRho;
		else if ( jmass > 0.9 && jmassNeutral>0.8 )                  recoDecay[j]=tauUtils::decayA1_1p;
	      } else if ( ngam==3 ) {
		if      ( jmass < 0.95 && jmassNeutral<0.4 )  recoDecay[j]=tauUtils::decayRho;
		else                                          recoDecay[j]=tauUtils::decayA1_1p;
	      } else if ( ngam>=4 ) {
		recoDecay[j]=tauUtils::decayA1_1p;
	      } else {
		recoDecay[j]=tauUtils::decayOthersingleprong;
	      } 
	    } else {
	      recoDecay[j]=tauUtils::decayMultiprong;
	    }

	    // polarimeter from just the energies
	    if ( recoDecay[j]==tauUtils::decayChPi ) {
	      //	      float polEstimator = ench / 250. ;
	      float polEstimator = tauUtils::getLongPolarimeterFromEn_pi(  tauUtils::getTLV( highestPtChargedPFO[j].second ), 250. );
	      hSEL_rec_pi_mcdec[isample]->Fill( matchedTauByDirDecay[j] );

	      hSEL_rec_pi_pol[isample]->Fill( polEstimator );

	      if (  mctauhelicity[j] > 0.5 ) {
		hSEL_rec_pi_pol_MCpos[isample]->Fill( polEstimator );
	      } else if (  mctauhelicity[j] < -0.5 ) {
		hSEL_rec_pi_pol_MCneg[isample]->Fill( polEstimator );
	      } else {
		hSEL_rec_pi_pol_MCoth[isample]->Fill( polEstimator );
	      }


	      hSEL_rec_pi_coneMass[isample]->Fill( jmass );

	      float dpol = mcapproxPolarimeter[j] > -98 ? polEstimator - mcapproxPolarimeter[j] : -1.5 ;


	      hSEL_recpi_mcall_pol[isample]->Fill( mcapproxPolarimeter[j] , polEstimator );
	      hSEL_recpi_mcall_dpol[isample]->Fill( dpol );
	      if ( matchedTauByDirDecay[j] == recoDecay[j] ) {
		hSEL_recpi_mcpi_pol[isample]->Fill( mcapproxPolarimeter[j] , polEstimator );
		hSEL_recpi_mcpi_dpol[isample]->Fill( dpol );
	      }

	    } else if  ( recoDecay[j]==tauUtils::decayRho ) {
	      //	      float polEstimator = ( ench - engam ) / 250. ;
	      float polEstimator = tauUtils::getLongPolarimeterFromEn_rho( tauUtils::getTLV( highestPtChargedPFO[j].second ),  trimneutraljet4mom[j] , 250. );
	      hSEL_rec_rho_mcdec[isample]->Fill( matchedTauByDirDecay[j] );
	      hSEL_rec_rho_pol[isample]->Fill( polEstimator );

	      if (  mctauhelicity[j] > 0.5 ) {
		hSEL_rec_rho_pol_MCpos[isample]->Fill( polEstimator );
	      } else if (  mctauhelicity[j] < -0.5 ) {
		hSEL_rec_rho_pol_MCneg[isample]->Fill( polEstimator );
	      } else {
		hSEL_rec_rho_pol_MCoth[isample]->Fill( polEstimator );
	      }

	      hSEL_rec_rho_coneMass[isample]->Fill( jmass  );

	      float dpol = mcapproxPolarimeter[j] > -98 ? polEstimator - mcapproxPolarimeter[j] : -1.5 ;

	      hSEL_recrho_mcall_pol[isample]->Fill( mcapproxPolarimeter[j] , polEstimator );
	      hSEL_recrho_mcall_dpol[isample]->Fill( dpol );
	      if ( matchedTauByDirDecay[j] == recoDecay[j] ) {
		hSEL_recrho_mcrho_pol[isample]->Fill( mcapproxPolarimeter[j] , polEstimator );
		hSEL_recrho_mcrho_dpol[isample]->Fill( dpol );
	      }


	    } else if (recoDecay[j]==tauUtils::decayA1_1p ) {
	      hSEL_rec_a1p_mcdec[isample]->Fill( matchedTauByDirDecay[j] );
	      hSEL_rec_a1p_coneMass[isample]->Fill( jmass  );
	    }


	    





	  } // two cones loop


	  // try to get the polarimeters for events in which both taus reconstructed to decay to pi or rho
	  //	  if ( ( recoDecay[0]==tauUtils::decayChPi || recoDecay[0]==tauUtils::decayRho ) && ( recoDecay[1]==tauUtils::decayChPi || recoDecay[1]==tauUtils::decayRho ) ) {

	  //// try to get the polarimeters for events in which a tau reconstructed to decay to pi or rho
	  //if ( ( recoDecay[0]==tauUtils::decayChPi || recoDecay[0]==tauUtils::decayRho || recoDecay[1]==tauUtils::decayChPi || recoDecay[1]==tauUtils::decayRho ) &&
	  //     ( recoDecay[0]!=tauUtils::decayMultiprong &&  recoDecay[1]!=tauUtils::decayMultiprong) ) {
	  //
	  //  int evtType;
	  //  if      ( recoDecay[0]==tauUtils::decayChPi && recoDecay[1]==tauUtils::decayChPi ) evtType=0; // both pi
	  //  else if ( recoDecay[0]==tauUtils::decayRho  && recoDecay[1]==tauUtils::decayRho  ) evtType=2; // both rho
	  //  else evtType=1; // mixed
	  //
	  //
	  //  //---------------------------------------------------------------
	  //  // back-to-back fit (cone intersections)
	  //  //---------------------------------------------------------------
	  //  
	  //  int coneFitFlag;
	  //  std::vector < std::pair <TVector3, TVector3> > nuMom = getNuMom( trimmed_jet4mom[0], trimmed_jet4mom[1], 500. , coneFitFlag);
	  //
	  //  cout << "cone fit output flag: " << coneFitFlag << " number of cone solutions : " << nuMom.size() << endl;
	  //
	  //  _hFit_coneFit_fitFlag[isample][evtType]->Fill( coneFitFlag );
	  //
	  //  TVector3 polarimeterVec[2];
	  //  float polLong[2] = {-2,-2};
	  //  for ( size_t k=0; k<nuMom.size(); k++) {
	  //    for (int i=0; i<2; i++) {
	  //	TVector3 numom= i==0? nuMom[k].first : nuMom[k].second;
	  //	TLorentzVector chargedtlv =  tauUtils::getTLV(highestPtChargedPFO[i].second);
	  //	TLorentzVector neutrinotlv; neutrinotlv.SetVectM( numom, 0 );
	  //
	  //
	  //	// TLorentzVector totalTauMom = trimmed_jet4mom[i] + neutrinotlv;
	  //	// 
	  //	// TVector3 dvec = tauUtils::getDvector( highestPtChargedPFO[i].second->getTracks()[0]->getTrackStates()[TrackState::AtIP] );
	  //	// 
	  //	// cout << "sol" << k << " tau" << i << " : " << totalTauMom.Angle( dvec ) << " " << totalTauMom.M() << endl;
	  //
	  //
	  //	if ( recoDecay[i]==tauUtils::decayChPi ) {
	  //	  polarimeterVec[i] = tauUtils::getPolarimeter_pi( chargedtlv, neutrinotlv );
	  //	  polLong[i] = cos( polarimeterVec[i].Angle( ( chargedtlv+neutrinotlv).Vect() ) );
	  //	} else if (  recoDecay[i]==tauUtils::decayRho ) {
	  //	  TLorentzVector pi0tlv = trimmed_jet4mom[i] - chargedtlv;
	  //	  polarimeterVec[i] = tauUtils::getPolarimeter_rho( chargedtlv, pi0tlv, neutrinotlv );
	  //	  polLong[i] = cos( polarimeterVec[i].Angle( ( trimmed_jet4mom[i]+neutrinotlv).Vect() ) );
	  //	} else {
	  //	  cout << "ERROR unexpected decay mode!!!" << endl;
	  //	}
	  //	if ( matchedTauByDirIndex[i]>=0 ) {
	  //	  int idecc= recoDecay[i]==tauUtils::decayChPi ? 0 : 1 ;
	  //	  _hFit_coneFit_polRecMc[isample][idecc]->Fill( mcexactPolarimeter[ matchedTauByDirIndex[i] ], polLong[i] );
	  //	  _hFit_coneFit_polRecMcDiff[isample][idecc]->Fill( polLong[i] - mcexactPolarimeter[ matchedTauByDirIndex[i] ] );
	  //	  _hFit_coneFit_polRec[isample][idecc]->Fill( polLong[i] );
	  //	}
	  //    }
	  //    _hFit_coneFit_pol[isample][evtType]->Fill( polLong[0] );
	  //    _hFit_coneFit_pol[isample][evtType]->Fill( polLong[1] );
	  //    _hFit_coneFit_pol12[isample][evtType]->Fill( polLong[0], polLong[1] );
	  //  }


	    // give up on fit for now
	    // //---------------------------------------------------------------
	    // // now try the impact param based fit
	    // //---------------------------------------------------------------
	    // 
	    // _2tauEvtFitter->reset();
	    // _2tauEvtFitter->setEcom(500.);
	    // 
	    // float avez = ( highestPtChargedPFO[0].second->getTracks()[0]->getZ0() + highestPtChargedPFO[1].second->getTracks()[0]->getZ0() )/2.;
	    // 
	    // TVector3 assumedIP (0,0, avez ) ;
	    // 
	    // //cout << "WARNING cheating PV pos!!!!" << endl;
	    // //    assumedIP = pvposMC;
	    // 
	    // _2tauEvtFitter->setIP( assumedIP );
	    // 
	    // for (int i=0; i<2; i++) {
	    //   TLorentzVector pchg( highestPtChargedPFO[i].second->getMomentum()[0],
	    // 			   highestPtChargedPFO[i].second->getMomentum()[1],
	    // 			   highestPtChargedPFO[i].second->getMomentum()[2],
	    // 			   highestPtChargedPFO[i].second->getEnergy() );
	    //   _2tauEvtFitter->setChargedHadronMomentum(i, pchg);
	    //   _2tauEvtFitter->setNeutralHadronMomentum(i, trimneutraljet4mom[i] );
	    //   TVector3 trkmom, dVec, eVec;
	    //   tauUtils::calculateTrackWRTip( highestPtChargedPFO[i].second->getTracks()[0]->getTrackState( TrackState::AtIP ), trkmom, dVec, eVec , assumedIP );
	    //   _2tauEvtFitter->setChargedHadronImpactParameterVector(i, dVec);
	    // }
	    // int strategy=0;
	    // _2tauEvtFitter->fitIt_single_single(strategy);
	    // cout << "good IP fit ? " << _2tauEvtFitter->goodFit() << endl;
	    // int nReasonableIPsolutions(0);
	    // 
	    // int besti1(-1);
	    // int besti2(-1);
	    // float bestpt(99999999);
	    // 
	    // for (int i1=0; i1<2; i1++) {
	    //   for (int i2=0; i2<2; i2++) {
	    // 	cout << "solution " << i1 << " " << i2 << endl;
	    // 	cout << "totE, pt, pz= " << _2tauEvtFitter->getEvent4Momentum(i1,i2).E() << " " << _2tauEvtFitter->getPt(i1,i2) << " " << _2tauEvtFitter->getEvent4Momentum(i1,i2).Z() << endl;
	    // 	cout << "decay lengths " << _2tauEvtFitter->getDecayLength(0, i1) << " " << _2tauEvtFitter->getDecayLength(1, i2 ) << endl;
	    // 
	    // 	_hFit_ipFit_allSols_totE[isample][evtType]->Fill( _2tauEvtFitter->getEvent4Momentum(i1,i2).E() );
	    // 	_hFit_ipFit_allSols_totPt[isample][evtType]->Fill( _2tauEvtFitter->getPt(i1,i2) );
	    // 	_hFit_ipFit_allSols_totPz[isample][evtType]->Fill( _2tauEvtFitter->getEvent4Momentum(i1,i2).Pz() );
	    // 	_hFit_ipFit_allSols_decayLength[isample][evtType]->Fill( _2tauEvtFitter->getDecayLength(0, i1) );
	    // 	_hFit_ipFit_allSols_decayLength[isample][evtType]->Fill( _2tauEvtFitter->getDecayLength(1, i2) );
	    // 
	    // 	if ( _2tauEvtFitter->getDecayLength(0, i1) > 0 && _2tauEvtFitter->getDecayLength(0, i1) < 100 &&
	    // 	     _2tauEvtFitter->getDecayLength(1, i2) > 0 && _2tauEvtFitter->getDecayLength(1, i2) < 100 &&
	    // 	     _2tauEvtFitter->getEvent4Momentum(i1,i2).E() > 450. &&  _2tauEvtFitter->getEvent4Momentum(i1,i2).E() < 550. &&
	    // 	     _2tauEvtFitter->getPt(i1,i2) < 25. && 
	    // 	     fabs(_2tauEvtFitter->getEvent4Momentum(i1,i2).Pz()) < 25. ) {
	    // 
	    // 	  if ( _2tauEvtFitter->getPt(i1,i2) < bestpt ) {
	    // 	    bestpt= _2tauEvtFitter->getPt(i1,i2);
	    // 	    besti1=i1;
	    // 	    besti2=i2;
	    // 	  }
	    // 
	    // 	  _hFit_ipFit_reasonableSols_totE[isample][evtType]->Fill( _2tauEvtFitter->getEvent4Momentum(i1,i2).E() );
	    // 	  _hFit_ipFit_reasonableSols_totPt[isample][evtType]->Fill( _2tauEvtFitter->getPt(i1,i2) );
	    // 	  _hFit_ipFit_reasonableSols_totPz[isample][evtType]->Fill( _2tauEvtFitter->getEvent4Momentum(i1,i2).Pz() );
	    // 	  _hFit_ipFit_reasonableSols_decayLength[isample][evtType]->Fill( _2tauEvtFitter->getDecayLength(0, i1) );
	    // 	  _hFit_ipFit_reasonableSols_decayLength[isample][evtType]->Fill( _2tauEvtFitter->getDecayLength(1, i2) );
	    // 
	    // 	  TVector3 polarimeterVecIP[2];
	    // 	  float polLongIP[2];
	    // 	  for (int i=0; i<2; i++) {
	    // 	    TLorentzVector chargedtlv =  tauUtils::getTLV(highestPtChargedPFO[i].second);
	    // 	    TLorentzVector neutrinotlv = i==0 ? 
	    // 	      _2tauEvtFitter->getNeutrino4Momentum( i, i1 ) :
	    // 	      _2tauEvtFitter->getNeutrino4Momentum( i, i2 );
	    // 	    if ( recoDecay[i]==tauUtils::decayChPi ) {
	    // 	      polarimeterVecIP[i] = tauUtils::getPolarimeter_pi( chargedtlv, neutrinotlv );
	    // 	      polLongIP[i] = cos( polarimeterVecIP[i].Angle( ( chargedtlv+neutrinotlv).Vect() ) );
	    // 	    } else if (  recoDecay[i]==tauUtils::decayRho ) {
	    // 	      TLorentzVector pi0tlv = trimmed_jet4mom[i] - chargedtlv;
	    // 	      polarimeterVecIP[i] = tauUtils::getPolarimeter_rho( chargedtlv, pi0tlv, neutrinotlv );
	    // 	      polLongIP[i] = cos( polarimeterVecIP[i].Angle( ( trimmed_jet4mom[i]+neutrinotlv).Vect() ) );
	    // 	    } else {
	    // 	      cout << "ERROR unexpected decay mode!!!" << endl;
	    // 	    }
	    // 	    if ( matchedTauByDirIndex[i]>=0 ) {
	    // 	      int idecc= recoDecay[i]==tauUtils::decayChPi ? 0 : 1 ;
	    // 	      _hFit_ipFit_polRecMc[isample][idecc]->Fill( mcexactPolarimeter[ matchedTauByDirIndex[i] ], polLongIP[i] );
	    // 	      _hFit_ipFit_polRecMcDiff[isample][idecc]->Fill( polLongIP[i] - mcexactPolarimeter[ matchedTauByDirIndex[i] ] );
	    // 	      _hFit_ipFit_polRec[isample][idecc]->Fill( polLongIP[i] );
	    // 	    }
	    // 	  }
	    // 	  _hFit_ipFit_pol12[isample][evtType]->Fill( polLongIP[0], polLongIP[1] );
	    // 	  _hFit_ipFit_pol[isample][evtType]->Fill( polLongIP[0] );
	    // 	  _hFit_ipFit_pol[isample][evtType]->Fill( polLongIP[1] );
	    // 
	    // 	  nReasonableIPsolutions++;
	    // 	}
	    //   }
	    // }
	    // cout << "nReasonableIPsolutions " << nReasonableIPsolutions << endl;
	    // 
	    // if ( besti1>=0 && besti2>=0 ) {
	    //   _hFit_ipFit_bestSol_totE[isample][evtType]->Fill( _2tauEvtFitter->getEvent4Momentum(besti1,besti2).E() );
	    //   _hFit_ipFit_bestSol_totPt[isample][evtType]->Fill( _2tauEvtFitter->getPt(besti1,besti2) );
	    //   _hFit_ipFit_bestSol_totPz[isample][evtType]->Fill( _2tauEvtFitter->getEvent4Momentum(besti1,besti2).Pz() );
	    //   _hFit_ipFit_bestSol_decayLength[isample][evtType]->Fill( _2tauEvtFitter->getDecayLength(0, besti1) );
	    //   _hFit_ipFit_bestSol_decayLength[isample][evtType]->Fill( _2tauEvtFitter->getDecayLength(1, besti2) );
	    // }
	    // 
	    // 
	    // 
	    // _hFit_nSolutions_cone_ip[isample][evtType]->Fill(  nuMom.size(), nReasonableIPsolutions );
	    // _hFit_nSolutions_cone[isample][evtType]->Fill(  nuMom.size() );
	    // _hFit_nSolutions_ip[isample][evtType]->Fill( nReasonableIPsolutions );



	  //	  }


	} // select


 if (verboseFF) cout << "hello 10" << endl;

	if ( cumulSel ) {

	  // write out the taus
	  LCCollectionVec* newcol = new LCCollectionVec();
	  for ( int j=0; j<2; j++) {

	    ReconstructedParticleImpl* tau = new ReconstructedParticleImpl();
	    for ( size_t k=0; k<coneparticles[j].size(); k++) {
	      tau->addParticle( coneparticles[j][k] );
	    }
	    tau->setType(15);
	    TLorentzVector mom = tauUtils::getTLV( coneparticles[j] );
	    float threemom[3];
	    threemom[0]=mom.X();
	    threemom[1]=mom.Y();
	    threemom[2]=mom.Z();
	    tau->setMomentum(threemom);
	    tau->setEnergy( mom.E() );
	    tau->setMass( mom.M() );
	    newcol->addElement( tau );
	  }
	  evt->addCollection( newcol, "danielKeitaRecoTaus" );
	}

      } // second seeds found
    } // first seed

  } catch(DataNotAvailableException &e) {};

  return;
}


void danielKeitaTauFinderProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void danielKeitaTauFinderProcessor::end(){

  std::ofstream resultsfile;
  resultsfile.open("results_" + _outfile + ".txt");

  resultsfile << "danielKeitaTauFinderProcessor SUMMARY" << endl;

  _fout->cd();

//  TH1F* eff_ttmass = (TH1F*) h_ttmass -> Clone("seleff_ttmass");
//  eff_ttmass->Sumw2();
//  eff_ttmass->Divide( hSEL_ttmass, h_ttmass, 1, 1, "B" );
//  
//  for (int i=0; i<NCLASS; i++) {
//    TH1D* hh1 = h_mcTau_costh[i]->ProjectionX();
//    TH1D* hh2 = hSEL_mcTau_costh[i]->ProjectionX();
//    TH1D* heff = (TH1D*) hh1->Clone( "seleff_"+TString(hh1->GetName()) );
//    heff->Sumw2();
//    heff->Divide( hh2, hh1, 1, 1, "B" );
//
//    TH1F* eff_dec = (TH1F*) h_dec[i] -> Clone("seleff_"+TString(h_dec[i]->GetName()) );
//    eff_dec->Sumw2();
//    eff_dec->Divide( hSEL_dec[i], h_dec[i], 1, 1, "B" );
//
//  }

  _fout->Write(0);
  _fout->Close();

  resultsfile << setw(20) << "";
  for (int i=0; i<NCLASS; i++)
    resultsfile << setw(8) << classLabels[i] << " ";
  resultsfile << " : " <<  setw(8) << "TOTAL";
  resultsfile << "  |  ";
  for (int i=0; i<NCLASS; i++)
    resultsfile << setw(8) << classLabels[i] << " ";

  resultsfile << " : " <<  setw(8) << "TOTAL";
  resultsfile << endl;

  resultsfile << setw(20) << "_nOrig"; printEff(_nOrig, NULL, resultsfile);
  resultsfile << setw(20) << "_nPresel"; printEff(_nPresel, _nOrig, resultsfile);
  resultsfile << "----------------------" << endl;
  resultsfile << setw(20) << "_nTwoSeeds";   printEff( _nTwoSeeds, _nOrig , resultsfile);
  resultsfile << setw(20) << "_nGoodSeedDir";  printEff( _nGoodSeedDir, _nOrig , resultsfile);
  resultsfile << "----------------------" << endl;
  resultsfile <<  setw(20) << "_nSel_oocone";  printEff( _nSel_outofcone, _nOrig , resultsfile);
  resultsfile <<  setw(20) << "_nSel_acoLin"; printEff( _nSel_acoLin, _nOrig , resultsfile);
  resultsfile <<  setw(20) << "_nSel_acoPlan"; printEff( _nSel_acoPlan, _nOrig , resultsfile);
  resultsfile <<  setw(20) << "_nSel_isr"; printEff( _nSel_isr, _nOrig , resultsfile);
  resultsfile <<  setw(20) << "_nSel_leptop"; printEff( _nSel_lepton, _nOrig , resultsfile);
  resultsfile <<  setw(20) << "_nSel_seedcl"; printEff( _nSel_seedCluster, _nOrig , resultsfile);
  resultsfile <<  setw(20) << "_nSel_tmass";  printEff( _nSel_tmass, _nOrig , resultsfile);
  resultsfile <<  setw(20) << "_nSel_chg"; printEff( _nSel_chg, _nOrig , resultsfile);
  //  resultsfile <<  setw(20) << "_nSel"; printEff( _nSel, _nOrig , resultsfile);
  resultsfile << "----------------------" << endl;
  resultsfile <<  setw(20) << "_nCumulSel_oocone"; printEff( _nCumulSel_outofcone, _nOrig , resultsfile);
  resultsfile <<  setw(20) << "_nCumulSel_acoLin"; printEff( _nCumulSel_acoLin, _nOrig , resultsfile);
  resultsfile <<  setw(20) << "_nCumulSel_acoPlan"; printEff( _nCumulSel_acoPlan, _nOrig , resultsfile);
  resultsfile <<  setw(20) << "_nCumulSel_isr"; printEff( _nCumulSel_isr, _nOrig , resultsfile);
  resultsfile <<  setw(20) << "_nCumulSel_lepton"; printEff( _nCumulSel_lepton, _nOrig , resultsfile);
  resultsfile <<  setw(20) << "_nCumulSel_seedcl"; printEff( _nCumulSel_seedCluster, _nOrig , resultsfile);
  resultsfile <<  setw(20) << "_nCumulSel_tmass"; printEff( _nCumulSel_tmass, _nOrig , resultsfile);
  resultsfile <<  setw(20) << "_nCumulSel_chg"; printEff( _nCumulSel_chg, _nOrig , resultsfile);
  resultsfile <<  setw(20) << "_nCumulSel"; printEff( _nCumulSel, _nOrig , resultsfile);
  
  resultsfile.close();

  std::cout << "n gamma converted / not converted in trk = " << _nconvertedInTrk << " " << _nNotconvertedInTrk << endl;


  std::cout << "danielKeitaTauFinderProcessor::end()  " << name() 
 	    << std::endl ;
}

void danielKeitaTauFinderProcessor::printEff(int* nsel, int* norig, std::ofstream & resultsfile) {

  if ( !nsel ) return;

  int ntot(0);
  int ntotOrig(0);
  for (int i=0; i<NCLASS; i++) {
    ntot+=nsel[i];
    if ( norig ) ntotOrig+=norig[i];
  }

  resultsfile << std::showpoint;
  resultsfile << std::setprecision(5);

  for (int i=0; i<NCLASS; i++) resultsfile << setw(8) << nsel[i] << " "; 
  resultsfile << " : " <<  setw(8) << ntot;
  resultsfile << "  |  ";
  if ( norig ) {
    for (int i=0; i<NCLASS; i++) resultsfile << setw(8) << 100.*nsel[i]/norig[i] << " ";
    resultsfile << " : " <<  setw(8) << 100.*ntot/ntotOrig;
  }
  resultsfile << endl;
  return;
}

std::vector < std::pair <TVector3, TVector3> > danielKeitaTauFinderProcessor::getNuMom( TLorentzVector jet1,  TLorentzVector jet2, float ecom, int & iflag ) {

  std::vector < std::pair <TVector3, TVector3> > neutrino_solutions;

  //  cout << "hello from getNuMom" << endl;

  // assume back-to-back all energy in COM
  // jet1, 2 are 4-momenta of hadronic momentum
  // ecom is assumed COM energy
  const double mtau=1.777; // GeV/c2
  double ebeam = ecom/2.;

  double nuen1 = ebeam - jet1.E();
  double nuen2 = ebeam - jet2.E();

  if ( nuen1<0 || nuen2<0 ) {

    cout << "WARNING, negative neutrino energy: " << ebeam << " " << jet1.E() << " " << nuen1 << " ; " << jet2.E() << " " << nuen2 << endl;
    iflag=1;

  } else {

    // angles between neutrino1 direction and jet1 dir
    double cosTheta1 = ( 1./( 2.*nuen1*jet1.Vect().Mag()) ) * ( pow(ebeam,2) - pow( mtau,2 ) - jet1.Vect().Mag2() - pow(nuen1, 2) );
    // angles between neutrino1 direction and -jet2 dir
    double cosTheta2 = ( 1./( 2.*nuen1*jet2.Vect().Mag()) ) * ( pow(ebeam,2) - pow( mtau,2 ) + jet2.Vect().Mag2() + 2.*jet1.Vect().Dot(jet2.Vect()) - pow(nuen2,2) );

    if ( fabs( cosTheta1 ) > 1. || fabs( cosTheta2 ) > 1. ) {

      iflag=2; // this means that we cannot make a tau wuth required energy and mass by adding a single neutrino

    } else {


      double theta1=acos(cosTheta1);
      double theta2=acos(cosTheta2);

      // angle between jet1 and -jet2
      double theta12 = jet1.Vect().Angle( -1.*jet2.Vect() );

      if ( theta12 > theta1+theta2 ) {

	iflag=3; // this means the 2 cones do not intersect


      } else {

	// rotate so that jet1 is along z
	TVector3 zaxis(0,0,1);

	TVector3 d1Rot = jet1.Vect(); d1Rot*= 1./d1Rot.Mag(); // jet1 direction
	TVector3 d2Rot = -jet2.Vect(); d2Rot*= 1./d2Rot.Mag(); // -jet2 direction

	TVector3 axis = d1Rot.Cross( zaxis ); // axis and angle of rotation to bring jet1 along z axis
	double angle = d1Rot.Angle( zaxis );
      
	d1Rot.Rotate( angle, axis ); // rotate both vectors
	d2Rot.Rotate( angle, axis );

	// in rotated frame, write intersection direction n as
	//   n = ( sin(theta1) cos(phi), sin(theta1) sin(phi), cos(theta1) )   [angle phi is unknown]
	// know that n Dot -jet2dir = cosTheta2
	// cosTheta2 = sin(theta1) cos(phi) * d2RotX + sin(theta1) sin(phi) * d2RotY + cos(theta1) * d2RotZ
	// sin(phi) * sin(theta1)*d2RotY + cos(phi) * sin(theta1)*d2RotX = cosTheta2 - cos(theta1)*d2RotZ
	// A sin(phi) + B cos(phi) = C
	
	double A = sin(theta1)*d2Rot.Y();
	double B = sin(theta1)*d2Rot.X();
	double C = cosTheta2 - cosTheta1*d2Rot.Z();


	if ( fabs( C/sqrt( A*A + B*B ) ) > 1. ) {

	  iflag=4; // what does this mean???


	} else {

	  iflag=0;

	  double aa = asin( C/sqrt( A*A + B*B ) );
	  double bb = atan2 ( B, A );

	  // 2 solutions
	  for (int isol=0; isol<2; isol++) {
	    double phi = isol==0 ? 
	      aa - bb : 
	      (TMath::Pi() - aa)-bb;

	    TVector3 n1 ( sin(theta1)*cos(phi), sin(theta1)*sin(phi), cos(theta1) );

	    n1.Rotate( -angle, axis );

	    TVector3 calculatedNeutrino1=n1; calculatedNeutrino1*=nuen1/calculatedNeutrino1.Mag();
	    TVector3 calculatedNeutrino2 = - ( calculatedNeutrino1 + jet1.Vect() + jet2.Vect() );
	    TLorentzVector nu1; nu1.SetVectM( calculatedNeutrino1, 0 );
	    TLorentzVector nu2; nu2.SetVectM( calculatedNeutrino2, 0 );

	    neutrino_solutions.push_back( std::pair <TVector3, TVector3> ( calculatedNeutrino1, calculatedNeutrino2 ) );

	  }
	  
	}

      }

    }

  }

  return neutrino_solutions; // std::pair <TVector3, TVector3> ( TVector3(0,0,0), TVector3(0,0,0) );
    
}
