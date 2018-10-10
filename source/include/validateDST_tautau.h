#ifndef validateDST_TauTauProcessor_h
#define validateDST_TauTauProcessor_h 1
#include "marlin/Processor.h"

#include "EVENT/MCParticle.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TVector3.h"

#include "tautau2fEventFitter.h"

using namespace marlin ;


class validateDST_TauTauProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new validateDST_TauTauProcessor ; }

  validateDST_TauTauProcessor( const validateDST_TauTauProcessor&) = delete ;
  validateDST_TauTauProcessor& operator=(const validateDST_TauTauProcessor& ) = delete;

  validateDST_TauTauProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
    
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 

  virtual void check( LCEvent * evt ) ; 
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;


 protected:

  TVector3 getTauolaPolMC(MCParticle* tau);

  int getTauDecayMode_MC( MCParticle* mctau );
  int gettype( int pdg );
  float getTauPiDecayAngle_MC( MCParticle* mctau );
  float getTauPiEnFrac_MC( MCParticle* mctau );
  float getTauPiEn_MC( MCParticle* mctau );
  std::vector <MCParticle*> getTauDaughters_byPdg_MC( MCParticle* mctau, int pdg=0, bool bothCharges=true );
  std::vector <MCParticle*> getTauGrandDaughters_byPdg_MC( MCParticle* mctau, int pdg=0, bool bothCharges=true );
  std::vector <const MCParticle*> getStableTauDescendents( const MCParticle* mctau );
  bool dau_created_in_sim(  const MCParticle* mcp );


  void mcana( std::vector <MCParticle*> cheatPart );


  enum {tauDec_mu=0, tauDec_el, tauDec_a1chg_3p, tauDec_a1chg_1p, tauDec_rho, tauDec_hadW, tauDec_kchg, tauDec_pi, NTAUDEC};
  enum {type_mu=0, type_el, type_gam, type_chad, type_nhad, type_kshort, type_lambda, NTYPE};


  enum {pfoType_Photon=0, pfoType_ChHad, pfoType_NeuHad, pfoType_Mu, pfoType_El, pfoType_Lambda, pfoType_KShort, NPFOTYPE };


  void printDecayMode( int imode );

  int _isPythiaTauDecay;
  int _isDBD;
  int _mcskim;
  int iplot;

  tautau2fEventFitter* _eventFitter;

  StringVec _infiles{};
  std::string _outfile{};
  TFile* _fout = nullptr ;

  TH1F* _h_eventCounter = nullptr ;

  TH1F* _h_mcHas94 = nullptr ;
  TH2F*  _h_mcTauPlusVtxRZ =  nullptr;
  TH2F*  _h_mcTauPlusVtxXY =  nullptr;
  TH2F*  _h_mcTauMinusVtxRZ = nullptr;
  TH2F*  _h_mcTauMinusVtxXY = nullptr;

  TH1F* _h_mcTauPlusDecl = nullptr;
  TH1F* _h_mcTauMinusDecl = nullptr;
  TH1F* _h_mcTauPlusLife = nullptr;
  TH1F* _h_mcTauMinusLife = nullptr;

  TH1F* _h_mcTauTauMass = nullptr ;
  TH1F* _h_mcMuMuMass = nullptr ;
  TH1F* _h_mcTauTauMassSel[2] = {nullptr} ;
  TH1F* _h_mcDecMode = nullptr ;
  TH2F* _h_mcDecMode2d = nullptr ;

  TH1F* _h_mctaunu = nullptr ;

  TH1F* _h_mctoten_mumu = nullptr ;
  TH1F* _h_mctoten_tautau = nullptr ;

  TH1F* _h_pfototen_mumu = nullptr ;
  TH1F* _h_pfototen_tautau = nullptr ;

  TH2F* _h_mcPiPiDecAng[2] = {nullptr} ;
  TH1F* _h_mcPiPiEnFrac[2] = {nullptr} ;
  
  TH2F* _h_mcTauEn[2] = {nullptr} ;

  TH2F* _h_mc_aColin = nullptr ;
  TH2F* _h_mc_aCoplan = nullptr ;

  TH2F* _h_mcVis_aColin = nullptr ;
  TH2F* _h_mcVis_aCoplan = nullptr ;
  TH1F* _h_mcVis_invMass = nullptr ;

  TH1F* _h_mcCheckTauMass = nullptr ;

  TH2F* _h_mcPiPi_PiEnFrac[2] = {nullptr} ;
  TH2F* _h_mcPiPi_PiEn[2] = {nullptr} ;

  TH2F* _h_mcPiPi_polThTh[2] = {nullptr} ;
  TH2F* _h_mcPiRho_polThTh[2] = {nullptr} ;
  TH2F* _h_mcRhoPi_polThTh[2] = {nullptr} ;
  TH2F* _h_mcRhoRho_polThTh[2] = {nullptr} ;

  TH2F* _h_mcRhoRho_colin_polThTh[2] = {nullptr} ;

  TH1F* _h_mcPiPi_polComb1[2] = {nullptr} ;
  TH1F* _h_mcPiRho_polComb1[2] = {nullptr} ;
  TH1F* _h_mcRhoPi_polComb1[2] = {nullptr} ;
  TH1F* _h_mcRhoRho_polComb1[2] = {nullptr} ;

  TH1F* _h_mcPiPi_polComb2[2] = {nullptr} ;
  TH1F* _h_mcPiRho_polComb2[2] = {nullptr} ;
  TH1F* _h_mcRhoPi_polComb2[2] = {nullptr} ;
  TH1F* _h_mcRhoRho_polComb2[2] = {nullptr} ;


  TH2F* _h_mcPfoMatch_tauPi_nChNeu = nullptr;
  TH1F* _h_mcPfoMatch_tauPi_nChpi =  nullptr;
  TH1F* _h_mcPfoMatch_tauPi_nNeuh =  nullptr;
  TH1F* _h_mcPfoMatch_tauPi_nGamma =  nullptr;
  TH1F* _h_mcPfoMatch_tauPi_nEl =  nullptr;
  TH1F* _h_mcPfoMatch_tauPi_nMu =  nullptr;
  TH1F* _h_mcPfoMatch_tauPi_nV0 =  nullptr;

  TH2F*  _h_mcPfoMatch_tauRho_nChNeu = nullptr;
  TH1F*  _h_mcPfoMatch_tauRho_nChpi = nullptr;
  TH1F*  _h_mcPfoMatch_tauRho_nNeuh = nullptr;
  TH1F*  _h_mcPfoMatch_tauRho_nGamma = nullptr;
  TH1F*  _h_mcPfoMatch_tauRho_nEl = nullptr;
  TH1F*  _h_mcPfoMatch_tauRho_nMu = nullptr;
  TH1F*  _h_mcPfoMatch_tauRho_nV0 = nullptr;

  TH1F*  _h_mcPfoMatch_chgz0 = nullptr;
  TH1F*  _h_mcPfoMatch_chgz0WrtPv = nullptr;

  TH2F* _h_evnPt[4];
  TH2F* _h_evnPz[4];
  TH2F* _h_evnISRMass[4];
  TH2F* _h_evnE[4];

  TH1F* _h_mcfitall_ttpz = nullptr;
  TH1F* _h_mcfitall_ttE = nullptr;
  TH1F* _h_mcfitall_ttpt = nullptr;


  TH2F* _h_mcTauCosth[2] ;
  TH1F* _h_nPFO[2] ;
  TH1F* _h_pfo_id[2] ;
  TH2F* _h_pfos_by_type[2] ;

  TH1F* _h_pfos_by_dec_type[2][NTAUDEC][NTYPE] ;

  TH1F* _h_visMass[2][NTAUDEC];
  TH1F* _h_ggMass[2][NTAUDEC];
  TH1F* _h_npfos[2][NTAUDEC];

  TH2F* _h_pfotauang[2];
  TH1F* _h_ntrk[2];
  TH1F* _h_trkD0[2];
  TH1F* _h_trkD0err[2];
  TH1F* _h_trkD0sig[2];

  TH1F* _h_mumu_nPFO        ;
  TH1F* _h_mumu_nPFOChg	;
  TH1F* _h_mumu_nPFONeu     ;
  TH1F* _h_mumu_pfoNeu_costh;
  TH1F* _h_mumu_pfoChg_costh;
  TH1F* _h_mumu_pfoNeu_en   ;
  TH1F* _h_mumu_pfoChg_en   ;

  TH1F* _h_mumu_pfoMCchad_costh;
  TH1F* _h_mumu_pfoMCnhad_costh;
  TH1F* _h_mumu_pfoMCgam_costh;
  TH1F* _h_mumu_pfoMCele_costh;
  TH1F* _h_mumu_pfoMCmuo_costh;
  TH1F* _h_mumu_pfoMCoth_costh;
  TH1F* _h_mumu_trk_costh;
  TH1F* _h_mumu_mcele_costh;

  TH1F* _mc_taujet_np;
  TH1F* _mc_taujet_m;
  TH1F* _mc_taujet_e;


  TH1F* _h_cpcheck;
  TH1F* _h_cpcheck_leplep;
  TH1F* _h_cpcheck_pipi;
  TH1F* _h_cpcheck_rhorho;
  TH1F* _h_cpcheck_a3a3;
  TH1F* _h_cpcheck_a1a1;

  TH2F* _h_cpcheckPolAngle_a3a3;
  TH1F* _h_cpcheckDPhiDecPlane_a3a3;
  TH1F* _h_cpcheckDPhiDecPlane_a3a3_r1;
  TH1F* _h_cpcheckDPhiDecPlane_a3a3_r2;
  TH1F* _h_cpcheckDPhiDecPlane_a3a3_r3;

  TH2F* _h_cpcheckPolAngle_pipi;
  TH1F* _h_cpcheckDPhiDecPlane_pipi;
  TH1F* _h_cpcheckDPhiDecPlane_pipi_r1;
  TH1F* _h_cpcheckDPhiDecPlane_pipi_r2;
  TH1F* _h_cpcheckDPhiDecPlane_pipi_r3;

  TH2F* _h_cpcheckPolAngle_rhorho;
  TH1F* _h_cpcheckDPhiDecPlane_rhorho;
  TH1F* _h_cpcheckDPhiDecPlane_rhorho_r1;
  TH1F* _h_cpcheckDPhiDecPlane_rhorho_r2;
  TH1F* _h_cpcheckDPhiDecPlane_rhorho_r3;

  TH2F* _h_cpcheck2_leplep;
  TH2F* _h_cpcheck2_pipi;
  TH2F* _h_cpcheck2_rhorho;
  TH2F* _h_cpcheck2_a3a3;
  TH2F* _h_cpcheck2_a1a1;

};


#endif



