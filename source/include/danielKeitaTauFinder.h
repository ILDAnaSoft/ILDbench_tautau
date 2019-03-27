#ifndef danielKeitaTauFinderProcessor_h
#define danielKeitaTauFinderProcessor_h 1
#include "marlin/Processor.h"

#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

#include "tauUtils.h"
#include "tautau2fEventFitter.h"

#include <fstream>

using namespace marlin ;


class danielKeitaTauFinderProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new danielKeitaTauFinderProcessor ; }
  
  danielKeitaTauFinderProcessor() ;
  
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

  bool isNeutral( MCParticle* mcp ) { return fabs( mcp->getCharge() ) < 0.1; }
  bool isNeutral( ReconstructedParticle* rp ) { return fabs( rp->getCharge() ) < 0.1; }
  bool isCharged( MCParticle* mcp ) { return !isNeutral(mcp) ; }
  bool isCharged( ReconstructedParticle* rp ) { return !isNeutral( rp ); }

  void printEff(int* nsel, int* norig, std::ofstream & resultsfile);
  std::vector < std::pair <TVector3, TVector3> > getNuMom( TLorentzVector jet1,  TLorentzVector jet2, float ecom, int & iflag );

  float _conesize;

  int _nNotconvertedInTrk, _nconvertedInTrk;

  bool _highestPt;
  bool _useDistilled;

  tautau2fEventFitter* _2tauEvtFitter;

  std::string _outfile;
  TFile* _fout;
  TH1F* h_gamgam_mass;

  TH1F* h_ttmass;
  TH1F* hSEL_ttmass;

  TH2F* h_conversionPos;
  TH2F* h_v0Pos;

  TH2F* h_conversionPos2;
  TH2F* h_v0Pos2;

  TH1F* h_neuPandoraCompound_mass;
  TH1F* h_neuPandoraCompound_ntrk;
  TH2F* h_neuPandoraCompound_vtxPos;
  TH2F* h_neuPandoraCompound_vtxPos2;
  TH1F* h_neuPandoraCompound_vtxChisq;

  TH2F* h_neutronClus_eOnP_beforeAfter;

  TH2F* h_neutralHadronPFO_mainMCcontrib_pdgEn;
  TH1F* h_rhoDecaySinglePhoClus_reason;
  TH2F* h_rhoDecaySinglePhoClus_mergedGammaClusterEvals;
  TH2F* h_rhoDecaySinglePhoClus_singleGammaClusterEvals;

  enum {NCLASS=7, NDEC=8};

  std::string classLabels[NCLASS];

  TH2F* h_mctautau_ecom_tauMinusCosth;

  TH2F* h_pirho_tauMinusCosth_mcPolarExact[NCLASS];
  TH2F* h_pirho_mcPolarExactPlusMinus[NCLASS];
  TH2F* h_pirho_mcPolarApproxPlusMinus[NCLASS];

//  TH1F*  h_rho_mcapprox_cospsi_helPos[NCLASS];
//  TH1F*  h_rho_mcapprox_costheta_helPos[NCLASS];
//  TH1F*  h_rho_mcapprox_cosbeta_helPos[NCLASS];
//
//  TH1F*  h_rho_mcapprox_cospsi_helNeg[NCLASS];
//  TH1F*  h_rho_mcapprox_costheta_helNeg[NCLASS];
//  TH1F*  h_rho_mcapprox_cosbeta_helNeg[NCLASS];


  TH1F* h_nIsoMuons[NCLASS];
  TH1F* h_nIsoElectrons[NCLASS];

  TH2F* h_rho_mcPolar[NCLASS];
  TH2F* h_rho_mcPolar2[NCLASS];
  TH1F* h_pi_mcPolar[NCLASS];

  TH1F* h_rho_mcPolarExact[NCLASS];
  TH1F* h_pi_mcPolarExact[NCLASS];


  TH1F* h_rho_mcPolarApprox[NCLASS];
  TH1F* h_pi_mcPolarApprox[NCLASS];
  
  TH1F* h_rho_mcPolarApproxCheatEn[NCLASS];
  TH1F* h_pi_mcPolarApproxCheatEn[NCLASS];
  
  TH2F* h_rho_mcPolarExactApprox[NCLASS];

  TH1F* h_rho_mcPolarExact_helNeg[NCLASS];
  TH1F* h_rho_mcPolarExact_helPos[NCLASS];
  TH1F*  h_pi_mcPolarExact_helNeg[NCLASS];
  TH1F*  h_pi_mcPolarExact_helPos[NCLASS];

  TH1F* hSEL_rho_mcPolarExact_helNeg[NCLASS];
  TH1F* hSEL_rho_mcPolarExact_helPos[NCLASS];
  TH1F* hSEL_pi_mcPolarExact_helNeg[NCLASS];
  TH1F* hSEL_pi_mcPolarExact_helPos[NCLASS];

  TH1F* h_rho_mcPolarApprox_helNeg[NCLASS];
  TH1F* h_rho_mcPolarApprox_helPos[NCLASS];
  TH1F*  h_pi_mcPolarApprox_helNeg[NCLASS];
  TH1F*  h_pi_mcPolarApprox_helPos[NCLASS];

  TH1F* hSEL_rho_mcPolarApprox_helNeg[NCLASS];
  TH1F* hSEL_rho_mcPolarApprox_helPos[NCLASS];
  TH1F* hSEL_pi_mcPolarApprox_helNeg[NCLASS];
  TH1F* hSEL_pi_mcPolarApprox_helPos[NCLASS];

  TH1F* h_rho_mcPolarApproxCheatEn_helNeg[NCLASS];
  TH1F* h_rho_mcPolarApproxCheatEn_helPos[NCLASS];
  TH1F*  h_pi_mcPolarApproxCheatEn_helNeg[NCLASS];
  TH1F*  h_pi_mcPolarApproxCheatEn_helPos[NCLASS];

  TH2F* h_rho_mcPolarComp[NCLASS];
  TH2F* h_pi_mcPolarComp[NCLASS];

  TH1F* h_mc_tauMinus_costh[NCLASS];
  TH1F* hSEL_mc_tauMinus_costh[NCLASS];
  TH2F* h_mc_tauSpin[NCLASS];

  TH2F* h_npfo_chg_nneu[NCLASS];


  TH1F*  h_gammaClus_nMajor;
  TH1F*  h_gammaClus_nSignf;
  TH2F*  h_gammaClus_nSignf_bMass;
  TH2F*  h_gammaClus_nSignf_bDist;
  TH1F*  h_neutronClus_nMajor;
  TH1F*  h_neutronClus_nSignf;
  TH2F*  h_neutronClus_nSignf_bMass;
  TH2F*  h_neutronClus_nSignf_bDist;



  TH2F* h_tchg_tmass[NCLASS];

  TH2F* h_rawmass_trimmass[NCLASS];

  TH2F* h_ngam_tmass[NCLASS];
  TH2F* h_nchg_tmass[NCLASS];
  TH2F* h_nnhad_tmass[NCLASS];
  TH2F* h_ngam_nchg[NCLASS];
  TH1F* h_prongAngle[NCLASS];
  TH1F* h_outsideEnergy[NCLASS];
  TH1F* h_insideEnergy[NCLASS];
  TH1F* h_outsidePt[NCLASS];

  TH2F* h_mcTau_costh[NCLASS];
  TH2F* hSEL_mcTau_costh[NCLASS];

  TH2F* h_dec_tmass[NCLASS];
  TH2F* h_dec_ngam[NCLASS];
  TH2F* h_dec_nchg[NCLASS];
  TH1F* h_dec[NCLASS];
  TH1F* hSEL_dec[NCLASS];

  TH2F* hSEL_tchg_tmass[NCLASS];
  TH2F* hSEL_ngam_tmass[NCLASS];
  TH2F* hSEL_nchg_tmass[NCLASS];
  TH2F* hSEL_nnhad_tmass[NCLASS];
  TH2F* hSEL_ngam_nchg[NCLASS];
  TH1F* hSEL_prongAngle[NCLASS];
  TH1F* hSEL_outsideEnergy[NCLASS];
  TH1F* hSEL_insideEnergy[NCLASS];
  TH1F* hSEL_outsidePt[NCLASS];

  //  TH2F* h_a3p_cone_chgMC_pfo[NCLASS];
  //  TH2F* h_a3p_cone_chgMC_trk[NCLASS];
  //  TH2F* h_a3p_cone_trk_pfo  [NCLASS];

  TH1F* h_dec_cone_nchpfo[NCLASS][NDEC];
  TH1F* h_dec_cone_ntracks[NCLASS][NDEC];
  TH1F* h_dec_cone_trackNTPC[NCLASS][NDEC];
  TH1F* h_dec_cone_nnhadpfo[NCLASS][NDEC];
  TH1F* h_dec_cone_ngammapfo[NCLASS][NDEC];
  TH1F* h_dec_cone_npi0pfo[NCLASS][NDEC];
  TH1F* h_dec_cone_ncompoundpfo[NCLASS][NDEC];
  TH1F* h_dec_cone_nhaden[NCLASS][NDEC];
  TH1F* h_dec_cone_gammaen[NCLASS][NDEC];
  TH1F* h_dec_cone_nhadenFrac[NCLASS][NDEC];
  TH1F* h_dec_cone_gammaenFrac[NCLASS][NDEC];

  TH1F* h_dec_seed_clusterWidth1[NCLASS][NDEC];
  TH1F* h_dec_seed_clusterWidth2[NCLASS][NDEC];
  TH1F* h_dec_seed_clusterLength[NCLASS][NDEC];
  TH1F* h_dec_seed_energy[NCLASS][NDEC];

  TH1F* hSEL_dec_seed_clusterWidth1[NCLASS][NDEC];
  TH1F* hSEL_dec_seed_clusterWidth2[NCLASS][NDEC];
  TH1F* hSEL_dec_seed_clusterLength[NCLASS][NDEC];
  TH1F* hSEL_dec_seed_energy[NCLASS][NDEC];

  TH1F* h_dec_cone_visMass[NCLASS][NDEC];
  TH1F* h_dec_cone_neutralvisMass[NCLASS][NDEC];
  TH1F* h_dec_cone_visMassDiff[NCLASS][NDEC];
  TH1F* h_dec_cone_neutralvisMassDiff[NCLASS][NDEC];
  TH1F* h_dec_cone_visEnergyDiff[NCLASS][NDEC];
  TH1F* h_dec_cone_neutralvisEnergyDiff[NCLASS][NDEC];

//  TH2F* h_dec_coneTRIM_vis_neutral_Mass[NCLASS][NDEC];
//  TH2F* h_dec_coneTRIM_vis_neutral_Mass_0g[NCLASS][NDEC];
//  TH2F* h_dec_coneTRIM_vis_neutral_Mass_1g[NCLASS][NDEC];
//  TH2F* h_dec_coneTRIM_vis_neutral_Mass_2g[NCLASS][NDEC];
//  TH2F* h_dec_coneTRIM_vis_neutral_Mass_3g[NCLASS][NDEC];
//  TH2F* h_dec_coneTRIM_vis_neutral_Mass_4g[NCLASS][NDEC];

  TH2F* h_dec_cone_vis_neutral_Mass[NCLASS][NDEC];
  TH2F* h_dec_cone_vis_neutral_Mass_0g[NCLASS][NDEC];
  TH2F* h_dec_cone_vis_neutral_Mass_1g[NCLASS][NDEC];
  TH2F* h_dec_cone_vis_neutral_Mass_2g[NCLASS][NDEC];
  TH2F* h_dec_cone_vis_neutral_Mass_3g[NCLASS][NDEC];
  TH2F* h_dec_cone_vis_neutral_Mass_4g[NCLASS][NDEC];

  TH2F* hSEL_dec_coneTRIM_vis_neutral_Mass[NCLASS][NDEC];
  TH2F* hSEL_dec_coneTRIM_vis_neutral_Mass_0g[NCLASS][NDEC];
  TH2F* hSEL_dec_coneTRIM_vis_neutral_Mass_1g[NCLASS][NDEC];
  TH2F* hSEL_dec_coneTRIM_vis_neutral_Mass_2g[NCLASS][NDEC];
  TH2F* hSEL_dec_coneTRIM_vis_neutral_Mass_3g[NCLASS][NDEC];
  TH2F* hSEL_dec_coneTRIM_vis_neutral_Mass_4g[NCLASS][NDEC];

  TH2F* hSEL_dec_cone_vis_neutral_Mass[NCLASS][NDEC];
  TH2F* hSEL_dec_cone_vis_neutral_Mass_0g[NCLASS][NDEC];
  TH2F* hSEL_dec_cone_vis_neutral_Mass_1g[NCLASS][NDEC];
  TH2F* hSEL_dec_cone_vis_neutral_Mass_2g[NCLASS][NDEC];
  TH2F* hSEL_dec_cone_vis_neutral_Mass_3g[NCLASS][NDEC];
  TH2F* hSEL_dec_cone_vis_neutral_Mass_4g[NCLASS][NDEC];

  TH1F* h_dec_mc_ngamma               [NCLASS][NDEC];
  TH1F* h_dec_mc_vismass              [NCLASS][NDEC];
  TH1F* h_dec_mc_visneutralmass       [NCLASS][NDEC];
  TH1F* h_dec_mc_cone_ngamma          [NCLASS][NDEC];
  TH1F* h_dec_mc_cone_vismass         [NCLASS][NDEC];
  TH1F* h_dec_mc_cone_visneutralmass  [NCLASS][NDEC];

  TH1F* h_dec_mcall_cone_ngamma          [NCLASS][NDEC];
  TH1F* h_dec_mcall_cone_vismass         [NCLASS][NDEC];
  TH1F* h_dec_mcall_cone_visneutralmass  [NCLASS][NDEC];
 
//  TH1F* h_dec_coneTRIM_visMass[NCLASS][NDEC];
//  TH1F* h_dec_coneTRIM_neutralvisMass[NCLASS][NDEC];
//  TH1F* h_dec_coneTRIM_visMassDiff[NCLASS][NDEC];
//  TH1F* h_dec_coneTRIM_neutralvisMassDiff[NCLASS][NDEC];
//  TH1F* h_dec_coneTRIM_visEnergyDiff[NCLASS][NDEC];
//  TH1F* h_dec_coneTRIM_neutralvisEnergyDiff[NCLASS][NDEC];


  TH1F* hSEL_dec_cone_ntracks[NCLASS][NDEC];
  TH1F* hSEL_dec_cone_trackNTPC[NCLASS][NDEC];
  TH1F* hSEL_dec_cone_nchpfo[NCLASS][NDEC];
  TH1F* hSEL_dec_cone_nnhadpfo[NCLASS][NDEC];
  TH1F* hSEL_dec_cone_ngammapfo[NCLASS][NDEC];
  TH1F* hSEL_dec_cone_ncompoundpfo[NCLASS][NDEC];
  TH1F* hSEL_dec_cone_nhaden[NCLASS][NDEC];
  TH1F* hSEL_dec_cone_gammaen[NCLASS][NDEC];
  TH1F* hSEL_dec_cone_nhadenFrac[NCLASS][NDEC];
  TH1F* hSEL_dec_cone_gammaenFrac[NCLASS][NDEC];

  TH1F* hSEL_dec_cone_visMass[NCLASS][NDEC];
  TH1F* hSEL_dec_cone_neutralvisMass[NCLASS][NDEC];
  TH1F* hSEL_dec_cone_visMassDiff[NCLASS][NDEC];
  TH1F* hSEL_dec_cone_neutralvisMassDiff[NCLASS][NDEC];
  TH1F* hSEL_dec_cone_visEnergyDiff[NCLASS][NDEC];
  TH1F* hSEL_dec_cone_neutralvisEnergyDiff[NCLASS][NDEC];

  TH1F* hSEL_dec_coneTRIM_visMass[NCLASS][NDEC];
  TH1F* hSEL_dec_coneTRIM_neutralvisMass[NCLASS][NDEC];
  TH1F* hSEL_dec_coneTRIM_visMassDiff[NCLASS][NDEC];
  TH1F* hSEL_dec_coneTRIM_neutralvisMassDiff[NCLASS][NDEC];
  TH1F* hSEL_dec_coneTRIM_visEnergyDiff[NCLASS][NDEC];
  TH1F* hSEL_dec_coneTRIM_neutralvisEnergyDiff[NCLASS][NDEC];


  TH1F* h_seed_dZ0[NCLASS];
  TH1F* h_seed_dZ0ave[NCLASS];


  TH1F* h_maxSeedEn[NCLASS];
  TH1F* h_minSeedEn[NCLASS];
  TH1F* h_maxSeedEn_eonp[NCLASS];
  TH1F* h_minSeedEn_eonp[NCLASS];
  TH1F* h_maxSeedEn_caloen[NCLASS];
  TH1F* h_minSeedEn_caloen[NCLASS];
  TH1F* h_maxSeedCosth[NCLASS];
  TH1F* h_minSeedCosth[NCLASS];

  TH1F* h_seedClusterWidth1[NCLASS];
  TH1F* h_seedClusterWidth2[NCLASS];
  TH1F* h_seedClusterLength[NCLASS];

  TH1F* h_seedEcalEn[NCLASS];
  TH1F* h_seedEcalEonP[NCLASS];
  TH1F* h_ooconeMaxGammaEn[NCLASS];
  TH1F* h_ooconeMaxPFOEn[NCLASS];
  TH1F* h_coneMass[NCLASS];
  TH1F* h_seedSeedAngle[NCLASS];
  TH1F* h_seedSeedDphi[NCLASS];
  TH1F* h_visjetAcolinearity[NCLASS];
  TH1F* h_visjetAcoplanarity[NCLASS];


  TH1F* hSEL_maxSeedEn[NCLASS];
  TH1F* hSEL_minSeedEn[NCLASS];
  TH1F* hSEL_maxSeedEn_eonp[NCLASS];
  TH1F* hSEL_minSeedEn_eonp[NCLASS];
  TH1F* hSEL_maxSeedEn_caloen[NCLASS];
  TH1F* hSEL_minSeedEn_caloen[NCLASS];
  TH1F* hSEL_maxSeedCosth[NCLASS];
  TH1F* hSEL_minSeedCosth[NCLASS];
  TH1F* hSEL_seedEcalEn[NCLASS];
  TH1F* hSEL_seedEcalEonP[NCLASS];
  TH1F* hSEL_ooconeMaxGammaEn[NCLASS];
  TH1F* hSEL_ooconeMaxPFOEn[NCLASS];
  TH1F* hSEL_coneMass[NCLASS];
  TH1F* hSEL_seedSeedAngle[NCLASS];
  TH1F* hSEL_seedSeedDphi[NCLASS];
  TH1F* hSEL_visjetAcolinearity[NCLASS];
  TH1F* hSEL_visjetAcoplanarity[NCLASS];

  TH1F* hSEL_tauMinusCosth[NCLASS];

  TH1F* hSEL_mcdec[NCLASS];

  TH1F* hSEL_rec_rho_mcdec[NCLASS];
  TH1F* hSEL_rec_rho_pol  [NCLASS];
  TH1F* hSEL_rec_rho_pol_MCpos[NCLASS];
  TH1F* hSEL_rec_rho_pol_MCneg[NCLASS];
  TH1F* hSEL_rec_rho_pol_MCoth[NCLASS];
  TH1F* hSEL_rec_rho_coneMass[NCLASS];

  TH1F* hSEL_rec_a1p_mcdec[NCLASS];
  TH1F* hSEL_rec_a1p_pol  [NCLASS];
  TH1F* hSEL_rec_a1p_coneMass[NCLASS];

  TH1F* hSEL_rec_pi_mcdec[NCLASS];
  TH1F* hSEL_rec_pi_pol  [NCLASS];
  TH1F* hSEL_rec_pi_pol_MCpos  [NCLASS];
  TH1F* hSEL_rec_pi_pol_MCneg  [NCLASS];
  TH1F* hSEL_rec_pi_pol_MCoth  [NCLASS];
  TH1F* hSEL_rec_pi_coneMass[NCLASS];

  TH2F* hSEL_recpi_mcall_pol[NCLASS];
  TH2F* hSEL_recpi_mcpi_pol[NCLASS];

  TH1F* hSEL_recpi_mcall_dpol[NCLASS];
  TH1F* hSEL_recpi_mcpi_dpol[NCLASS];

  TH2F* hSEL_recrho_mcall_pol[NCLASS];
  TH2F* hSEL_recrho_mcrho_pol[NCLASS];

  TH1F* hSEL_recrho_mcall_dpol[NCLASS];
  TH1F* hSEL_recrho_mcrho_dpol[NCLASS];


//  TH2F* _hFit_nSolutions_cone_ip[NCLASS][3];
//  TH1F* _hFit_nSolutions_cone[NCLASS][3];
//  TH1F* _hFit_nSolutions_ip[NCLASS][3];
//
//  TH2F* _hFit_coneFit_pol12[NCLASS][3];
//  TH1F* _hFit_coneFit_pol[NCLASS][3];
//  TH1F* _hFit_coneFit_fitFlag[NCLASS][3];
//  TH2F* _hFit_coneFit_polRecMc[NCLASS][2];
//  TH1F* _hFit_coneFit_polRec[NCLASS][2];
//  TH1F* _hFit_coneFit_polRecMcDiff[NCLASS][2];
//
//  TH2F* _hFit_ipFit_pol12[NCLASS][3];
//  TH1F* _hFit_ipFit_pol[NCLASS][3];
//  TH2F* _hFit_ipFit_polRecMc[NCLASS][2];
//  TH1F* _hFit_ipFit_polRecMcDiff[NCLASS][2];
//  TH1F* _hFit_ipFit_polRec[NCLASS][2];
//
//  TH1F* _hFit_ipFit_allSols_totE[NCLASS][3];
//  TH1F* _hFit_ipFit_allSols_totPt[NCLASS][3];
//  TH1F* _hFit_ipFit_allSols_totPz[NCLASS][3];
//  TH1F* _hFit_ipFit_allSols_decayLength[NCLASS][3];
//
//  TH1F* _hFit_ipFit_reasonableSols_totE[NCLASS][3];
//  TH1F* _hFit_ipFit_reasonableSols_totPt[NCLASS][3];
//  TH1F* _hFit_ipFit_reasonableSols_totPz[NCLASS][3];
//  TH1F* _hFit_ipFit_reasonableSols_decayLength[NCLASS][3];
//
//  TH1F* _hFit_ipFit_bestSol_totE[NCLASS][3];
//  TH1F* _hFit_ipFit_bestSol_totPt[NCLASS][3];
//  TH1F* _hFit_ipFit_bestSol_totPz[NCLASS][3];
//  TH1F* _hFit_ipFit_bestSol_decayLength[NCLASS][3];


  int _nOrig[NCLASS];
  int _nPresel[NCLASS];
  int _nTwoSeeds[NCLASS];
  int _nTwoWellMatchedSeeds[NCLASS];
  int _nGoodSeedDir[NCLASS];
  int _nSel_secondseenen[NCLASS];
  int _nSel_tmass[NCLASS];
  int _nSel_outofcone[NCLASS];
  int _nSel_acoLin[NCLASS];
  int _nSel_acoPlan[NCLASS];
  int _nSel_chg[NCLASS];
  int _nSel_lepton[NCLASS];
  //  int _nSel_seedcaloen[NCLASS];
  int _nSel_seedCluster[NCLASS];
  int _nSel_isr[NCLASS];
  int _nSel[NCLASS];

  int _nCumulSel_secondseenen[NCLASS];
  int _nCumulSel_tmass[NCLASS];
  int _nCumulSel_outofcone[NCLASS];
  int _nCumulSel_acoLin[NCLASS];
  int _nCumulSel_acoPlan[NCLASS];
  int _nCumulSel_chg[NCLASS];
  int _nCumulSel_lepton[NCLASS];
  //  int _nCumulSel_seedcaloen[NCLASS];
  int _nCumulSel_seedCluster[NCLASS];
  int _nCumulSel_isr[NCLASS];
  int _nCumulSel[NCLASS];

  LCRelationNavigator* _relNavi;
  LCRelationNavigator* _relNavi2;


};


#endif



