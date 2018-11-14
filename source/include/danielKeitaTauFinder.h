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

  float _conesize;

  int _nNotconvertedInTrk, _nconvertedInTrk;

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


  TH2F* h_tchg_tmass[4];
  TH2F* h_ngam_tmass[4];
  TH2F* h_nchg_tmass[4];
  TH2F* h_nnhad_tmass[4];
  TH2F* h_ngam_nchg[4];
  TH2F* h_ttmass_prongAngle[4];
  TH2F* h_ttmass_outsideEnergy[4];
  TH2F* h_ttmass_insideEnergy[4];
  TH2F* h_ttmass_outsidePt[4];

  TH2F* h_mcTau_costh[4];
  TH2F* hSEL_mcTau_costh[4];

  TH2F* h_dec_tmass[4];
  TH2F* h_dec_ngam[4];
  TH2F* h_dec_nchg[4];
  TH1F* h_dec[4];
  TH1F* hSEL_dec[4];


  TH2F* hSEL_tchg_tmass[4];
  TH2F* hSEL_ngam_tmass[4];
  TH2F* hSEL_nchg_tmass[4];
  TH2F* hSEL_nnhad_tmass[4];
  TH2F* hSEL_ngam_nchg[4];
  TH2F* hSEL_ttmass_prongAngle[4];
  TH2F* hSEL_ttmass_outsideEnergy[4];
  TH2F* hSEL_ttmass_insideEnergy[4];
  TH2F* hSEL_ttmass_outsidePt[4];

  TH2F* h_a3p_cone_chgMC_pfo[4];
  TH2F* h_a3p_cone_chgMC_trk[4];
  TH2F* h_a3p_cone_trk_pfo  [4];

  int _nOrig[4];
  int _nTwoSeeds[4];
  int _nTwoWellMatchedSeeds[4];
  int _nGoodSeedDir[4];
  int _nSel_tmass[4];
  int _nSel_outofcone[4];
  int _nSel_angle[4];
  int _nSel_chg[4];
  int _nSel_dilep[4];
  int _nSel[4];

};


#endif



