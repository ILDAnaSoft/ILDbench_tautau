#ifndef _TAUUTILS_H_
#define _TAUUTILS_H_

#include "TH1F.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "EVENT/TrackState.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Cluster.h"
#include "EVENT/LCCollection.h"

#include "UTIL/PIDHandler.h"

#include "UTIL/LCRelationNavigator.h"

#include <cassert>

class tauUtils {

 public:

  static TLorentzVector getFourMomentum( const EVENT::Cluster* cl );
  static TLorentzVector getFourMomentum( const EVENT::TrackState* tst , double mass );
  static TLorentzVector getFourMomentum( const EVENT::MCParticle* rp );
  //  static TVector3 getDvector           ( const TrackState* tst , TVector3 IP );  // from IP to 2d PCA
  static TVector3 getDvector           ( const EVENT::TrackState* tst );
  static TVector3 getDvector           ( const EVENT::MCParticle* mc );
  static void calculateTrackWRTip( const EVENT::TrackState* tst, TVector3& momentum, TVector3& dVec, TVector3& eVec , TVector3 IP );

  static TVector3 getPolarimeter_rho ( TLorentzVector charged, TLorentzVector neutral, TLorentzVector neutrino );
  static TVector3 getPolarimeter_pi  ( TLorentzVector charged, TLorentzVector neutrino );

  static TLorentzVector getPolarimeterTLV_rho ( TLorentzVector charged, TLorentzVector neutral, TLorentzVector neutrino );
  static TLorentzVector getPolarimeterTLV_pi  ( TLorentzVector charged, TLorentzVector neutrino );

  static float getPhiStar ( TLorentzVector tau0, TVector3 pol0, TVector3 pol1 );

  static float getPhiStar ( TLorentzVector tau0, TLorentzVector tau1, TLorentzVector pol0, TLorentzVector pol1 );

  static float getPhiStarRhoRho ( TLorentzVector charged0, TLorentzVector neutral0, TLorentzVector neutrino0,
				  TLorentzVector charged1, TLorentzVector neutral1, TLorentzVector neutrino1 );


  static void getElectronVars( const EVENT::ReconstructedParticle* rp, 
			       float & EonP, // energy / momentum
			       float & ecalFrac, // fraction in ECAL
			       float & epDiffRel, // enery/mom diff in terms of ecal resolution
			       float & eCalo, // total calo energy
			       float & eMuon // yoke energy
			       );


  static int isMuon( const EVENT::ReconstructedParticle* rp );
  static int isElectron( const EVENT::ReconstructedParticle* rp );
  static int isPhoton( const EVENT::ReconstructedParticle* rp );


  static bool isPdg( EVENT::ReconstructedParticle* rp , const int pdg , UTIL::PIDHandler& pidhand );

  enum {NONE=0, LOOSE, TIGHT};


  static float getEcalEnergyFraction( const EVENT::ReconstructedParticle* rp );

  static std::map < const EVENT::ReconstructedParticle* , std::vector < const EVENT::ReconstructedParticle* > > 
    find_FSR_brems( const std::vector < const EVENT::ReconstructedParticle* > & chargedpfos , 
		    const std::vector < const EVENT::ReconstructedParticle* > & photonpfos    ) ;


  static const double bfield;
  static const double crossingAngle;

  static const double omega_pt_conv;
  static const double m_z;
  static const double m_tau;
  static const double m_mu;
  static const double m_el;
  static const double m_piChg;
  static const double m_pi0;
  static const double m_rho;
  static const double g_rho;
  static const double Pi    ;
  static const double twoPi ;
  static const double halfPi;


  static const double ecal_stochastic;
  static const double ecal_constant;
  static const double ctau_tau;

  static void rouge_getThetaDphi( TLorentzVector charged0, TLorentzVector neutral0, TLorentzVector neutrino0,
				  TLorentzVector charged1, TLorentzVector neutral1, TLorentzVector neutrino1,
				  float & theta0, float& theta1, float& dPhi );

  static void gunion_getThetaDphi( TLorentzVector charged0, TLorentzVector neutral0, TLorentzVector neutrino0,
				   TLorentzVector charged1, TLorentzVector neutral1, TLorentzVector neutrino1,
				   float & theta0, float& theta1, float& dPhi );

  static TVector3 getGunionPolarimeter_rho(  TLorentzVector tlchg, TLorentzVector tlpi0, TLorentzVector tlneu);

  static TVector3 getGunionPolarimeter_pi(  TLorentzVector tlchg, TLorentzVector tlneu);

  static TVector3 getGunionAxis( TLorentzVector t1, TLorentzVector t2 );


  static float getRougeLogLikelihood( float theta0, float theta1, float dPhi, float psi , bool useFullFunction=true );


  static void FillLogLike ( TH1F* hLike,
			    TLorentzVector charged0, TLorentzVector neutral0, TLorentzVector neutrino0,
			    TLorentzVector charged1, TLorentzVector neutral1, TLorentzVector neutrino1,
			    bool fullLike=true);

  static void FillLogLike ( TH1F* hLike, float theta0, float theta1, float dPhi,
			    bool useFullFunction=true );

  static TLorentzVector getPolarimeterTLV (  TLorentzVector charged, TLorentzVector neutral, TLorentzVector neutrino, int decayMode );
  static TLorentzVector getPolarimeterTLV (  TLorentzVector charged, TLorentzVector neutral, TLorentzVector neutrino );

  enum { decayChPi=0, decayRho, decayA1_1p, decayA1_3p , decayEl, decayMu , decayW , decayK , decayMultiprong , decayOthersingleprong, decayUnrecognised, NDECAYS};

  //  enum {TDEC_E=0,TDEC_M,TDEC_PI,TDEC_RHO,TDEC_A3P,TDEC_A1P,TDEC_HAD,TDEC_K, NDEC};
  //  TString decLab[NDEC];


  static TLorentzVector getTLV( EVENT::ReconstructedParticle* rp );
  static TLorentzVector getTLV( std::vector < EVENT::ReconstructedParticle* > rps );
  static TLorentzVector getTLV( EVENT::MCParticle* mc );
  static std::vector < EVENT::MCParticle* > getstablemctauDaughters( EVENT::MCParticle* mctau );
  static std::vector < EVENT::MCParticle* > getmctauDaughters( EVENT::MCParticle* mctau, bool onlyStable=false );
  static EVENT::MCParticle* getBestTrackMatch( EVENT::ReconstructedParticle* pfo, UTIL::LCRelationNavigator* relNavi );
  static EVENT::MCParticle* getBestCaloMatch ( EVENT::ReconstructedParticle* pfo, UTIL::LCRelationNavigator* relNavi ) ;

  static int getMCdecayMode( std::vector < EVENT::MCParticle* > mcps );
  static int getMCdecayMode( EVENT::MCParticle* mcp ) {
    std::vector < EVENT::MCParticle* > sd = getstablemctauDaughters( mcp );
    return getMCdecayMode ( sd );
  }

  static TString getTauDecLab( int code );

  static std::vector < EVENT::MCParticle* > findFinalTaus( EVENT::LCCollection* mccol );


};


struct rhoInfo {
  int trackIndex;
  int photonIndex[2];
  float mass_ggt;
  float mass_gg;
  float pi0_en;
  float pi0_prob;
  float pi0_imbalan;
  TLorentzVector track4mom;
  TLorentzVector pi0fitted4mom;  
};

class zInfo {
 public:
  zInfo() {}
  ~zInfo() {}

  void addTrack( const EVENT::ReconstructedParticle* rp ) { 
    assert ( rp->getCharge() != 0 );
    tracks.push_back(rp);
  }
  void addPhoton(  const EVENT::ReconstructedParticle* rp ) {
    assert ( rp->getCharge() == 0 );
    FSRcandidates.push_back( rp );
  }

  void reset() {tracks.clear(); FSRcandidates.clear();}

  float getMass() {
    TLorentzVector totmom(0,0,0,0);
    TLorentzVector amom(0,0,0,0);
    for (size_t i=0; i<tracks.size(); i++) {
      float mass(0);
      if      ( tauUtils::isElectron(tracks[i]) ) mass = tauUtils::m_el;
      else if ( tauUtils::isMuon(tracks[i]) )     mass = tauUtils::m_mu;
      else                                        mass = tauUtils::m_piChg;
      
      amom.SetVectM( TVector3 (tracks[i]->getMomentum()[0], tracks[i]->getMomentum()[1], tracks[i]->getMomentum()[2] ), mass );
      totmom+=amom;
    }
    for (size_t i=0; i<FSRcandidates.size(); i++) {
      float mass(0); // assume photons
      amom.SetVectM( TVector3 (FSRcandidates[i]->getMomentum()[0], FSRcandidates[i]->getMomentum()[1], FSRcandidates[i]->getMomentum()[2] ), mass );
      totmom+=amom;
    }
    return totmom.M();
  }

  const std::vector < const EVENT::ReconstructedParticle* > & getTracks() { return tracks; }
  const std::vector < const EVENT::ReconstructedParticle* > & getPhotons() { return FSRcandidates; }

 private:
  std::vector < const EVENT::ReconstructedParticle* > tracks;
  std::vector < const EVENT::ReconstructedParticle* > FSRcandidates;

};



#endif
