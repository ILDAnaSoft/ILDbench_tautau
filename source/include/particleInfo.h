#ifndef _PARTICLEINFO_H_
#define _PARTICLEINFO_H_

#include "TLorentzVector.h"
#include "TVector3.h"

#include "EVENT/Cluster.h"
#include "EVENT/Track.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"
#include "EVENT/TrackState.h"

#include <iostream>

class particleInfo{

 public:
  particleInfo(const EVENT::MCParticle* rp) { init(); set(rp); 
    //    std::cout << " mc particleInfo" << std::endl;
  }
  particleInfo(EVENT::Track* rp) { init(); set(rp); }
  particleInfo(const EVENT::ReconstructedParticle* rp) { init(); set(rp); }
  particleInfo( int pdg, TLorentzVector tl, int itype ) { init(); set(pdg, tl, itype); }
  ~particleInfo();

  void set( const EVENT::MCParticle* rp );
  void set( const EVENT::ReconstructedParticle* rp);
  void set( EVENT::Track* rp);
  void set( const EVENT::Cluster* cl);
  void set( int pdg, TLorentzVector tl, int itype);
  void set( TLorentzVector tl, int itype);
  TLorentzVector getFourMomentum(int itype);

  TVector3 getDvec(int itype ); // d vector is the vector from IP to 2d-PCA
  TVector3 getEvec(int itype ); // e vector is the vector from IP to 3d-PCA

  void Print();
  bool hasMC() {return _hasMc;}
  bool hasRECO() {return _hasReco;}
  int getPDG(int itype) {return _pdg[itype];}
  int getCharge( int imcrec );

  const EVENT::ReconstructedParticle* getRP() { return _rp; }

  void forceHasReco(bool b=true) {_hasReco=b;}

  const EVENT::TrackState* getTrackState() {return _trkSt;}
  EVENT::Track* getTrack() {return _trk;}

  bool isGoodRecoTrack();
  bool isGoodRecoPhoton();

  enum { imc=0, irec, irecIP, icalc, icalc2, ivtxcalc, ivtxcalc2, NN };

  void setIP(TVector3 pp) {
    _ip=pp;
    updateIP();
  }

  const EVENT::MCParticle* getMCP() {return _mcp;}

 private:

  void init();

  const EVENT::MCParticle* _mcp;
  const EVENT::ReconstructedParticle* _rp;
  const EVENT::Cluster* _cl;
  EVENT::Track* _trk;
  const EVENT::TrackState* _trkSt;

  bool _hasMc;
  bool _hasReco;
  bool _hasCalc;

  int _pdg[NN];
  TLorentzVector _4vec[NN];
  TVector3       _dVec[NN];
  TVector3       _eVec[NN];

  TVector3 _ip;

  void updateIP();

};

#endif
