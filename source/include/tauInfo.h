#ifndef _TAUINFO_H_
#define _TAUINFO_H_

#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"

#include "particleInfo.h"
#include "vertexInfo.h"

#include <iostream>

// this class is to store MC and measured 4momenta and vertices of tau and its products

class tauInfo{

 public:

  tauInfo() {

    //    std::cout << "constructing tauinfo" << std::endl;

    _ip.SetXYZ(0,0,0);
    _mcVertex.SetXYZ(0,0,0);
    for (int i=0; i<particleInfo::NN; i++) {
      _decayVertex[i]=NULL;
      _decayChannel[i]=-1;
    }
    _multiP_centralMom=-999;
    tau=0;

    for (int i=0; i<4; i++) {
      _totExtraEnInCone[i]=999;
      _chgExtraEnInCone[i]=999;
    }

  }

  ~tauInfo() {

    //    std::cout << "destructing tauinfo" << std::endl;

    cleanup();
  }

  void addCharged      ( const EVENT::ReconstructedParticle* rp );
  void addPhoton       ( const EVENT::ReconstructedParticle* rp );
  void addPizero       ( const EVENT::ReconstructedParticle* rp );
  void addConversion   ( const EVENT::ReconstructedParticle* rp );
  

  void addCharged      ( const EVENT::MCParticle* rp );
  void addPhoton       ( const EVENT::MCParticle* rp );
  void addPizero       ( const EVENT::MCParticle* rp );
  void addNeutrino     ( const EVENT::MCParticle* rp );
  void addNeutralHadron( const EVENT::MCParticle* rp );
  void addFsrPhoton    ( const EVENT::MCParticle* rp );

  void setTau          ( const EVENT::MCParticle* rp );

  void setConeTotEn( int i, float x ) {
    assert ( i>=0 && i<4 );
    _totExtraEnInCone[i]=x;
  }
  void setConeChgEn( int i, float x ) {
    assert ( i>=0 && i<4 );
    _chgExtraEnInCone[i]=x;
  }

  //  bool setPsi ( double psi , bool verbose=false );

  std::vector <particleInfo*> getCharged();
  std::vector <particleInfo*> getPhoton ();
  std::vector <particleInfo*> getFsrPhoton ();
  std::vector <particleInfo*> getPizero ();
  std::vector <particleInfo*> getNeutrino();
  std::vector <particleInfo*> getNeutralHadron();
  particleInfo* getTau() {return tau;}

  int getCharge( int imeas );

  TLorentzVector getVisible4mom( int imeas );
  TLorentzVector getChargedVisible4mom( int imeas );
  TLorentzVector getNeutralVisible4mom( int imeas );

  vertexInfo* getDecayVertex( int imccalc );
  TVector3 getMCDecayVertex();

  void setIP( TVector3 p ); // {_ip=p;}

  float get_multiP_centralMom() {return _multiP_centralMom;}

  void Print();

  bool isGoodTau();

  void setDecayMode( int imc, int idd );
  int getDecayMode( int imeas );

  TVector3 getMultiEvec(int irec);

  float getConeTotEn( int i ) {
    assert ( i>=0 && i<4 );
    return _totExtraEnInCone[i];
  }
  float getConeChgEn( int i ) {
    assert ( i>=0 && i<4 );
    return _chgExtraEnInCone[i];
  }


 private:

  std::vector <particleInfo*> charged;
  std::vector <particleInfo*> pizeros;
  std::vector <particleInfo*> photons;
  std::vector <particleInfo*> fsrphotons;
  std::vector <particleInfo*> conversions;
  std::vector <particleInfo*> neutrinos;
  std::vector <particleInfo*> neutralHadrons;

  particleInfo* tau;

  // the calculated/fitted parameters of the neutrino
  double    _psi; // the neutrino angle within the plane

  int _decayChannel[particleInfo::NN];

  void calculateMcDecayLength();

  TVector3 _decayVtxPos[particleInfo::NN];
  double _decayLength[particleInfo::NN];
  double _decayTime[particleInfo::NN];
  vertexInfo* _decayVertex[particleInfo::NN];
  TVector3 _mcVertex;
  float _multiP_centralMom;
  TVector3 _ip;

  float _totExtraEnInCone[4];
  float _chgExtraEnInCone[4];

  void cleanup();
};

#endif
