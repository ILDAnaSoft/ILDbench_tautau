#ifndef _MINMALTAUTAU_H_
#define _MINMALTAUTAU_H_

#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TF1.h"

#include <cassert>

class tautau2fEventFitter {

 public:

  tautau2fEventFitter() {init(); functor=0;minimizer=0;}
  ~tautau2fEventFitter() {
    if (functor) {delete functor; functor=0;}
    if (minimizer) {delete minimizer; minimizer=0;}
  }

  void reset() {init();}
  void setEcom(float ecom) {_ecom=ecom;}

  // event quantitiy setters
  void setRecoil(TLorentzVector p) {recoil=p; prepared=false;}

  void setIP(TVector3 p) {_ip = p; prepared=false;}

  // tau properties
  //  void setVertex(int itau, TVector3 pos) {assert( itau>=0 && itau<ntau ); multiProng[itau]=true; tauVertices[itau]=pos;}

  void setNeutralHadronMomentum(int itau, TLorentzVector p) {assert( itau>=0 && itau<ntau ); neutralHadronMomentum[itau]=p; prepared=false;}
  void setChargedHadronMomentum(int itau, TLorentzVector p) {assert( itau>=0 && itau<ntau ); chargedHadronMomentum[itau]=p; prepared=false;}

  // this is some point on the trajectory, wrt some nominal IP
  void setChargedHadronImpactParameterVector(int itau, TVector3 d) {assert( itau>=0 && itau<ntau ); chargedHadronImpactParameterVector[itau]=d; prepared=false;}

  void setVertex(int itau, TVector3 pos) {assert( itau>=0 && itau<ntau ); multiProng[itau]=true; tauVertices[itau]=pos;}

  bool fitIt_single_single(int bestStrategy);
  //  void fitIt_multi_multi();
  //  void fitIt_multi_multi_noIPZ();
  
  bool setIPz( double z );
  bool setPsi( int itau, double psi );
  double getPt_plusISR(bool j0, bool j1);
  double getPt(bool j0, bool j1);
  double getPt();

  double getDecayLength(int itau, bool jsol);
  double getDecayLength(int itau);

  double getDecayTimePDF(int itau, bool jsol);
  double getDecayTimePDF(int itau);

  TLorentzVector getNeutrino4Momentum(int itau, bool jsol);
  TLorentzVector getNeutrino4Momentum(int itau);

  TLorentzVector getTau4Momentum(int itau, bool jsol);
  TLorentzVector getTau4Momentum(int itau);

  TLorentzVector getEvent4Momentum(bool j0, bool j1);
  TLorentzVector getEvent4Momentum();

  TLorentzVector getEventVisible4Momentum();

  void setVerbose(bool verb=true) {verbose=verb;}

  void calculate_multi_pdf(float ipz);
  float get_multiNoIP_lifePdf(int itau, int jsol){return multi_noip_life_pdf[itau][jsol];}
  float get_multiNoIP_nuenPdf(int itau, int jsol){return multi_noip_nuen_pdf[itau][jsol];}

  TF1* f_a1_nuE;
  float getA1_Enu_pdf(float nuen);

  double getPsi(int itau) { assert(itau>=0 && itau<ntau); return tauPsi[itau]; }

  bool goodFit() {return goodfit;}

 private:

  float _ecom{};

  TLorentzVector recoil{};
  enum { ntau=2 };

  bool multiProng[ntau];
  TVector3 tauVertices[ntau];

  TLorentzVector neutralHadronMomentum[ntau];
  TLorentzVector chargedHadronMomentum[ntau];
  TVector3 chargedHadronImpactParameterVector[ntau];
  TVector3 chargedHadronEVector[ntau];

  TVector3 normalToTrackPlane[ntau];
  TVector3 kPerp[ntau];
  TLorentzVector hadMom[ntau];

  TVector3 _ip;
  TVector3 hPlane[ntau];
  TVector3 fPlane[ntau];

  double tauPsi[ntau];
  bool isolution[ntau];

  bool  multi_noip_valid_sol[ntau][2];
  float multi_noip_life_pdf[ntau][2];
  float multi_noip_nuen_pdf[ntau][2];

  TLorentzVector neutrinoMomentum[ntau][2]; // two solutions per psi

  //  TF1* f_a1_nuE;
  //  float getA1_Enu_pdf(float nuen);

  void init();

  void prepareForFit();
  bool calculateNeutrino(int itau);
  double fitfn_single_single(const double* pars);
  bool verbose{};
  bool prepared{};
  bool goodfit{};

  TVector3 getEVector( TVector3 momentum, TVector3 impactVector );


  static const float m_tau;
  static const float ctau_tau;

  ROOT::Math::Functor* functor{}; // ( this, & tautau2fEventFitter::fitfn_single_single, 2);
  ROOT::Math::Minimizer* minimizer{};

  
};



#endif
