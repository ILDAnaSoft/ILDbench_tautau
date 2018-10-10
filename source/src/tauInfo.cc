#include "tauInfo.h"
#include "tauUtils.h"

#include <iostream>
using std::cout;
using std::endl;

using namespace EVENT;

void tauInfo::addCharged      ( const ReconstructedParticle* rp ) {
  charged    .push_back(new particleInfo(rp));
  //  cout << "adding charged" << endl;
}

void tauInfo::addPizero       ( const ReconstructedParticle* rp ) {
  pizeros    .push_back(new particleInfo(rp));
  //  cout << "adding pi0" << endl;
}

void tauInfo::addPhoton       ( const ReconstructedParticle* rp ) {
  photons    .push_back(new particleInfo(rp));
  //  cout << "adding photon" << endl;
}


void tauInfo::addConversion   ( const ReconstructedParticle* rp ) {
  conversions.push_back(new particleInfo(rp));
  //  cout << "adding conversion" << endl;
}

void tauInfo::addCharged      ( const MCParticle* rp ) {
  charged.push_back(new particleInfo(rp));
  //  cout << "adding mc charged" << endl;
}

void tauInfo::addPhoton       ( const MCParticle* rp ) {
  photons.push_back(new particleInfo(rp));
  //  cout << "adding mc photon" << endl;
}

void tauInfo::addFsrPhoton    ( const MCParticle* rp ) {
  fsrphotons.push_back(new particleInfo(rp));
  //  cout << "adding mc fsr" << endl;
}

void tauInfo::addPizero       ( const MCParticle* rp ) {
  //  cout << "adding mc pi0" << endl;
  pizeros.push_back(new particleInfo(rp));
  for (size_t i=0; i<rp->getDaughters().size(); i++) {
    MCParticle* dd = rp->getDaughters()[i];
    if ( dd->getPDG()==22 ) {
      addPhoton( dd );
    } else {
      // cout << "warning, non-photon daughter of pi0 " << dd->getPDG() << endl;
    }
  }
}

void tauInfo::addNeutrino     ( const MCParticle* rp ) { neutrinos.push_back(new particleInfo(rp)); }
void tauInfo::addNeutralHadron( const MCParticle* rp ) {neutralHadrons.push_back(new particleInfo(rp)); }
void tauInfo::setTau          ( const MCParticle* rp ) { tau=new particleInfo(rp); }

std::vector <particleInfo*> tauInfo::getCharged () {return charged;}
std::vector <particleInfo*> tauInfo::getPhoton  () {return photons;}
std::vector <particleInfo*> tauInfo::getFsrPhoton  () {return fsrphotons;}
std::vector <particleInfo*> tauInfo::getPizero  () {return pizeros;}
std::vector <particleInfo*> tauInfo::getNeutrino() {return neutrinos;}
std::vector <particleInfo*> tauInfo::getNeutralHadron() {return neutralHadrons;}

TLorentzVector tauInfo::getChargedVisible4mom( int imeas ) {
  assert( imeas==particleInfo::imc || imeas==particleInfo::irec );
  TLorentzVector total4mom(0,0,0,0);
  for ( size_t i=0; i<getCharged().size(); i++)
    total4mom+=getCharged()[i]->getFourMomentum(imeas);
  return total4mom;
}


int tauInfo::getCharge( int imeas ) {
  assert( imeas==particleInfo::imc || imeas==particleInfo::irec );
  int charge(0);

  if ( imeas==particleInfo::imc ) {
    charge = getTau()->getMCP()->getCharge();
  } else {
    for ( size_t i=0; i<getCharged().size(); i++) {
      if ( imeas==particleInfo::irec && getCharged()[i]->hasRECO() ) charge+=getCharged()[i]->getRP()->getCharge();
      // if ( imeas==particleInfo::imc  && getCharged()[i]->hasMC()   ) charge+=getCharged()[i]->getMCP()->getCharge();
    }
  }

  return charge;
}

TLorentzVector tauInfo::getNeutralVisible4mom( int imeas ) {
  assert( imeas==particleInfo::imc || imeas==particleInfo::irec );
  TLorentzVector total4mom(0,0,0,0);

  //  if ( getPizero().size()>0 ) { // take info from pi zeros
  //for ( size_t i=0; i<getPizero().size(); i++)
  //  total4mom+=getPizero()[i]->getFourMomentum(imeas);
  //  } else { // take from photons
  //for ( size_t i=0; i<getPhoton().size(); i++) {
  //  total4mom+=getPhoton()[i]->getFourMomentum(imeas);
  // }
  //  }

  if ( imeas==particleInfo::imc ) {

    if ( getPizero().size()>0 ) { // take info from pi zeros
      for ( size_t i=0; i<getPizero().size(); i++)
        total4mom+=getPizero()[i]->getFourMomentum(imeas);
    } else { // take from photons
      for ( size_t i=0; i<getPhoton().size(); i++) {
        total4mom+=getPhoton()[i]->getFourMomentum(imeas);
      }
    }

  } else { // reco

    // pizeros and photons do not overlap
    //  cout << "hello from getNeutralVisible4mom nPi0, nGam = " << getPizero().size() << " " << getPhoton().size() << endl;

    for ( size_t i=0; i<getPizero().size(); i++) {
      //    cout << "pizero mass: " << getPizero()[i]->getFourMomentum(imeas).M() << " energy " << getPizero()[i]->getFourMomentum(imeas).E() << endl;
      //    getPizero()[i]->getFourMomentum(imeas).Print();
      total4mom+=getPizero()[i]->getFourMomentum(imeas);
    }

    //  cout << "sum pi0s:" ; total4mom.Print();
    //  cout << " mass " << total4mom.M() << endl;

    for ( size_t i=0; i<getPhoton().size(); i++) {
      //    cout << "photon mass: " << getPhoton()[i]->getFourMomentum(imeas).M() << endl;
      total4mom+=getPhoton()[i]->getFourMomentum(imeas);
    }

  }

  return total4mom;
}

TLorentzVector tauInfo::getVisible4mom( int imeas ) {
  return getChargedVisible4mom(imeas) + getNeutralVisible4mom(imeas);
}

vertexInfo* tauInfo::getDecayVertex( int imccalc ) {

  //  cout << "hello from getDecayVertex " << imccalc << " " << _decayVertex[imccalc] << endl;

  assert( imccalc>=0 && imccalc<particleInfo::NN);
  assert( imccalc==particleInfo::irec );

  if ( ! _decayVertex[imccalc] ) {
    vertexInfo* newvtx = new vertexInfo();
    for (size_t i=0; i<charged.size(); i++) {
      if ( charged[i]->hasRECO() && charged[i]->getTrack() ) {

        //      cout << "adding track " << i << " " << charged[i]->getTrack() << endl;

        newvtx->addTrack( charged[i]->getTrack() );
      }
    }
    _decayVertex[imccalc] = newvtx;
  }

  return _decayVertex[imccalc];
}

void tauInfo::Print() {
  cout << "printing tau info..." << endl;
  for (size_t i=0; i<charged.size(); i++) charged[i]->Print();
  for (size_t i=0; i<pizeros.size(); i++) pizeros[i]->Print();
  for (size_t i=0; i<photons.size(); i++) photons[i]->Print();
  for (size_t i=0; i<conversions.size(); i++) conversions[i]->Print();
  for (size_t i=0; i<neutrinos.size(); i++) neutrinos[i]->Print();
  tau->Print();
}

bool tauInfo::isGoodTau() {
  if ( charged.size()==0 ) return false;
  for ( size_t i=0; i<charged.size(); i++) {
    if ( ! charged[i]->isGoodRecoTrack() ) return false;
  }
  for (size_t i=0; i<pizeros.size(); i++)
    if ( ! pizeros[i]->hasRECO() ) return false;
  return true;
}

void tauInfo::setIP( TVector3 p ) {

  //  cout << "hello from tauInfo::setIP" << endl;

  for (size_t i=0; i<charged.size(); i++) {
    if ( charged[i]->hasRECO() ) charged[i]->setIP(p);
  }
  _ip = p;
  return;
}


TVector3 tauInfo::getMultiEvec(int irec) {

  //  cout << "hello from getMultiEvec " << irec << endl;

  assert( getCharged().size()>1 );
  assert( irec==particleInfo::irec );
  vertexInfo* thisvtx = getDecayVertex(irec);
  //  std::cout << "thisvtx " << thisvtx << std::endl;
  TVector3 lVector = thisvtx->getEigenVector(0);
  //  lVector.Print();
  //  thisvtx->getVertexPosition().Print();
  TVector3 eVector = lVector.Cross( thisvtx->getVertexPosition().Cross( lVector ) );
  //  eVector.Print();
  return eVector;
}


void tauInfo::calculateMcDecayLength() {
  particleInfo* ptau = getTau();
  float decayDist(0);
  if ( ptau->getMCP()->getDaughters().size()>0 ) {
    decayDist=0;
    for (int j=0; j<3; j++) decayDist += pow( ptau->getMCP()->getDaughters()[0]->getVertex()[j] , 2 );
    decayDist=sqrt(decayDist);
  } else {
    cout << "ERROR: MC tau with no daughters???" << endl;
    assert (0);
  }
  _decayLength[particleInfo::imc]=decayDist;
}

TVector3 tauInfo::getMCDecayVertex() {
  if ( _mcVertex.Mag()<=0 ) {
    if ( getNeutrino().size()>0 ) {
      _mcVertex.SetXYZ(
                       getNeutrino()[0]->getMCP()->getVertex()[0],
                       getNeutrino()[0]->getMCP()->getVertex()[1],
                       getNeutrino()[0]->getMCP()->getVertex()[2] );
    }
  }
  return _mcVertex;
}

void tauInfo::setDecayMode( int imc, int idd ) {
  _decayChannel[imc]=idd;
  return;
}

int tauInfo::getDecayMode( int imeas ) {

  int decMode(-1);

  if ( imeas==particleInfo::imc ) decMode = _decayChannel[imeas];

  else if ( imeas==particleInfo::irec ) {

    bool electron(false);
    bool muon(false);

    int nchrec(0);
    for (size_t i=0; i<charged.size(); i++) {
      if ( charged[i]->hasRECO() ) nchrec++;
      if ( tauUtils::isMuon( charged[i]->getRP() )==tauUtils::TIGHT ) muon=true;
      if ( tauUtils::isElectron( charged[i]->getRP() )==tauUtils::TIGHT ) electron=true;
    }
    int npiz(0);
    for (size_t i=0; i<pizeros.size(); i++) {
      if ( pizeros[i]->hasRECO() ) npiz++;
    }
    int npho(0);
    for (size_t i=0; i<photons.size(); i++) {
      if ( photons[i]->hasRECO() ) npho++;
    }

    //    cout << "number of hasRECO ch, piz, pho = " << nchrec << " " << npiz << " " << npho ;

    float vismass = getVisible4mom( particleInfo::irec ).M() ;
    //    cout << " ; visible mass = " << vismass << endl;

    if ( nchrec==1 ) {

      // if      ( muon          ) decMode = tauUtils::decayMu;
      // else if ( electron      ) decMode = tauUtils::decayEl;
      // else if ( vismass < 0.4 ) decMode = tauUtils::decayChPi;
      // else if ( vismass < 1.2 ) decMode = tauUtils::decayRho;
      // else decMode = tauUtils::decayA1_1p;

      // daniel updated - may10, 2017
      // use lepton id only for cases which look like single prong without pi0
      // in other cases, which look like rho, etc, assume prong is hadronic

      if ( vismass < 0.2 ) {
        if      ( muon          ) decMode = tauUtils::decayMu;
        else if ( electron      ) decMode = tauUtils::decayEl;
        else                      decMode = tauUtils::decayChPi;
      }

      else if ( vismass < 1.2 ) decMode = tauUtils::decayRho;
      else                      decMode = tauUtils::decayA1_1p;

      // keep this selection loose at this stage
      //      else if ( vismass > 0.5 && vismass < 1.0 ) decMode = tauUtils::decayRho; // daniel tightened mass req 29.05.2017
      //      else if ( vismass > 0.999 ) decMode = tauUtils::decayA1_1p;
      //      else                      decMode = tauUtils::decayUnrecognised; // below rho peak
      //    } else if ( nchrec==3 && !muon && !electron ) {

    } else if ( nchrec==3 ) { // daniel updated may2017: assume hadronic prongs
      decMode = tauUtils::decayA1_3p;
    } else {
      cout << "can't decide decay mode!" << endl;
    }

  }

  return decMode;

  // int mode(-1);
  // if ( getCharged().size()==1 ) {
  //   if ( getPizero().size()==0 ) mode=decayChPi;
  //   else if ( getPizero().size()==1 ) mode=decayRho;
  //   else if ( getPizero().size()==2 ) mode=decayA1_1p;
  // } else if ( getCharged().size()==3 ) {
  //   if ( getPizero().size()==0 ) mode=decayA1_3p;
  // }
  // return mode;
}

void tauInfo::cleanup() {

  for ( size_t i=0; i<getCharged().size(); i++) {
    //    std::cout << "deleting charged" << std::endl;
    delete getCharged()[i];
  }
  getCharged().clear();

  for ( size_t i=0; i<getPhoton().size(); i++) {
    //    cout << "deleting photon" << endl;
    delete getPhoton()[i];
  }
  getPhoton().clear();

  for ( size_t i=0; i<getFsrPhoton().size(); i++) {
    //    cout << "deleting fsr" << endl;
    delete getFsrPhoton()[i];
  }
  getFsrPhoton().clear();

  for ( size_t i=0; i<getPizero().size(); i++) {
    //    cout << "deleting pi0" << endl;
    delete getPizero()[i];
  }
  getPizero().clear();

  for ( size_t i=0; i<getNeutrino().size(); i++) {
    //    cout << "deleting neutrino" << endl;
    delete getNeutrino()[i];
  }
  getNeutrino().clear();

  for ( size_t i=0; i<getNeutralHadron().size(); i++) {
    //    cout << "deleting neuhad" << endl;
    delete getNeutralHadron()[i];
  }
  getNeutralHadron().clear();

  for (int i=0; i<particleInfo::NN; i++) {
    if ( _decayVertex[i] ) {
      delete _decayVertex[i];
      _decayVertex[i]=NULL;
    }
  }

  if ( tau ) delete tau;

}
