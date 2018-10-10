#include "particleInfo.h"
#include "tauUtils.h"

#include <iostream>
using std::cout;
using std::endl;

using namespace EVENT;


particleInfo::~particleInfo() {
  //  std::cout << "deleting particle info object" << std::endl;
}

void particleInfo::set( const MCParticle* rp ) {
  assert(rp);
  _hasMc=true;
  _mcp=rp;
  _pdg[imc] = rp->getPDG();
  _4vec[imc] = tauUtils::getFourMomentum(_mcp);
  _dVec[imc] = tauUtils::getDvector(_mcp);

  _eVec[imc] = _4vec[imc].Vect().Cross(_dVec[imc].Cross(_4vec[imc].Vect()) );
  // now make sure that size of eVec corresponds to distance from IP to real 3d PCA
  double theta = _dVec[imc].Angle(_eVec[imc]); // angle between e and d
  double modE = _dVec[imc].Mag()*cos(theta);   // dist from ip to 3d PCA
  _eVec[imc]*=modE/_eVec[imc].Mag();
}

void particleInfo::set( const ReconstructedParticle* rp) {
  assert(rp);
  _hasReco=true;
  _rp=rp;

  if ( fabs( _rp->getCharge() ) > 0.5 ) { // treat as single charged particle

    if ( _rp->getTracks().size()==0 ) {
      cout << "ERROR charged RP with no track!" << endl;
    } else if ( _rp->getTracks().size()==1 ) {
      _trk = _rp->getTracks()[0];
    } else {

      // cout << "input momentum " << _rp->getMomentum()[0] << " " << _rp->getMomentum()[1] << " " << _rp->getMomentum()[2] << endl;
      // cout << "more than one track for charged reco part: " << _rp->getTracks().size() << endl;
      // for (size_t kk=0; kk<_rp->getTracks().size(); kk++) {
      // 	cout << "track " << kk << endl;
      // 	set(_rp->getTracks()[kk]);
      // 	_4vec[irec].Print();
      // }

      // take one with hit at smaller radius
      float minRad=9999;
      EVENT::Track* chosenTrack(0);
      for (size_t it=0; it<_rp->getTracks().size(); it++) {
	EVENT::Track* trk = _rp->getTracks()[it];
	//      cout << "track " << it << "inner radius " << trk-> getRadiusOfInnermostHit () << endl;
	if ( trk-> getRadiusOfInnermostHit () < minRad ) {
	  minRad = trk-> getRadiusOfInnermostHit ();
	  chosenTrack=trk;
	}
      }
      assert(chosenTrack);
      _trk = chosenTrack;
    }

    _pdg[irec] = int(_rp->getCharge()) * 211; // assume pi+-
    assert (_trk);
    set( _trk );

  } else { // neutral

    if ( _rp->getClusters().size()>=0 || _rp->getParticles().size()>0 ) {
      
      _4vec[irec].SetXYZT( _rp->getMomentum()[0], _rp->getMomentum()[1], _rp->getMomentum()[2], _rp->getEnergy() ); //    = tauUtils::getFourMomentum(_cl);

      if ( _rp->getType()==22 || _4vec[irec].M() < 0.1 ) { // assume photon, adjust energy
	_pdg[irec] = 22;
	double pmag = _4vec[irec].Vect().Mag();
	double e = _4vec[irec].E();
	double scale = e/pmag;
	_4vec[irec].SetVectM( scale*_4vec[irec].Vect(), 0.0 );
      } else { // call it a pizero
	_pdg[irec] = 111;
      }	
      _dVec[irec].SetXYZ(0,0,0);
      _eVec[irec].SetXYZ(0,0,0);
    } else {
      cout << "WEIRD, neutral reco particle with " <<  _rp->getClusters().size() << " clusters" << endl;
      assert(0);
    }
  } 

  return;
}

void particleInfo::set( EVENT::Track* rp) {
  assert(rp);
  _hasReco=true;
  _trk = rp;
  assert (_trk);
  _trkSt = _trk->getTrackState (EVENT::TrackState::AtIP) ;
  assert (_trkSt);
  _4vec[irec] = tauUtils::getFourMomentum(_trkSt, 0.14);
  _dVec[irec] = tauUtils::getDvector( _trkSt );
  _eVec[irec] = _4vec[irec].Vect().Cross(_dVec[irec].Cross(_4vec[irec].Vect()) );

  double theta = _dVec[irec].Angle(_eVec[irec]); // angle between e and d
  double modE = _dVec[irec].Mag()*cos(theta);   // dist from ip to 3d PCA
  _eVec[irec]*=modE/_eVec[irec].Mag();

  int sign = rp->getOmega()>0 ? 1 : -1;
  _pdg[irec] = sign * 211;

  updateIP();
}


void particleInfo::set( const Cluster* cl ) {

  _cl = cl;
  assert (_cl);

  _hasReco=true;

  float mag(0);
  for (int i=0; i<3; i++) mag+=pow(_cl->getPosition()[i],2);
  mag=sqrt(mag);
  float norm=_cl->getEnergy()/mag;

  _4vec[irec].SetXYZT( norm*_cl->getPosition()[0], norm*_cl->getPosition()[1], norm*_cl->getPosition()[2], _cl->getEnergy() );
  _dVec[irec].SetXYZ(0,0,0);
  _eVec[irec].SetXYZ(0,0,0);

  return;
}



void particleInfo::set( int pdg, TLorentzVector tl, int itype) {
  // usually for neutrino...
  assert ( itype<NN );
  _4vec[itype] = tl;
  _pdg[itype] = pdg;
  return;
}


void particleInfo::set( TLorentzVector tl, int itype) {
  assert ( itype<NN );
  _4vec[itype] = tl;
  return;
}


int particleInfo::getCharge( int imcrec ) {
  switch ( abs( _pdg[imcrec] ) ) {
  case 22:
  case 221: // eta
  case 2112:
  case 111:
  case 130:
  case 310:
  case 12:
  case 14:
  case 16:
    return 0;
  case 211:
  case 321:
  case 323:
  case 2212:
  case 14122: // need to add this one
    return _pdg[imcrec]<0 ? -1 : 1;
  case 11:
  case 13:
  case 15:
    return _pdg[imcrec]<0 ? 1 : -1;
  case 1:
  case 3:
  case 5:
    return _pdg[imcrec]<0 ? 1 : -1; // ahh, we have it as integer!!
  case 2:
  case 4:
  case 6:
    return _pdg[imcrec]<0 ? 1 : -1;
  default:
    cout << " particleInfo::getCharge doesn't know charge of " << _pdg[imcrec] << endl;
    assert(0);
  }
  return 0;
}

TLorentzVector particleInfo::getFourMomentum(int itype) {
  assert(itype<NN);
  return _4vec[itype];
}

TVector3 particleInfo::getDvec(int itype) {
  assert(itype<NN);
  return _dVec[itype];
}

void particleInfo::updateIP() {

  if ( !_trkSt ) {
    cout << "ERROR asking to update IP of particle with no trackState..." << endl;
    cout << "MCParticle " <<  _mcp << endl;
    cout << "ReconstructedParticle* " << _rp << endl;
    cout << "Cluster* " << _cl << endl;
    cout << "Track* " << _trk << endl;
    cout << "TrackState* " << _trkSt << endl;
    cout << "PDG " << _pdg[imc] << endl;
    cout << _rp->getCharge( ) << endl;
    cout << _rp->getTracks().size() << " " << _rp->getClusters().size() << endl;
  }

  assert (_trkSt);

  //  cout << "hello from updateIP" << endl;

  TVector3 momdir;
  tauUtils::calculateTrackWRTip( _trkSt, momdir, _dVec[irecIP], _eVec[irecIP] , _ip );
  _4vec[irecIP].SetVectM( (_4vec[irec].Vect().Mag())*momdir , _4vec[irec].M() );

  //  cout << "dvecs:  angle = " << _dVec[irec].Angle( _dVec[irecIP] ) << endl;
  //  _dVec[irec].Print();
  //  _dVec[irecIP].Print();
  //  cout << "4vecs:  angle = " << _4vec[irec].Vect().Angle( _4vec[irecIP].Vect() ) << endl;
  //  _4vec[irec].Print();
  //  _4vec[irecIP].Print();
  //  cout << "d-4vec angle = " << _4vec[irec].Vect().Angle( _dVec[irec] ) << " " << _4vec[irecIP].Vect().Angle( _dVec[irecIP] ) << endl;

  return;
}

TVector3 particleInfo::getEvec(int itype) {
  assert(itype<NN);
  return _eVec[itype];
}

void particleInfo::Print() {
  cout << "  printing particleInfo: ";
  //    cout << "  mc, rp, cl, trk, trkSt = " << _mcp << " " << _rp << " " << _cl << " " << _trk << " " << _trkSt  << endl;
  cout << "  mcPDG " << _pdg[0] << ", recoPDG " << _pdg[1] << " hasMC " << _hasMc << " hasRECO " << _hasReco << endl;
  // for (int in=0; in<NN; in++) {

  cout << "  rec    " << _4vec[irec].E() << " " << _4vec[irec].X() << " " << _4vec[irec].Y() << " " << _4vec[irec].Z() << " : " << _4vec[irec].M() << endl;
  cout << "  calc   " << _4vec[icalc].E() << " " << _4vec[icalc].X() << " " << _4vec[icalc].Y() << " " << _4vec[icalc].Z() << " : " << _4vec[icalc].M() << endl;
  cout << "  calc2  " << _4vec[icalc2].E() << " " << _4vec[icalc2].X() << " " << _4vec[icalc2].Y() << " " << _4vec[icalc2].Z() << " : " << _4vec[icalc2].M() << endl;
  cout << "  mc     " << _4vec[imc].E() << " " << _4vec[imc].X() << " " << _4vec[imc].Y() << " " << _4vec[imc].Z() << " : " << _4vec[imc].M() << endl;

  //   _dVec[in].Print();
  //   _eVec[in].Print();
  // }
  return;
}

void particleInfo::init() {
  _mcp=0;
  _rp=0;
  _cl=0;
  _trk=0;
  _trkSt=0;
  _hasMc=false;
  _hasReco=false;
  _pdg[0]=0; _pdg[1]=0;

  for (int n=0; n<NN; n++) {
    _4vec[n].SetXYZT(0,0,0,0);
    _dVec[n].SetXYZ(0,0,0);
    _eVec[n].SetXYZ(0,0,0);
  }

  //  std::cout << "initialising particleInfo object" << std::endl;

}

bool particleInfo::isGoodRecoTrack() {
  return hasRECO();
}

bool particleInfo::isGoodRecoPhoton() {
  return hasRECO();
}
