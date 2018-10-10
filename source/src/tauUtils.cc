#include "tauUtils.h"

#include "tauInfo.h"
//#include "MySimpleHelix.h"
#include "SimpleHelix.h"
#include <CalorimeterHitType.h>

using std::cout;
using std::endl;

using namespace EVENT;
using namespace UTIL;

const double tauUtils::crossingAngle = 0.007;

const double tauUtils::omega_pt_conv = 2.99792458e-4; // for mm, Tesla, GeV
const double tauUtils::bfield = 3.5; // Tesla

const double tauUtils::ecal_stochastic = 0.17;
const double tauUtils::ecal_constant = 0.01;

const double tauUtils::m_tau = 1.777;
const double tauUtils::ctau_tau = 87.03e-3;
const double tauUtils::m_z     = 91.19;
const double tauUtils::m_piChg = 0.139568;
const double tauUtils::m_pi0 = 0.135;
const double tauUtils::m_mu    = 0.105;
const double tauUtils::m_el    = 0.000511;
const double tauUtils::m_rho   = 0.775;
const double tauUtils::g_rho   = 0.150;
const double tauUtils::Pi     = acos(-1);
const double tauUtils::twoPi  = 2*Pi;
const double tauUtils::halfPi = 0.5*Pi;

TLorentzVector tauUtils::getFourMomentum( const Cluster* cl ) {
  double en = cl->getEnergy();
  double pos[3];
  double mag(0);
  for (int i=0; i<3; i++) {
    pos[i] = cl->getPosition()[i];
    mag+=pow(pos[i],2);
  }
  mag=sqrt(mag);
  return TLorentzVector( en*pos[0]/mag, en*pos[1]/mag, en*pos[2]/mag, en );
}

TLorentzVector tauUtils::getFourMomentum( const TrackState* tst , double mass ) {
  assert(tst);
  double pt = omega_pt_conv*fabs(bfield/tst->getOmega());
  double phi0 = tst->getPhi();
  double tanl = tst->getTanLambda();
  TLorentzVector fv;
  fv.SetXYZM( pt*cos(phi0), pt*sin(phi0), pt*tanl, mass );
  return fv;
}

TLorentzVector tauUtils::getFourMomentum( const MCParticle* rp ) {
  assert(rp);
  return TLorentzVector( rp->getMomentum()[0],
                         rp->getMomentum()[1],
                         rp->getMomentum()[2],
                         rp->getEnergy() );
}

// from IP to 2d PCA
// this is the closest point to the ref point in 2d
TVector3 tauUtils::getDvector( const TrackState* tst ) {
  assert(tst);
  return TVector3 ( tst->getReferencePoint() [0] - tst->getD0()*sin(tst->getPhi()),
                    tst->getReferencePoint() [1] + tst->getD0()*cos(tst->getPhi()),
                    // tst->getReferencePoint() [1] + tst->getZ0() );
                    tst->getReferencePoint() [2] + tst->getZ0() ); // bug fix djeans 26June
}


void tauUtils::calculateTrackWRTip( const TrackState* tst, TVector3& momentumDir, TVector3& dVec, TVector3& eVec , TVector3 IP ) {

  assert(tst);

  LCVector3D refpt(  tst->getReferencePoint()[0],  tst->getReferencePoint()[1],  tst->getReferencePoint()[2] );
  //  MySimpleHelix shel( tst->getD0(),
  SimpleHelix shel( tst->getD0(),
                    tst->getPhi(),
                    tst->getOmega(),
                    tst->getZ0(),
                    tst->getTanLambda(),
                    refpt
		    );

  //  cout << "New IP " << IP.X() << " " << IP.Y() << " " << IP.Z() << endl;

  LCVector3D newrefpt( IP.X(), IP.Y(), IP.Z() );
  double new_s = shel.getPathAt( newrefpt );

  LCVector3D new_pca = shel.getPosition( new_s );
  LCVector3D new_direction = shel.getDirection ( new_s );

  //  cout << "new PCA " << new_pca.getX() << " " << new_pca.getY() << " " << new_pca.getZ() << endl;
  //  cout << "new direction " << new_direction.getX() << " " << new_direction.getY() << " " << new_direction.getZ() << endl;

  momentumDir.SetXYZ( new_direction.getX(), new_direction.getY(), new_direction.getZ() );
  dVec.SetXYZ( new_pca.getX(), new_pca.getY(), new_pca.getZ() );

  eVec = momentumDir.Cross( dVec.Cross( momentumDir ) );

  // now make sure that size of eVec corresponds to distance from IP to real 3d PCA
  double theta = dVec.Angle(eVec); // angle between e and d
  double modE = dVec.Mag()*cos(theta);   // dist from ip to 3d PCA
  eVec*=modE/eVec.Mag();

  return;
}



TVector3 tauUtils::getDvector( const MCParticle* mc ) {
  assert (mc);
  double vtx[3] = { mc->getVertex()[0], mc->getVertex()[1], mc->getVertex()[2] };
  double mom[3] = { mc->getMomentum()[0], mc->getMomentum()[1], mc->getMomentum()[2] };
  // d = vtx + alpha*mom
  double alpha = -(vtx[0]*mom[0]+vtx[1]*mom[1])/(mom[0]*mom[0]+mom[1]*mom[1]); // closest approach in 2d
  return TVector3( vtx[0]+alpha*mom[0], vtx[1]+alpha*mom[1], vtx[2]+alpha*mom[2] );
}


TVector3 tauUtils::getPolarimeter_rho ( TLorentzVector charged, TLorentzVector neutral, TLorentzVector neutrino ) {
  // calculate polarimeter vector in tau RF for rho decays
  TLorentzVector tvtau = charged + neutral + neutrino;
  charged.Boost( - tvtau.BoostVector() );
  neutral.Boost( - tvtau.BoostVector() );
  neutrino.Boost( - tvtau.BoostVector() );
  TLorentzVector q = charged-neutral;
  TLorentzVector Q = charged+neutral;
  TVector3 polarimeterVector = 2*(q.Dot(neutrino))*q.Vect() - q.M2()*neutrino.Vect();
  // for hadronic states, polarimeter size = 1
  polarimeterVector *= -1./polarimeterVector.Mag();
 return polarimeterVector;
}

TVector3 tauUtils::getPolarimeter_pi ( TLorentzVector charged, TLorentzVector neutrino ) {
  // calculate polarimeter vector in tau RF for pi decays
  TLorentzVector tvtau = charged + neutrino;
  charged.Boost( - tvtau.BoostVector() );
  neutrino.Boost( - tvtau.BoostVector() );
  TVector3 polarimeterVector = charged.Vect();
  polarimeterVector *= -1./polarimeterVector.Mag();
  return polarimeterVector;
}

TLorentzVector tauUtils::getPolarimeterTLV (  TLorentzVector charged, TLorentzVector neutral, TLorentzVector neutrino ) {
  TLorentzVector pol;
  if ( neutral.E() > 0 ) { // if there is any neutral energy
    pol = getPolarimeterTLV_rho ( charged, neutral, neutrino );
  } else {
    pol = getPolarimeterTLV_pi ( charged, neutrino );
  }
  return pol;
}

TLorentzVector tauUtils::getPolarimeterTLV (  TLorentzVector charged, TLorentzVector neutral, TLorentzVector neutrino, int decayMode ) {
  TLorentzVector pol;
  if ( decayMode==decayRho ) {
    pol = getPolarimeterTLV_rho ( charged, neutral, neutrino );
  } else if ( decayMode==decayChPi ) {
    pol = getPolarimeterTLV_pi ( charged, neutrino );
  } else {
    cout << "tauUtils::getPolarimeterTLV WARNING, do not know how to calculate polarimeter for decay mode " << decayMode << endl;
  }
  return pol;
}


TLorentzVector tauUtils::getPolarimeterTLV_rho ( TLorentzVector charged, TLorentzVector neutral, TLorentzVector neutrino ) {
  // calculate polarimeter vector in given frame for rho decays
  //  TLorentzVector tvtau = charged + neutral + neutrino;
  // charged.Boost( - tvtau.BoostVector() );
  // neutral.Boost( - tvtau.BoostVector() );
  // neutrino.Boost( - tvtau.BoostVector() );
  TLorentzVector q = charged-neutral;
  TLorentzVector Q = charged+neutral;
  TLorentzVector polarimeterVector = 2*(q.Dot(neutrino))*q - q.M2()*neutrino;
  // for hadronic states, polarimeter size = 1
  //  polarimeterVector *= -1./polarimeterVector.Mag();
  return polarimeterVector;
}

TLorentzVector tauUtils::getPolarimeterTLV_pi ( TLorentzVector charged, TLorentzVector neutrino ) {
  // calculate polarimeter vector in tau RF for pi decays
  // TLorentzVector tvtau = charged + neutrino;
  // charged.Boost( - tvtau.BoostVector() );
  // neutrino.Boost( - tvtau.BoostVector() );
  // TLorentzVector polarimeterVector = -charged;
  //  polarimeterVector *= -1./polarimeterVector.Mag();
  return charged;
}

float tauUtils::getPhiStar ( TLorentzVector tau0, TLorentzVector tau1, TLorentzVector pol0, TLorentzVector pol1 ) {
  // boost pols to the tau-tau RF
  TVector3 boost = -(tau0+tau1).BoostVector();

  TLorentzVector tt0(tau0); tt0.Boost(boost);
  TVector3 refDir = tt0.Vect(); refDir=refDir.Unit();

  pol0.Boost(boost);
  pol1.Boost(boost);

  TVector3 pp0 = pol0.Vect();
  TVector3 pp1 = pol1.Vect();

  pp0 -= ( pp0.Dot( refDir ) )*refDir;
  pp1 -= ( pp1.Dot( refDir ) )*refDir;

  float phi = pp0.Angle( pp1 );
  if ( ( pp0.Cross(pp1) ).Dot( refDir ) < 0 )  phi*=-1;

  return phi;
}

float tauUtils::getPhiStar ( TLorentzVector tau0, TVector3 pol0, TVector3 pol1 ) {

  TVector3 refDir = tau0.Vect(); refDir*=1./refDir.Mag(); // tau0 momentum direction in lab

  pol0 -= (pol0.Dot(refDir))*refDir;
  pol1 -= (pol1.Dot(refDir))*refDir;

  float phi = pol0.Angle(pol1);
  if ( ( pol0.Cross(pol1) ).Dot(refDir) < 0 ) phi*=-1;

  return phi;
}

float tauUtils::getPhiStarRhoRho ( TLorentzVector charged0, TLorentzVector neutral0, TLorentzVector neutrino0,
				   TLorentzVector charged1, TLorentzVector neutral1, TLorentzVector neutrino1 ) {
  TLorentzVector tau0 = charged0+neutral0+neutrino0;
  TLorentzVector tau1 = charged1+neutral1+neutrino1;


  TLorentzVector pol0 = neutral0.E()>0 ? getPolarimeterTLV_rho( charged0, neutral0, neutrino0 ): getPolarimeterTLV_pi(charged0, neutrino0);
  TLorentzVector pol1 = neutral1.E()>0 ? getPolarimeterTLV_rho( charged1, neutral1, neutrino1 ): getPolarimeterTLV_pi(charged1, neutrino1);
  //  TLorentzVector pol1 = getPolarimeterTLV_rho( charged1, neutral1, neutrino1 );

  return getPhiStar( tau0, tau1, pol0, pol1 );
}

void tauUtils::rouge_getThetaDphi( TLorentzVector charged0, TLorentzVector neutral0, TLorentzVector neutrino0,
				   TLorentzVector charged1, TLorentzVector neutral1, TLorentzVector neutrino1,
				   float & theta0, float& theta1, float& dPhi ) {

  // tau 4-momenta
  TLorentzVector tau0 = charged0+neutral0+neutrino0;
  TLorentzVector tau1 = charged1+neutral1+neutrino1;

  // polarimeters in lab frame
  TLorentzVector pol0 = neutral0.E()>0 ? getPolarimeterTLV_rho( charged0, neutral0, neutrino0 ) : getPolarimeterTLV_pi ( charged0, neutrino0 );
  TLorentzVector pol1 = neutral1.E()>0 ? getPolarimeterTLV_rho( charged1, neutral1, neutrino1 ) : getPolarimeterTLV_pi ( charged1, neutrino1 );

  // reference directions for polar angle of polarisations
  TLorentzVector z1(tau0); z1.Boost( -tau1.BoostVector() ); // direction of flight of tau0 in tau1 RF
  TLorentzVector z0(tau1); z0.Boost( -tau0.BoostVector() ); // vice versa

  TLorentzVector p0(pol0); p0.Boost( -tau0.BoostVector() ); // polarimeter of tau0 in its own rest frame
  TLorentzVector p1(pol1); p1.Boost( -tau1.BoostVector() ); // polarimeter of tau1 in its own rest frame

  theta0 = p0.Vect().Angle( z0.Vect() );
  theta1 = p1.Vect().Angle( -z1.Vect() );

  // for deltaPhi, boost into "higgs" frame
  TLorentzVector tt=tau0+tau1;
  tau0.Boost( -tt.BoostVector() );
  TVector3 refdir = tau0.Vect(); refdir=refdir.Unit();
  p0=pol0; p0.Boost( -tt.BoostVector() );
  p1=pol1; p1.Boost( -tt.BoostVector() );

  // in H frame, components of polarisation perp to reference direction
  TVector3 norm0 = p0.Vect() - (p0.Vect().Dot(refdir))*refdir;
  TVector3 norm1 = p1.Vect() - (p1.Vect().Dot(refdir))*refdir;

  // signed angle betweem these
  dPhi = norm0.Angle(norm1);
  if ( (norm0.Cross(norm1)).Dot(refdir) < 0 ) dPhi = -dPhi;

  return;
}

void tauUtils::gunion_getThetaDphi( TLorentzVector charged0, TLorentzVector neutral0, TLorentzVector neutrino0,
				    TLorentzVector charged1, TLorentzVector neutral1, TLorentzVector neutrino1,
				    float & theta0, float& theta1, float& dPhi ) {

  // tau 4-momenta
  //  TLorentzVector tau0 = charged0+neutral0+neutrino0;
  //  TLorentzVector tau1 = charged1+neutral1+neutrino1;

  TLorentzVector tlchg[2], tlpi0[2], tlneu[2];
  tlchg[0]=charged0;
  tlpi0[0]=neutral0;
  tlneu[0]=neutrino0;
  tlchg[1]=charged1;
  tlpi0[1]=neutral1;
  tlneu[1]=neutrino1;

  TLorentzVector tau[2];
  for (int i=0; i<2; i++) {
    tau[i] = tlchg[i]+tlpi0[i]+tlneu[i];
  }

  // z axis is direction of tau0 in tau-tau rest frame
  // TLorentzVector tlax = tau[0];
  // tlax.Boost( - (tau[0]+tau[1]).BoostVector() );
  // TVector3 axis = tlax.Vect();
  // axis*=1./axis.Mag();

  TVector3 axis = getGunionAxis( tau[0], tau[1] );

  // polarimeters defined in respective tau rest frames
  TVector3 S[2]; // this is the 3-d polarimeter
  TVector3 Sn[2]; // this is the component normal to the tau flight direction

  for (int i=0; i<2; i++) {
    //// boost to tau's rest frame
    //tlchg[i].Boost( - tau[i].BoostVector() );
    //tlneu[i].Boost( - tau[i].BoostVector() );
    //tlpi0[i].Boost( - tau[i].BoostVector() );

    if ( tlpi0[i].E()>0 ) { // treat as rho
      // m_tau * ( tlchg[i].Vect() - tlpi0[i].Vect() )*( tlchg[i].E() - tlpi0[i].E() ) + tlneu[i].Vect() * ( (tlchg[i]+tlpi0[i]).M2()/2.);
      S[i] = getGunionPolarimeter_rho(tlchg[i], tlpi0[i], tlneu[i]);
    } else {
      //      S[i] = tlchg[i].Vect();
      S[i] = getGunionPolarimeter_pi( tlchg[i], tlneu[i] );
    }

    S[i]*=1./S[i].Mag();
    // get component perpendicular to axis
    Sn[i] = S[i] - axis*(S[i].Dot(axis));
  }

  // these are longitudinal components
  theta0 = S[0].Angle(axis);
  theta1 = S[1].Angle(axis);

  // signed angle between perpendicular components
  dPhi = Sn[0].Angle(Sn[1]);
  if ( (Sn[0].Cross(Sn[1])).Dot(axis) < 0 ) dPhi = -dPhi;
  
  return;
}

TVector3 tauUtils::getGunionAxis( TLorentzVector t1, TLorentzVector t2 ) {
  // z axis is direction of tau0 in tau-tau rest frame
  TLorentzVector tlax = t1;
  tlax.Boost( - (t1+t2).BoostVector() );
  TVector3 axis = tlax.Vect();
  axis*=1./axis.Mag();
  return axis;
}

TVector3 tauUtils::getGunionPolarimeter_pi(  TLorentzVector tlchg, TLorentzVector tlneu) {
  TLorentzVector tau = tlchg+tlneu;
  // boost everything to tau's rest frame
  tlchg.Boost( - tau.BoostVector() );
  tlneu.Boost( - tau.BoostVector() );
  // the polarimeter
  return tlchg.Vect();
}

TVector3 tauUtils::getGunionPolarimeter_rho(  TLorentzVector tlchg, TLorentzVector tlpi0, TLorentzVector tlneu) {
  TLorentzVector tau = tlchg+tlpi0+tlneu;
  // boost everything to tau's rest frame
  tlchg.Boost( - tau.BoostVector() );
  tlneu.Boost( - tau.BoostVector() );
  tlpi0.Boost( - tau.BoostVector() );
  // the polarimeter
  return m_tau * ( tlchg.Vect() - tlpi0.Vect() )*( tlchg.E() - tlpi0.E() ) + tlneu.Vect() * ( (tlchg+tlpi0).M2()/2.);
}


float tauUtils::getRougeLogLikelihood( float theta0, float theta1, float dPhi, float psi , bool useFullFunction ) {
  float cc = cos(theta0)*cos(theta1);
  float ss = sin(theta0)*sin(theta1);
  float loglikelihood;
  if ( useFullFunction ) {
    loglikelihood = log( (1./(8*acos(-1)))*( 1 + cc - ss*cos(dPhi - 2*psi) ) );
  } else {
    loglikelihood = log( (1./(2*acos(-1)))*( 1 - (acos(-1)*acos(-1)/16.)*cos(dPhi - 2*psi) ) );
  }
  return loglikelihood;
}

void tauUtils::FillLogLike ( TH1F* hLike,
			     TLorentzVector charged0, TLorentzVector neutral0, TLorentzVector neutrino0,
			     TLorentzVector charged1, TLorentzVector neutral1, TLorentzVector neutrino1,
			     bool useFullFunction ) {
  float theta0, theta1, dPhi;
  rouge_getThetaDphi( charged0, neutral0, neutrino0, charged1, neutral1, neutrino1, theta0, theta1, dPhi );
  FillLogLike( hLike, theta0, theta1, dPhi, useFullFunction );
  return;
}

void tauUtils::FillLogLike ( TH1F* hLike, float theta0, float theta1, float dPhi,
			       bool useFullFunction ) {
  // add to the likelihood histogram
  for (int i=1; i<=hLike->GetNbinsX(); i++) {
    float bincent = hLike->GetXaxis()->GetBinCenter(i);
    float psi = acos(-1)*bincent; // normalise by pi
    float loglikelihood = getRougeLogLikelihood( theta0, theta1, dPhi, psi, useFullFunction );
    hLike->Fill( bincent, -loglikelihood );
  }
  return;
}


int tauUtils::isMuon( const ReconstructedParticle* rp ) {

  //// may 2016: use Pandora PID
  //if ( abs(rp->getType())==13 ) {
  //  return TIGHT;
  //} else {
  //  return NONE;
  //}


  float e_on_p, ecalFrac, ediff_sig, ecalo, emuon;
  getElectronVars( rp, e_on_p, ecalFrac, ediff_sig, ecalo, emuon );

  float p(0);
  for (int i=0; i<3; i++) p+=pow( rp->getMomentum()[i], 2 );
  p=sqrt(p);

  if ( emuon > 1.2 && 
       ecalo < 15. &&
       e_on_p < 0.3 ) {
    return TIGHT;
  } else if ( abs(rp->getType())==13 || // pandora pid
	     ( ( emuon > 1.2 || p<15 ) && 
	       ecalo < 25. &&
	       e_on_p < 0.4 ) ) {
    return LOOSE;
  } else {
    return NONE;
  }



  /*
  int nl_ecal(0);
  int nh_ecal(0);

  int nl_hcal(0);
  int nh_hcal(0);
  int nl_yoke(0);
  int nh_yoke(0);
  int nl_oth(0);
  int nh_oth(0);

  ClusterVec clv = rp->getClusters();
  for (size_t i=0; i<clv.size(); i++) {
    const Cluster* cl = clv[i];
    std::map < std::vector < int > , int > hitSummary;
    for (size_t f=0; f<cl->getCalorimeterHits().size(); f++) {
      CalorimeterHit* ch = cl->getCalorimeterHits()[f];
      CHT cht(ch->getType());
      std::vector < int > id_layout_layer;
      id_layout_layer.push_back( cht.caloID() );
      id_layout_layer.push_back( cht.layout() );
      id_layout_layer.push_back( cht.layer() );
      if ( hitSummary.find( id_layout_layer ) != hitSummary.end() ) hitSummary[id_layout_layer]++;
      else hitSummary[id_layout_layer]=1;
    }

    for ( std::map < std::vector < int > , int >::iterator jtt=hitSummary.begin(); jtt!=hitSummary.end(); jtt++) {
      std::vector < int > ff = jtt->first;
      if ( ff[0]==CHT::ecal ) {
        nl_ecal++;
        nh_ecal+=jtt->second;
      } else if ( ff[0]==CHT::hcal ) {
        nl_hcal++;
        nh_hcal+=jtt->second;
      } else if ( ff[0]==CHT::yoke ) {
        nl_yoke++;
        nh_yoke+=jtt->second;
      } else {
        nl_oth++;
        nh_oth+=jtt->second;
      }
    }
  }
  bool ismuon=
    (nl_yoke > 3) ||
    (nl_hcal > 30 && float(nh_hcal)/float(nl_hcal) < 3);


  std::cout << "hello from tauUtils::isMuon : result " <<  ismuon << " recoPart = " << rp->getType() << std::endl;

  return ismuon;

  */

}

float tauUtils::getEcalEnergyFraction( const ReconstructedParticle* rp ) {
  float ecalEn(0);
  float calEn(0);
  ClusterVec clv = rp->getClusters();
  for (size_t i=0; i<clv.size(); i++) {
    const Cluster* cl = clv[i];
    ecalEn += cl->getSubdetectorEnergies()[0];
    calEn  += cl->getEnergy();
  }
  return ecalEn/calEn;
}



int tauUtils::isPhoton( const ReconstructedParticle* rp ) {
  // just use pandora ID
  //  return rp->getType()==22; //  && rp->getMass()<0.01; // why did I use this mass requirement?

  if ( rp->getType()==22 ) return TIGHT;
  else return NONE;
}

void tauUtils::getElectronVars( const ReconstructedParticle* rp, 
				float & EonP, // energy / momentum
				float & ecalFrac, // fraction in ECAL
				float & epDiffRel,
				float & eCalo,
				float & eMuon
				) { // enery/mom diff in terms of calo resolution

  float trkMom(0);
  if ( rp->getTracks().size()>0 ) {
    for (int k=0; k<3; k++) trkMom += pow(rp->getMomentum()[k], 2);
    trkMom=sqrt(trkMom);
  }

  //  std::cout << "ecal frac: " << getEcalEnergyFraction( rp ) << std::endl;

  float e_ecal(0);
  eCalo = 0;
  eMuon = 0;

  ClusterVec clv = rp->getClusters();
  for (size_t i=0; i<clv.size(); i++) {
    const Cluster* cl = clv[i];
    e_ecal+=cl->getSubdetectorEnergies()[0];
    eMuon+=cl->getSubdetectorEnergies()[2];
    for (size_t k=0; k<cl->getSubdetectorEnergies().size(); k++) {
      eCalo+=cl->getSubdetectorEnergies()[k];
    }
  }

  float de_ecal = sqrt( pow( sqrt(e_ecal)*0.17, 2 ) + pow( 0.01*e_ecal, 2 ) ); // expected EM calo resolution

  EonP = trkMom > 0 ? eCalo/trkMom : -1;
  ecalFrac = e_ecal/eCalo;
  epDiffRel = trkMom>0 && de_ecal>0 ? (e_ecal - trkMom)/de_ecal : -999;

  return;
}

int tauUtils::isElectron( const ReconstructedParticle* rp ) {

  float e_on_p, ecalFrac, ediff_sig, ecalo, emuon;
  getElectronVars( rp, e_on_p, ecalFrac, ediff_sig, ecalo, emuon );

  int electrontag=NONE;
  if ( ecalFrac>0.98 && // little HCAL contrib.                                                                             
       ( (e_on_p>0.85 && e_on_p<1.1) ||   // good E/p conpatibility (absolute terms)                                           
	 ediff_sig < 5 )                   //        compared to expected ecal en res                                           
       ) {
    electrontag=TIGHT;
  } else   if ( ecalFrac>0.95 && // little HCAL contrib.                                                                             
       ( (e_on_p>0.65 && e_on_p<1.3) ||   // good E/p conpatibility (absolute terms)    // daniel loosened from 0.75/1.2
	 ediff_sig < 8 )                   //        compared to expected ecal en res                                           
       ) {
    electrontag=LOOSE;
  };
  return electrontag;
}


bool tauUtils::isPdg( ReconstructedParticle* rp , const int pdg , UTIL::PIDHandler& pidhand ) {
  bool match=false;

//   cout << " -------------------- " << endl;
//   cout << "isPdg? energy " << rp->getEnergy() << " pandora pid " << rp->getType() << endl;
//   ParticleIDVec allpids = rp->getParticleIDs ();
//   cout << "number of pids: " << allpids.size() << endl;
//   for (size_t i=0; i< allpids.size(); i++) {
//     ParticleID* pid = allpids[i];
//     cout << i << " " << pid << " type " << pid->getType() << " PDG " << pid->getPDG() << " likelihood " << pid->getLikelihood () << " algotype " << pid->getAlgorithmType () << endl;
//   }

//  int lpdg = pidhand.getParticleID( rp, pidhand.getAlgorithmID( "LikelihoodPID" ) ).getPDG();
//
//  int showershapePid =  pidhand.getParticleID( rp, pidhand.getAlgorithmID( "ShowerShapesPID" ) ).getPDG();
//
//  int basicPid = pidhand.getParticleID( rp, pidhand.getAlgorithmID( "BasicVariablePID" ) ).getPDG();

//  if ( abs(lpdg)==abs(pdg) ) match=true;

  //  cout << "requested pdg: " << pdg << "  basic/shower/likelihoodPID result: " << basicPid << " " << showershapePid << " " << lpdg << " match? " << match << " pandora pid: " <<  rp->getType() << endl;

  return match;
}

std::map < const ReconstructedParticle* , std::vector < const ReconstructedParticle* > > tauUtils::find_FSR_brems( const std::vector < const ReconstructedParticle* > & chargedpfos , 
														   const std::vector < const ReconstructedParticle* > & photonpfos    ) {

  // find brems/FSR candidates near to charged tracks

  std::map < const ReconstructedParticle* , std::vector < const ReconstructedParticle* > > addedPhotons;
  std::vector < const ReconstructedParticle* > associatedPhotons;

  // look for brems/FSR photons
  for ( size_t i=0; i<chargedpfos.size(); i++) {

    TVector3 trkmom( chargedpfos[i]->getMomentum()[0], chargedpfos[i]->getMomentum()[1], chargedpfos[i]->getMomentum()[2] );

    bool possibleElectron = 
      chargedpfos[i]->getClusters().size()>0 && 
      getEcalEnergyFraction( chargedpfos[i] ) > 0.98;
    TVector3 clpos;
    float aveth(-999), phimin(-999), phimax(-999);
    if ( possibleElectron ) {
      // find biggest cluster
      float emax(-1);
      Cluster* clmax(0);
      for (size_t j=0; j<chargedpfos[i]->getClusters().size(); j++) {
	if ( chargedpfos[i]->getClusters()[j]->getEnergy()>emax ) {
	  emax = chargedpfos[i]->getClusters()[j]->getEnergy();
	  clmax = chargedpfos[i]->getClusters()[j];
	}
      }
      clpos.SetXYZ( clmax->getPosition()[0], clmax->getPosition()[1], clmax->getPosition()[2] );

      phimin = std::min( trkmom.Phi(), clpos.Phi());
      phimax = std::max( trkmom.Phi(), clpos.Phi());
      if ( fabs( trkmom.Theta() - clpos.Theta() ) < 0.01 ) { // consistent theta between track and cluster
	aveth = ( trkmom.Theta() + clpos.Theta() )/2.;
      } else {                                               // not such a clean case, reject this candidate
	possibleElectron = false;
      }
    }

    for (size_t k=0; k<photonpfos.size(); k++) {

      //      cout << "--- chg " << i << " gamma " << k << endl;

      TVector3 gamdir ( photonpfos[k]->getMomentum()[0], photonpfos[k]->getMomentum()[1], photonpfos[k]->getMomentum()[2] );

      bool matched=false;

      // first look for FSR, 
      // angle of photon direction to estimated momentum @ IP
      float angle = trkmom.Angle( gamdir );
      if ( cos(angle)>0.999 ) { // looks like FSR
	if ( find( associatedPhotons.begin(), associatedPhotons.end(), photonpfos[k] ) != associatedPhotons.end() ) {
	  //	  cout << "WARNING, this photon alread assigned to someone else!" << endl;
	} else {
	  //	  cout << "matched FSR " << angle << " " << cos(angle) << endl;
	  matched=true;
	  if ( addedPhotons.find( chargedpfos[i] ) != addedPhotons.end() ) 
	    addedPhotons[chargedpfos[i]].push_back(photonpfos[k]);
	  else {
	    std::vector < const ReconstructedParticle* > pp; pp.push_back(photonpfos[k]);
	    addedPhotons[chargedpfos[i]]=pp;
	  } 
	  associatedPhotons.push_back(photonpfos[k]);
	}
      } // FSR

      if ( !matched && possibleElectron ) { // track may be electron, look for brems
	bool matchesTheta = fabs( gamdir.Theta() - aveth ) < 0.01 ;
	bool matchesPhi(false);
	if ( phimax - phimin < acos(-1) ) {
	  if ( gamdir.Phi()>=phimin && gamdir.Phi()<=phimax ) matchesPhi=true;
	} else {
	  if ( gamdir.Phi()>=phimax || gamdir.Phi()<=phimin ) matchesPhi=true;
	}
	if ( matchesTheta && matchesPhi ) {
	  //	  cout << " good brems candidate! " << photonpfos[k]->getEnergy() << endl;
	  if ( find( associatedPhotons.begin(), associatedPhotons.end(), photonpfos[k] ) != associatedPhotons.end() ) {
	    //	    cout << "WARNING, this photon alread assigned to someone else!" << endl;
	  } else {
	    if ( addedPhotons.find( chargedpfos[i] ) != addedPhotons.end() )
	      addedPhotons[chargedpfos[i]].push_back(photonpfos[k]);
	    else {
	      std::vector < const ReconstructedParticle* > pp; pp.push_back(photonpfos[k]);
	      addedPhotons[chargedpfos[i]]=pp;
	    }
	    associatedPhotons.push_back(photonpfos[k]);
	  }
	}


      } // possible electron, brems

    } // photons

    // if ( addedPhotons.find(chargedpfos[i])!=addedPhotons.end() && addedPhotons[chargedpfos[i]].size()>0 ) {
    //   float addedEnergy(0);
    //   for (size_t kk=0; kk<addedPhotons[chargedpfos[i]].size(); kk++) {
    // 	addedEnergy+=addedPhotons[chargedpfos[i]][kk]->getEnergy();
    //   }
    //   cout << "BREMS: initial momentum, energy = " << trkmom.Mag() << " " << chargedpfos[i]->getEnergy() << " added energy " << addedEnergy << endl;
    // }
    
  } // charged

  return addedPhotons;
}


std::vector < MCParticle* > tauUtils::getstablemctauDaughters( MCParticle* mctau ) {
  return getmctauDaughters ( mctau, true );
}


std::vector < MCParticle* > tauUtils::getmctauDaughters( MCParticle* mctau, bool onlystable ) {
  assert( abs( mctau->getPDG())==15 );
  std::vector < MCParticle* > daughters;
  std::vector < MCParticle* > pp;
  pp.push_back(mctau);
  while ( pp.size()>0 ) {
    std::vector < MCParticle* > pp2;
    for ( size_t i=0; i< pp.size(); i++) {
      if ( pp[i]->getGeneratorStatus() == 1 ) {
	daughters.push_back( pp[i] );
      } else {
	if ( !onlystable ) daughters.push_back( pp[i] );
	for ( size_t k=0; k< pp[i]->getDaughters().size(); k++) {
	  pp2.push_back( pp[i]->getDaughters()[k] );
	}
      }
    }
    pp=pp2;
  }
  return daughters;
}

TLorentzVector tauUtils::getTLV( std::vector < ReconstructedParticle* > rps ) {
  TLorentzVector tot(0,0,0,0);
  for ( size_t k=0; k<rps.size(); k++) 
    tot+=getTLV(rps[k]);
  return tot;
}

TLorentzVector tauUtils::getTLV( ReconstructedParticle* rp ) {
  assert (rp );

  TLorentzVector fmom( rp->getMomentum()[0],
		       rp->getMomentum()[1],
		       rp->getMomentum()[2],
		       rp->getEnergy() );

  return fmom;

}

TLorentzVector tauUtils::getTLV( MCParticle* mc ) {
  assert ( mc );
  return TLorentzVector( mc->getMomentum()[0],
			 mc->getMomentum()[1],
			 mc->getMomentum()[2],
			 mc->getEnergy() );
}

MCParticle* tauUtils::getBestTrackMatch(  ReconstructedParticle* pfo, UTIL::LCRelationNavigator* relNavi ) {
  MCParticle* bestmatch(0);
  float bestTrkWeight(0);
  for (size_t k=0; k<relNavi->getRelatedToObjects( pfo ).size(); k++) {
    MCParticle* mcp = dynamic_cast <MCParticle*> (relNavi->getRelatedToObjects(pfo)[k]);
    float wgt = relNavi->getRelatedToWeights(pfo)[k];
    float trackwgt = (int(wgt)%10000)/1000.;
    if (trackwgt>bestTrkWeight) {
      bestmatch=mcp;
      bestTrkWeight=trackwgt;
    }
  }
  return bestmatch;
}

MCParticle* tauUtils::getBestCaloMatch(  ReconstructedParticle* pfo, UTIL::LCRelationNavigator* relNavi ) {
  MCParticle* bestmatch(0);
  float bestCalWeight(0);
  for (size_t k=0; k<relNavi->getRelatedToObjects( pfo ).size(); k++) {
    MCParticle* mcp = dynamic_cast <MCParticle*> (relNavi->getRelatedToObjects(pfo)[k]);
    float wgt = relNavi->getRelatedToWeights(pfo)[k];
    float clusterwgt = (int(wgt)/10000)/1000. ;
    if (clusterwgt>bestCalWeight) {
      bestmatch=mcp;
      bestCalWeight=clusterwgt;
    }
  }
  return bestmatch;
}



int tauUtils::getMCdecayMode( std::vector < MCParticle* > mcps ) {

  int npichg(0);
  int nel(0);
  int nmu(0);
  int nnu(0);
  int ngamma(0);
  int nK0L(0);
  int nKchg(0);

  int nChg(0);
  int nNeu(0);

  for ( size_t i=0; i<mcps.size(); i++) {
    switch ( abs( mcps[i]->getPDG() ) ) {
    case 11:
      nel++; 
      nChg++;
      break;
    case 13:
      nmu++; 
      nChg++;
      break;
    case 12:
    case 14:
    case 16:
      nnu++; 
      break;
    case 211:
      npichg++; 
      nChg++;
      break;
    case 22:
      ngamma++; 
      nNeu++;
      break;
    case 130:
      nK0L++;
      nNeu++;
      break;
    case 321:
      nKchg++; 
      nChg++;
      break;
    default:
      cout << "unknown pdg " << mcps[i]->getPDG() << endl;
      assert(0);
    }
  }



  int decayMode(-1);
  if ( nel==1 ) decayMode = decayEl;
  else if ( nmu==1 ) decayMode = decayMu;
  else if ( npichg==nChg && npichg==1 && nNeu==0 ) decayMode = decayChPi;
  else if ( npichg==nChg && npichg==1 && ngamma==2 && nNeu==2 ) decayMode = decayRho;
  else if ( npichg==nChg && npichg==3 && nNeu==0 ) decayMode = decayA1_3p;
  else if ( npichg==nChg && npichg==1 && ngamma==nNeu && ngamma==4 ) decayMode = decayA1_1p;
  else if ( nK0L>0 || nKchg>0 ) decayMode = decayK;
  else if ( npichg+nKchg>0 ) decayMode = decayW; // everything else...
  else {
    cout << "unrecognised tau decay!" << endl;
    cout << " pi " <<  npichg;
    cout << " e " <<  nel;
    cout << " m " <<  nmu;
    cout << " nu " <<  nnu;
    cout << " g " <<  ngamma;
    cout << " kl " <<  nK0L;
    cout << " k " <<  nKchg;
    cout << " chg " <<  nChg;
    cout << " neu " <<  nNeu;
    cout << endl;    
    assert(0);
  }

  return decayMode;
}

TString tauUtils::getTauDecLab( int code ) {
  //  enum { decayChPi=0, decayRho, decayA1_3p, decayA1_1p , decayEl, decayMu , decayW , decayK , decayMultiprong , decayOthersingleprong, decayUnrecognised};
  switch (code) {
  case decayChPi:
    return "PI";
  case decayRho:
    return "RHO";
  case decayA1_3p:
    return "A3P";
  case decayA1_1p:
    return "A1P";
  case decayEl:
    return "E";
  case decayMu:
    return "M";
  case decayW:
    return "HAD";
  case decayK:
    return "K";
  case decayMultiprong:
    return "MP";
  case decayOthersingleprong:
    return "SP";
  case decayUnrecognised:
  default:
    return "OTH";
  }

}


std::vector < EVENT::MCParticle* > tauUtils::findFinalTaus( EVENT::LCCollection* mccol ) {

  std::vector < EVENT::MCParticle* > finalmctaus;

  for (int j=0; j<mccol->getNumberOfElements(); j++) {
    MCParticle* mcp =  dynamic_cast<MCParticle*> (mccol->getElementAt(j));
    if ( abs( mcp->getPDG())==15 ) {
      std::vector < MCParticle* > daughters = mcp->getDaughters();
      bool hasfunnydaughter(false);
      for ( size_t i=0; i<daughters.size(); i++) {
	if ( abs(daughters[i]->getPDG())==15 || abs(daughters[i]->getPDG())==94 ) {
	  hasfunnydaughter=true;
	  break;
	}
      }
      if ( ! hasfunnydaughter ) finalmctaus.push_back( mcp );
    }
  }

  return finalmctaus;
}
