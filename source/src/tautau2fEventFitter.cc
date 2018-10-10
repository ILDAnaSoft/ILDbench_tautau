#include "tautau2fEventFitter.h"

#include "tauUtils.h"

#include <iostream>
using std::cout;
using std::endl;

const float tautau2fEventFitter::m_tau = 1.777;
const float tautau2fEventFitter::ctau_tau = 0.087;

void tautau2fEventFitter::prepareForFit() {

  assert ( _ecom > 0 );

  for (int itau=0; itau<2; itau++) {

    TVector3 ipToTraj = chargedHadronImpactParameterVector[itau]; // since in MC the IP is really @ 0,0,0, this is the true vector
    ipToTraj -= _ip; // make it wrt to IP (e.g. to investigate the effect of IP mismeasurement)

    //    chargedHadronEVector[itau] = getEVector( chargedHadronMomentum[itau].Vect(), chargedHadronImpactParameterVector[itau] );
    chargedHadronEVector[itau] = getEVector( chargedHadronMomentum[itau].Vect(), ipToTraj );

    normalToTrackPlane[itau] = chargedHadronMomentum[itau].Vect().Cross( chargedHadronEVector[itau] );
    normalToTrackPlane[itau]*=1./normalToTrackPlane[itau].Mag(); // normal to decay plane

    kPerp[itau] = (neutralHadronMomentum[itau].Vect().Dot(normalToTrackPlane[itau]))*normalToTrackPlane[itau];
    hadMom[itau] = chargedHadronMomentum[itau] + neutralHadronMomentum[itau]; // total hadronic 4mom

    // the two axes in the plane
    hPlane[itau] = hadMom[itau].Vect() - (hadMom[itau].Vect().Dot(normalToTrackPlane[itau]))*normalToTrackPlane[itau]; 
    hPlane[itau]*=1./hPlane[itau].Mag();

    fPlane[itau] = normalToTrackPlane[itau].Cross(hPlane[itau]); 
    fPlane[itau]*=1./fPlane[itau].Mag();

    if (verbose) {
      cout << "normalToTrackPlane = " << normalToTrackPlane[itau].X() << " " << normalToTrackPlane[itau].Y() << " " << normalToTrackPlane[itau].Z() << endl;
      cout << "kPerp = " << kPerp[itau].X() << " " << kPerp[itau].Y() << " " << kPerp[itau].Z() << endl;
    }
  }

  prepared=true;

  return;
}


bool tautau2fEventFitter::fitIt_single_single(int bestStrategy) {

  cout << "hello from fitIt_single_single, strategy " << bestStrategy << endl;

  // both taus are treated as single prong
  if (!functor) {
    functor = new ROOT::Math::Functor( this, & tautau2fEventFitter::fitfn_single_single, 3);
    minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "tautau_ptminISR");
    minimizer->SetMaxFunctionCalls(10000); // for Minuit/Minuit2
    minimizer->SetTolerance(0.1);
    minimizer->SetPrintLevel(-1);
    minimizer->SetFunction(*functor);
  }

  const double piOn2 = acos(-1)/2.;
  double bestPt(99999999);
  double bestPz(99999999);
  double bestMassDiff(9999999);

  bool bestSol[3]={0,0,0};

  //  int bestQuad=-1;
  double bestPsi[2]={-999,-999};
  double bestZ(-999);
  // loop over 4 energy solution combinations, and 4 quadrants
  for (int isol=0; isol<4; isol++) { // the two Energy solutions per tau at given psi
    isolution[0] = isol>1;
    isolution[1] = isol%2==1;
    for (int iquad=0; iquad<4; iquad++) { // the 4 quadrants bounded by psi1,2=0;

      // set the quadrant limits and starting points
      bool pos0 = iquad<2;
      bool pos1 = iquad%2==0;

      double psi0min = pos0 ? 0 : -piOn2;
      double psi0max = pos0 ? piOn2 : 0;

      double psi1min = pos1 ? 0 : -piOn2;
      double psi1max = pos1 ? piOn2 : 0;

      double step = 0.03;
      double start0 = pos0 ? step : -step;
      double start1 = pos1 ? step : -step;
    
      // Set the free variables to be minimized! start, step, min, max
      minimizer->SetLimitedVariable(0,"psi0",start0, step, psi0min, psi0max);
      minimizer->SetLimitedVariable(1,"psi1",start1, step, psi1min, psi1max);

      double avez = 0.5*(chargedHadronImpactParameterVector[0].Z() + chargedHadronImpactParameterVector[1].Z() );

      minimizer->SetLimitedVariable(2,"ipz",avez, 1., -1, 1);

      // do the minimization
      minimizer->Minimize();

      double thispt = minimizer->MinValue(); // not really pt any more...include ISR contribution (pt**2 + photon mass**2)

      const double *ffxs = minimizer->X();
      setPsi( 0, ffxs[0] );
      setPsi( 1, ffxs[1] );
      setIPz( ffxs[2] );

      cout << "SOLUTION " << isol << " quadrant " << iquad << " BEST FIT: pars " << ffxs[0] << " " << ffxs[1] << " " << ffxs[2] << " value " << thispt << endl;

      getEvent4Momentum( isolution[0], isolution[1] ).Print();

      cout << "decay lengths: " << getDecayLength( 0, isolution[0] ) << " " << getDecayLength( 1, isolution[1] ) << endl;


      cout << getTau4Momentum(0, isolution[0]).M() << " " << 
	getTau4Momentum(1, isolution[1]).M() << " " <<  
	(getTau4Momentum(0, isolution[0])+getTau4Momentum(1, isolution[1])).M() << endl;

      const float ptcutoff = 1.0;

      if ( getEvent4Momentum( isolution[0], isolution[1] ).E() < _ecom+10 && 
	   getDecayLength( 0, isolution[0] )>0 && getDecayLength( 1, isolution[1] )>0 ) { // positive decay length

	float thispz = fabs( getEvent4Momentum( isolution[0], isolution[1] ).Pz() );
	float thisMtt = (getTau4Momentum(0, isolution[0])+getTau4Momentum(1, isolution[1])).M();
	float thismassdiff = fabs( thisMtt - 125.0 );

	// cout << "event energy: sol " << isol << " quad " << iquad << " isol01 " << isolution[0] << " " << isolution[1] << 
	//   " totE " << getEvent4Momentum( isolution[0], isolution[1] ).E() << " pt " << thispt << " pz " << thispz << " thismassdiff " << thismassdiff << endl;

	bool better(false);
	if ( bestStrategy==0 ) { // smallest pt
	  if ( thispt<bestPt ) better=true;
	} else if ( bestStrategy==1 ) { //smallest pt if minpt>cutoff, smallest pz if minpt<cutoff
	  if ( bestPt>ptcutoff ) {
	    if ( thispt<bestPt ) better=true;
	  } else {
	    if ( thispt<ptcutoff && thispz<bestPz ) better=true;
	  }
	} else if ( bestStrategy==2 ) { //smallest pt if minpt>cutoff, closest mtt to mH if minpt<cutoff
	  if ( bestPt>ptcutoff ) {
	    if ( thispt<bestPt ) better=true;
	  } else {
	    if ( thispt<ptcutoff && thismassdiff < bestMassDiff ) better=true;
	  }
	}
	if ( better ) {

	  bestSol[0] = isolution[0];
	  bestSol[1] = isolution[1];
	  //	  bestQuad = iquad;
	  bestPt=thispt;
	  bestPz=thispz;
	  bestMassDiff = thismassdiff;
	  for (int j=0; j<2; j++) bestPsi[j]=ffxs[j];
	  bestZ = ffxs[2];
	  //	  cout << " better..." << bestSol[0] << " " << bestSol[1] << " " << bestQuad << " " << bestPsi[0] << " " << bestPsi[1] << endl;
	}
      }
    }
  }

  if (bestPt>_ecom) return false;

  // set the parameters to be best ones
  bool goodSolution=true;
  for (int i=0; i<2; i++) {
    isolution[i]=bestSol[i];
    goodSolution = goodSolution && setPsi(i, bestPsi[i]) && setIPz(bestZ);
  }

//  cout << "good soln? " << goodSolution << " ; " << bestPsi[0] << " " << bestPsi[1] << " " << getEvent4Momentum( bestSol[0], bestSol[1] ).E() << endl;  
//
//  cout << "bye from fitIt_single_single " << goodSolution << endl;

  goodfit = goodSolution;

  return goodSolution;
}

double tautau2fEventFitter::fitfn_single_single(const double* pars) {

  bool ok(true);
  for (int i=0; i<2; i++)
    ok = ok && setPsi(i,pars[i]);
  ok = ok && setIPz( pars[2] );
  double val = ok ? getPt_plusISR(isolution[0], isolution[1]) : 1000*(1 + sqrt(pow(pars[0],2) + pow(pars[1],2) + pow(pars[2],2) ) );

  //  cout << "SOL: " << isolution[0] << " " << isolution[1] << " pars: " << pars[0] << " " << pars[1] << " " << pars[2] << " OK? " << ok << " netmom: " << val << endl;


  return val;

}


bool tautau2fEventFitter::setIPz( double z) {
  _ip.SetXYZ(0,0,z);
  return true;
}

bool tautau2fEventFitter::setPsi(int itau, double psi) {
  assert( itau>=0 && itau<ntau ); 
  tauPsi[itau]=psi;
  return calculateNeutrino(itau);
}

bool tautau2fEventFitter::calculateNeutrino(int itau) {
  assert( itau>=0 && itau<ntau );

  //  if (!prepared) { // since ipz can change, we have to do this everytime
  prepareForFit(); 
  //  }

  if ( hadMom[itau].M() > m_tau ) return false;

  TVector3 qhat = cos(tauPsi[itau])*hPlane[itau] + sin(tauPsi[itau])*fPlane[itau]; // direction of neutrino in plane

  double A = ( pow(m_tau,2) - hadMom[itau].M2() ) - 2*hadMom[itau].Vect().Dot(kPerp[itau]);
  double B = 2*qhat.Dot(hadMom[itau].Vect()-kPerp[itau]);
  double C = 2*hadMom[itau].E();

  double a = (B*B-C*C);
  double b = 2*A*B;
  double c = A*A - C*C*kPerp[itau].Mag2();

  double X = b*b - 4*a*c;
  if ( X < 0 ) { // not physical
    if (verbose) cout << "unphysical..." << endl;
    for (int i=0; i<2; i++) {
      neutrinoMomentum[itau][i].SetXYZT(0,0,0,0);
    }
    return false;
  } else {
    X=sqrt(X);
    for (int i=0; i<2; i++) {
      double QQ = -b;
      QQ += i==0 ? -X : X;
      QQ/=2*a;
      neutrinoMomentum[itau][i].SetVectM( QQ*qhat - kPerp[itau], 0);
    }
  }
  return true;
}

TLorentzVector tautau2fEventFitter::getEventVisible4Momentum() {
  TLorentzVector ptot4(0,0,0,0);
  ptot4+=recoil;
  for (int i=0; i<ntau; i++) {
    ptot4+=neutralHadronMomentum[i];
    ptot4+=chargedHadronMomentum[i];
  }
  return ptot4;
}


TLorentzVector tautau2fEventFitter::getEvent4Momentum(bool j0, bool j1) {
  assert( j0>=0 && j0<2 );
  assert( j1>=0 && j1<2 );
  TLorentzVector ptot4 = getEventVisible4Momentum();
  ptot4+=neutrinoMomentum[0][(int) j0];
  ptot4+=neutrinoMomentum[1][(int) j1];
  return ptot4;
}

TLorentzVector tautau2fEventFitter::getEvent4Momentum() {
  return getEvent4Momentum(isolution[0], isolution[1]);
}


double tautau2fEventFitter::getPt(bool j0, bool j1) {
  return getEvent4Momentum(j0, j1).Pt();
}

double tautau2fEventFitter::getPt_plusISR(bool j0, bool j1) {

  double pt2 = getEvent4Momentum(j0, j1).Perp2() +        // event pt2
    pow( _ecom - getEvent4Momentum(j0, j1).E() - fabs( getEvent4Momentum(j0, j1).Pz() ), 2 ) ; // net pz**2 (including/assuming single ISR)

  return sqrt( pt2 );

}

double tautau2fEventFitter::getPt() {
  return getPt(isolution[0], isolution[1]);
}

TLorentzVector tautau2fEventFitter::getNeutrino4Momentum(int itau, bool jsol) {
  return neutrinoMomentum[itau][(int) jsol];
}

TLorentzVector tautau2fEventFitter::getNeutrino4Momentum(int itau) {
  return getNeutrino4Momentum( itau, isolution[itau]);
}

TLorentzVector tautau2fEventFitter::getTau4Momentum(int itau, bool jsol) {
  return neutralHadronMomentum[itau] + chargedHadronMomentum[itau] + neutrinoMomentum[itau][(int) jsol];
}

TLorentzVector tautau2fEventFitter::getTau4Momentum(int itau) {
  return getTau4Momentum(itau, isolution[itau]);
}


double tautau2fEventFitter::getDecayLength(int itau, bool jsol) {

  // decay length of tau in lab

  TVector3 tau3momDir = getTau4Momentum(itau, jsol).Vect();
  tau3momDir*=1./tau3momDir.Mag();

  TVector3 eVector = chargedHadronEVector[itau]; // chargedHadronImpactParameterVector[itau];
  double eMag = eVector.Mag();
  eVector*=1./eMag;

  // decay length = |e| / ( tauHat dot dHat );
  return eMag / ( tau3momDir.Dot( eVector ) );

}

double tautau2fEventFitter::getDecayLength(int itau) {
  return getDecayLength(itau, isolution[itau]);
}

double tautau2fEventFitter::getDecayTimePDF(int itau, bool jsol) {
  double beta = getTau4Momentum(itau, jsol).Beta();  
  double gamma = getTau4Momentum(itau, jsol).Gamma();  
  double L = getDecayLength(itau, jsol);
  return L>0 ? exp(-L/(beta*gamma*ctau_tau)) : 0;
}

double tautau2fEventFitter::getDecayTimePDF(int itau) {
  return getDecayTimePDF(itau, isolution[itau]);
}

void tautau2fEventFitter::init() {

  prepared=false;

  goodfit=false;

  _ecom=0;

  _ip.SetXYZ(0,0,0);
  recoil.SetXYZT(0,0,0,0);
  for (int i=0; i<2; i++) {
    isolution[i]=false;
    //    multiProng[i]=false;
    //    tauVertices[i].SetXYZ(0,0,0);
    neutralHadronMomentum[i].SetXYZT(0,0,0,0);
    chargedHadronMomentum[i].SetXYZT(0,0,0,0);
    chargedHadronImpactParameterVector[i].SetXYZ(0,0,0);
    chargedHadronEVector[i].SetXYZ(0,0,0);
  }
  //  f_a1_nuE=NULL;
  verbose=false;
}

TVector3 tautau2fEventFitter::getEVector( TVector3 momentum, TVector3 impactVector ) {
  // impactVector not necessarily perp to momentum in 3d (often 2-d impact vector)
  // this calculates the 3d impact vector, assuming linear trajectories
  TVector3 evec;
  evec = momentum.Cross(impactVector.Cross(momentum)); // guaranteed perp to momentum, and in plane of momentum and impactVector
  // now make sure that size of eVec corresponds to distance from IP to real 3d PCA
  double theta = impactVector.Angle(evec); // angle between e and d
  double modE = impactVector.Mag()*cos(theta);   // dist from ip to 3d PCA
  evec*=modE/evec.Mag();

  return evec;
}


void tautau2fEventFitter::calculate_multi_pdf(float ipz) {


  for (int i=0; i<ntau; i++) {
    for (int j=0; j<2; j++) {
      multi_noip_valid_sol[i][j]=false;
      multi_noip_life_pdf[i][j]=0;
      multi_noip_nuen_pdf[i][j]=0;
    }
  }

  TVector3 ipPos(0,0,ipz);

  for (int itau=0; itau<2; itau++) { // loop over taus
    TVector3 tauDir = tauVertices[itau] - ipPos; tauDir = tauDir.Unit();

    TLorentzVector hh = neutralHadronMomentum[itau] + chargedHadronMomentum[itau];

    double a = pow(hh.E(),2) - pow( tauDir.Dot(hh.Vect()), 2);
    double b = tauDir.Dot(hh.Vect()) * ( pow(tauUtils::m_tau,2) + hh.M2() );
    double c = pow( pow(tauUtils::m_tau,2) + hh.M2() , 2 )/4. + pow( tauUtils::m_tau*hh.E(), 2 );
    double ss = b*b - 4*a*c;

    if ( ss>=0 ) {
      for (int jsol=0; jsol<2; jsol++) { // the two possible solutions
	double t = (b + pow(-1, jsol)*sqrt(ss))/(2*a);
	if ( t > 0 ) { // positive decay length
	  multi_noip_valid_sol[itau][jsol]=true;

	  neutrinoMomentum[itau][jsol].SetVectM( t*tauDir - hh.Vect(), 0 );

	  TVector3 td = tauVertices[itau] - ipPos;
	  TLorentzVector tt = getTau4Momentum(itau, jsol);
	  float beta  = tt.Beta();  
	  float gamma = tt.Gamma();  
	  float L     = td.Mag();
	  multi_noip_life_pdf[itau][jsol] = exp(-L/(beta*gamma*tauUtils::ctau_tau));

          TLorentzVector nn = getNeutrino4Momentum(itau, jsol);
          nn.Boost( -tt.BoostVector() );
	  multi_noip_nuen_pdf[itau][jsol] = getA1_Enu_pdf(nn.E());
	} else {
	  neutrinoMomentum[itau][jsol].SetXYZT(0,0,0,0);
	}
      } // two solutions
    } // real solutions

  } // taus
  return;
}

float tautau2fEventFitter::getA1_Enu_pdf(float nuen) {

  if (!f_a1_nuE) {
    f_a1_nuE = new TF1("f_a1_nuE","gaus(0)+gaus(3)",-10,10);
    f_a1_nuE->SetParameter(0,  3.01986e-02);
    f_a1_nuE->SetParameter(1,  5.60681e-01);
    f_a1_nuE->SetParameter(2,  6.59832e-02);
    f_a1_nuE->SetParameter(3,  1.86331e-02);
    f_a1_nuE->SetParameter(4,  4.65307e-01);
    f_a1_nuE->SetParameter(5,  1.05306e-01);
  }
  return f_a1_nuE->Eval(nuen);
}
