#include "danielKeitaTauFinder.h"
#include <iostream>
#include <map>

#include <marlin/Global.h>
#include "lcio.h"

#include "EVENT/LCCollection.h"
#include "EVENT/Track.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/ReconstructedParticleImpl.h"

#include "TVector3.h"
#include "TMath.h"

#include "vertexInfo.h"

using namespace lcio ;
using namespace marlin ;
using std::cout;
using std::endl;

danielKeitaTauFinderProcessor adanielKeitaTauFinderProcessor ;


danielKeitaTauFinderProcessor::danielKeitaTauFinderProcessor() : Processor("danielKeitaTauFinderProcessor") {
  
  // modify processor description
  _description = "danielKeitaTauFinderProcessor does whatever it does ..." ;

  registerProcessorParameter("outputFilename",
                             "name of output file",
                             _outfile,
                             std::string( "findTau.root") );

  _conesize = 0.1; 

}


void danielKeitaTauFinderProcessor::init() { 
  cout << "hello from danielKeitaTauFinderProcessor::init" << endl;

  //  decLab[TDEC_E  ] = "E";
  //  decLab[TDEC_M  ] = "M";
  //  decLab[TDEC_PI ] = "PI";
  //  decLab[TDEC_RHO] = "RHO";
  //  decLab[TDEC_A3P] = "A3P";
  //  decLab[TDEC_A1P] = "A1P";
  //  decLab[TDEC_HAD] = "HAD";
  //  decLab[TDEC_K  ] = "K";  

  _fout = new TFile(_outfile.c_str(),"recreate");
  h_gamgam_mass = new TH1F( "ggmass", "ggmass", 200, 0, 0.4 );

  h_ttmass = new TH1F( "mcttmass", "mcttmass", 205, 0, 510 );
  hSEL_ttmass = new TH1F( "SEL_mcttmass", "SEL_mcttmass", 205, 0, 510 );

  for (int iss=0; iss<4; iss++) {
    TString samp="sample";
    samp+=iss; samp+="_";

    h_ttmass_prongAngle[iss] = new TH2F( samp+ "ttmass_prongAngle", samp+"ttmass_prongAngle", 200,0,600,300,0,3.2);
    h_ttmass_outsideEnergy[iss] = new TH2F( samp+ "ttmass_outsideEnergy", samp+"ttmass_outsideEnergy",200,0,600,100,0,100);
    h_ttmass_insideEnergy[iss] = new TH2F( samp+ "ttmass_insideEnergy", samp+"ttmass_insideEnergy",200,0,600,600,0,600);
    h_ttmass_outsidePt[iss] = new TH2F( samp+ "ttmass_outsidePt", samp+"ttmass_outsidePt",200,0,600,100,0,100);
    h_tchg_tmass[iss] = new TH2F( samp+ "tchg_tmass", samp+"tchg_tmass", 5,-2.5,2.5,100,0,5 );
    h_ngam_tmass[iss] = new TH2F( samp+ "ngam_tmass", samp+"ngam_tmass", 10,-0.5,9.5,100,0,5 );
    h_nchg_tmass[iss] = new TH2F( samp+ "nchg_tmass", samp+"nchg_tmass", 10,-0.5,9.5,100,0,5 );
    h_nnhad_tmass[iss] = new TH2F( samp+ "nnhad_tmass", samp+"nnhad_tmass", 10,-0.5,9.5,100,0,5 );
    h_ngam_nchg[iss] = new TH2F( samp+ "ngam_nchg", samp+"ngam_nchg", 10,-0.5,9.5,10,-0.5,9.5 );

    h_mcTau_costh[iss]  = new TH2F(  samp+ "MCtau_costh", samp+ "MCtau_costh", 100, 0, 1, 100, 0, 1 );

    hSEL_tchg_tmass[iss]           = new TH2F( samp+ "SEL_tchg_tmass", samp+"SEL_tchg_tmass", 5,-2.5,2.5,100,0,5 );
    hSEL_ttmass_prongAngle[iss]    = new TH2F( samp+ "SEL_ttmass_prongAngle",   samp+"SEL_ttmass_prongAngle", 200,0,600,300,0,3.2);
    hSEL_ttmass_outsideEnergy[iss] = new TH2F( samp+ "SEL_ttmass_outsideEnergy",samp+"SEL_ttmass_outsideEnergy",200,0,600,100,0,100);
    hSEL_ttmass_insideEnergy[iss]  = new TH2F( samp+ "SEL_ttmass_insideEnergy", samp+"SEL_ttmass_insideEnergy",200,0,600,600,0,600);
    hSEL_ttmass_outsidePt[iss]     = new TH2F( samp+ "SEL_ttmass_outsidePt",    samp+"SEL_ttmass_outsidePt",200,0,600,100,0,100);
    hSEL_ngam_tmass[iss]           = new TH2F( samp+ "SEL_ngam_tmass",          samp+"SEL_ngam_tmass", 10,-0.5,9.5,100,0,5 );
    hSEL_nchg_tmass[iss]           = new TH2F( samp+ "SEL_nchg_tmass",          samp+"SEL_nchg_tmass", 10,-0.5,9.5,100,0,5 );
    hSEL_nnhad_tmass[iss]          = new TH2F( samp+ "SEL_nnhad_tmass",         samp+"SEL_nnhad_tmass", 10,-0.5,9.5,100,0,5 );
    hSEL_ngam_nchg[iss]            = new TH2F( samp+ "SEL_ngam_nchg",           samp+"SEL_ngam_nchg", 10,-0.5,9.5,10,-0.5,9.5 );

    hSEL_mcTau_costh[iss]  = new TH2F(  samp+ "SEL_MCtau_costh", samp+ "SEL_MCtau_costh", 100, 0, 1, 100, 0, 1 );

    h_dec_tmass[iss]= new TH2F( samp+"h_dec_tmass", samp+"h_dec_tmass",tauUtils::NDECAYS,-.5,tauUtils::NDECAYS-0.5,100,0,5);
    h_dec_ngam[iss] = new TH2F( samp+"h_dec_ngam" , samp+"h_dec_ngam" ,tauUtils::NDECAYS,-.5,tauUtils::NDECAYS-0.5,10,-0.5,9.5);
    h_dec_nchg[iss] = new TH2F( samp+"h_dec_nchg" , samp+"h_dec_nchg" ,tauUtils::NDECAYS,-.5,tauUtils::NDECAYS-0.5,10,-0.5,9.5);
    h_dec[iss]      = new TH1F( samp+"h_dec"      , samp+"h_dec"      ,tauUtils::NDECAYS,-.5,tauUtils::NDECAYS-0.5);
    hSEL_dec[iss]   = new TH1F( samp+"hSEL_dec"   , samp+"hSEL_dec"   ,tauUtils::NDECAYS,-.5,tauUtils::NDECAYS-0.5);

    
    for (int i=0; i<tauUtils::NDECAYS; i++) {
      h_dec_tmass[iss]->GetXaxis()->SetBinLabel( i+1, tauUtils::getTauDecLab(i) );
      h_dec_ngam[iss]->GetXaxis()->SetBinLabel( i+1, tauUtils::getTauDecLab(i) );
      h_dec_nchg[iss]->GetXaxis()->SetBinLabel( i+1, tauUtils::getTauDecLab(i) );
      h_dec[iss]->GetXaxis()->SetBinLabel( i+1, tauUtils::getTauDecLab(i) );
      hSEL_dec[iss]->GetXaxis()->SetBinLabel( i+1, tauUtils::getTauDecLab(i) );
    }

    h_a3p_cone_chgMC_pfo[iss] = new TH2F( samp+"_a3p_cone_chgMC_pfo", samp+"_a3p_cone_chgMC_pfo", 8,-0.5,7.5,8,-0.5,7.5 );
    h_a3p_cone_chgMC_trk[iss] = new TH2F( samp+"_a3p_cone_chgMC_trk", samp+"_a3p_cone_chgMC_trk", 8,-0.5,7.5,8,-0.5,7.5 );
    h_a3p_cone_trk_pfo  [iss] = new TH2F( samp+"_a3p_cone_trk_pfo",   samp+"_a3p_cone_trk_pfo",   8,-0.5,7.5,8,-0.5,7.5 );

  }




  for (int i=0; i<4; i++) {
    _nOrig[i]=0;
    _nTwoSeeds[i]=0;
    _nTwoWellMatchedSeeds[i]=0;
    _nGoodSeedDir[i]=0;
    _nSel[i]=0;
    _nSel_tmass[i]=0;
    _nSel_outofcone[i]=0;
    _nSel_angle[i]=0;
    _nSel_chg[i]=0;
    _nSel_dilep[i]=0;
  }

  return;
}

void danielKeitaTauFinderProcessor::processRunHeader( LCRunHeader* run) { 
  cout << "hello from danielKeitaTauFinderProcessor::processRunHeader" << endl;
}

void danielKeitaTauFinderProcessor::processEvent( LCEvent * evt ) { 

  cout << "processEvent " << evt->getEventNumber() << endl;


  int isample(0); // MC sample - 0:high mass tt, 1:med mass tt, 2: lowmass tt, 3: mm
  float tautauInvMass(0);
  float taucosth[2]={-1,-1};
  std::vector <MCParticle*> stableMCtaudaughters[2];
  std::vector <MCParticle*> allMCtaudaughters[2];
  std::vector <MCParticle*> finalmctaus;
  int MCdecayMode[2]={-1,-1};

  try {
    LCCollection* mccol = evt->getCollection( "MCParticle" );


    finalmctaus = tauUtils::findFinalTaus( mccol );

    if ( finalmctaus.size()!=2 ) isample=3;
    else {
      tautauInvMass = ( tauUtils::getTLV(finalmctaus[0]) + tauUtils::getTLV(finalmctaus[1]) ).M();
      if ( tautauInvMass > 480 ) isample=0;
      else if ( tautauInvMass > 75 ) isample=1;
      else isample=2;

      taucosth[0] = fabs( tauUtils::getTLV(finalmctaus[0]).Vect().CosTheta() );
      taucosth[1] = fabs( tauUtils::getTLV(finalmctaus[1]).Vect().CosTheta() );
      if ( taucosth[0] > taucosth[1] ) {
	float temp = taucosth[0];
	taucosth[0] = taucosth[1];
	taucosth[1] = temp;
      }

    }

    _nOrig[isample]++;

    if ( finalmctaus.size()==2 ) {
      for (size_t i=0; i<finalmctaus.size(); i++) {
	stableMCtaudaughters[i] = tauUtils::getstablemctauDaughters( finalmctaus[i] );
	allMCtaudaughters[i] = tauUtils::getmctauDaughters( finalmctaus[i] );
	MCdecayMode[i] = tauUtils::getMCdecayMode( stableMCtaudaughters[i] );
      }
    }

  } catch(DataNotAvailableException &e) {};


  


  if ( tautauInvMass>0 ) {
    h_ttmass->Fill(tautauInvMass);
    h_mcTau_costh[isample]->Fill(taucosth[0], taucosth[1]);
  }

  //  cout << tautauInvMass << " " << isample << endl;

  //-----------------------------

  LCRelationNavigator* relNavi(0);
  try {
    LCCollection* linkcol = evt->getCollection( "RecoMCTruthLink" );
    relNavi = new LCRelationNavigator(linkcol);
  } catch(DataNotAvailableException &e) {};


  //-----------------------------


  try {
    LCCollection* pfocol = evt->getCollection( "DistilledPFOs" );

    
    //identify which really come from tau
    std::vector < ReconstructedParticle* > PFOfromtau[2];
    std::vector < ReconstructedParticle* > noMatchPFOs;
    std::vector < ReconstructedParticle* > underlyingPFOs;

    for (int j=0; j<pfocol->getNumberOfElements(); j++) {
      ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (pfocol->getElementAt(j));

      std::vector < ReconstructedParticle* > subpfos;
      subpfos.push_back( pfo );
      if ( pfo->getParticles().size()>0 ) {
	for ( size_t kk=0; kk<pfo->getParticles().size(); kk++ ) {
	  subpfos.push_back( pfo->getParticles()[kk] );
	}
      }

      //      cout << "matching " << pfo->getType() << " " << pfo->getEnergy() << " nsub " << subpfos.size() << endl;


      std::vector < MCParticle* > mcpmatch;
      for ( size_t kk=0; kk<subpfos.size(); kk++ ) {

	//	cout << " subpfo " << kk << " charge, en " <<  subpfos[kk]->getCharge() << " " << subpfos[kk]->getEnergy() << endl;

	if ( abs( subpfos[kk]->getCharge() ) > 0.1 ) {
	  MCParticle* mcp = tauUtils::getBestTrackMatch ( subpfos[kk], relNavi );
	  if ( mcp ) mcpmatch.push_back(mcp);
	} else {
	  MCParticle* mcp = tauUtils::getBestCaloMatch ( subpfos[kk], relNavi );
	  if ( mcp ) mcpmatch.push_back(mcp);
	}
      }

      //      for ( size_t k=0 ; k<mcpmatch.size(); k++) {
      //	cout << k << " " << mcpmatch[k]->getPDG() << " " <<  mcpmatch[k]->getEnergy() << endl;
      //      }
	  

      bool taumatched(false);
      if ( mcpmatch.size()==0 ) {
	noMatchPFOs.push_back( pfo );
      } else {
	for (int itau=0; itau<2; itau++) {
	  for (int jj=0; jj<mcpmatch.size(); jj++) {
	    if ( std::find( allMCtaudaughters[itau].begin(), allMCtaudaughters[itau].end(), mcpmatch[jj] ) != allMCtaudaughters[itau].end() ) {
	      PFOfromtau[itau].push_back( pfo );
	      taumatched=true;
	      break;
	    }
	  }
	  if ( taumatched ) break;
	}
	if (!taumatched) underlyingPFOs.push_back( pfo );
      }
    }



//     for ( int itau=0; itau<2; itau++) {
//       for ( size_t j=0; j<PFOfromtau[itau].size(); j++) {
// 	cout << "PFO matched to tau " << itau << " " << j << " " << PFOfromtau[itau][j]->getType() << " " <<  PFOfromtau[itau][j]->getEnergy() << endl;
//       }
//     }
//     for ( size_t j=0; j<noMatchPFOs.size(); j++) {
//       cout << "unmatched PFO " << j << " " << noMatchPFOs[j]->getType() << " " << noMatchPFOs[j]->getEnergy() << endl;
//     }
//     for ( size_t j=0; j<underlyingPFOs.size(); j++) {
//       cout << "PFO matched to underlying event " << j << " " << underlyingPFOs[j]->getType() << " " << underlyingPFOs[j]->getEnergy() << endl;
// 
//       if ( underlyingPFOs[j]->getEnergy() > 10. ) {
// 
// 	ReconstructedParticle* pfo =  underlyingPFOs[j] ;
// 	for (size_t k=0; k<relNavi->getRelatedToObjects( pfo ).size(); k++) {
// 	  MCParticle* mcp = dynamic_cast <MCParticle*> (relNavi->getRelatedToObjects(pfo)[k]);
// 	  float wgt = relNavi->getRelatedToWeights(pfo)[k];
// 	  int pdg = mcp->getPDG();
// 	  float trackwgt = (int(wgt)%10000)/1000.;
// 	  float clusterwgt = (int(wgt)/10000)/1000. ;
// 
// 	  cout << "     " << mcp->getPDG() << " " << mcp->getEnergy() << " trk/calo wt: " << trackwgt << " " << clusterwgt << endl;
// 	  for ( size_t jj=0; jj<mcp->getParents().size(); jj++) {
// 	    MCParticle* mcpp=mcp->getParents()[jj];
// 	    cout << "       mc parent " << mcpp->getPDG() << " " << mcpp->getEnergy() << endl;
// 
// 	    for ( size_t kk=0; kk<mcpp->getParents().size(); kk++) {
// 	      MCParticle* mcppp=mcpp->getParents()[kk];
// 	      cout << "         mc grandparent " << mcppp->getPDG() << " " << mcppp->getEnergy() << endl;
// 	    }
// 
// 
// 	  }
// 
// 	}
//       }
// 
// 
//     }
// 
//     cout << "CHECK N " << PFOfromtau[0].size()+PFOfromtau[1].size()+noMatchPFOs.size()+underlyingPFOs.size() << " " << pfocol->getNumberOfElements() << endl;


    // first find the best seed prongs
    std::pair< float , ReconstructedParticle* > highestPtChargedPFO[2];
    highestPtChargedPFO[0]=std::pair< float , ReconstructedParticle* >(-1, NULL);
    highestPtChargedPFO[1]=std::pair< float , ReconstructedParticle* >(-1, NULL);

    // first highest pt track
    for (int j=0; j<pfocol->getNumberOfElements(); j++) {
      ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (pfocol->getElementAt(j));
      if ( pfo->getCharge()!=0 ) {
	float pt = sqrt( pow( pfo->getMomentum()[0], 2 ) + pow( pfo->getMomentum()[1], 2 ) );
	if ( pt > highestPtChargedPFO[0].first ) {
	  highestPtChargedPFO[0] = std::pair< float , ReconstructedParticle* > ( pt, pfo );
	} 
      }
    }



    //if ( ! highestPtChargedPFO[0].second ) {
    //  cout << "did not find first seed track!" << endl;
    //} else {
    if ( highestPtChargedPFO[0].second ) {
      // find second cone seed
      TVector3 coneDir[2];
      coneDir[0] = tauUtils::getTLV(highestPtChargedPFO[0].second).Vect();
      // then highest pt at least pi/2 away from first in deltaPhi
      for (int j=0; j<pfocol->getNumberOfElements(); j++) {
	ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (pfocol->getElementAt(j));
	TVector3 mom( pfo->getMomentum()[0],pfo->getMomentum()[1],pfo->getMomentum()[2]);
	if ( pfo->getCharge()!=0 && 
	     abs(mom.DeltaPhi(coneDir[0]))>TMath::Pi()/2. ) {
	  float pt = mom.Pt();
	  if ( pt > highestPtChargedPFO[1].first ) {
	    highestPtChargedPFO[1] = std::pair< float , ReconstructedParticle* > ( pt, pfo );
	  } 
	}
      }





      // if ( ! highestPtChargedPFO[1].second ) {
      // 	cout << "did not find second seed track!" << endl;
      // } else {
      if ( highestPtChargedPFO[1].second ) {
	Track* trk1 = highestPtChargedPFO[0].second->getTracks()[0];
	Track* trk2 = highestPtChargedPFO[1].second->getTracks()[0];
	if (trk1 && trk2 ) {
	  vertexInfo* vtxInfo = new vertexInfo();
	  vtxInfo->addTrack(trk1);
	  vtxInfo->addTrack(trk2);
	  vtxInfo->setSeedPos( TVector3(0,0,trk1->getZ0()) );
	  TVector3 vtxpos = vtxInfo->getVertexPosition();
	  cout << "ntrk, valid?, vtxChisq " << vtxInfo->getNtrack() << " " << vtxInfo->isValid() << " " << vtxInfo->getVertexChisq() << 
	    " pos " << vtxpos[0] << " " << vtxpos[1] << " " << vtxpos[2] << endl;
	  // the 3 eigenvectors of the vtx ellipse
	  TVector3 evec0 = vtxInfo->getEigenVector(0);
	  TVector3 evec1 = vtxInfo->getEigenVector(1);
	  TVector3 evec2 = vtxInfo->getEigenVector(2);
	  delete vtxInfo;
	}




	coneDir[1] = tauUtils::getTLV(highestPtChargedPFO[1].second).Vect();

	float prongangle =  TMath::Pi() - tauUtils::getTLV( highestPtChargedPFO[0].second ).Angle( tauUtils::getTLV( highestPtChargedPFO[1].second ).Vect() );

	_nTwoSeeds[isample]++;

	// try to associate seeds with MC tau directions

	MCParticle* matchedTauByDir[2]={0,0};

	if ( finalmctaus.size()==2 ) {
	  for (int jseed=0; jseed<2; jseed++) { // the two seeds
	    for (int itau=0; itau<2; itau++) {
	      float tauSeedAngle =  tauUtils::getTLV( finalmctaus[itau] ).Angle( coneDir[jseed] );
	      if ( tauSeedAngle < 0.1 ) {
		matchedTauByDir[jseed]=finalmctaus[itau];
	      }
	    }
	  }
	  assert( !matchedTauByDir[0] || matchedTauByDir[0] != matchedTauByDir[1] );
	}

	bool goodSeedDir = matchedTauByDir[0] && matchedTauByDir[1];

	if (goodSeedDir) {
	  _nGoodSeedDir[isample]++;
	}

	// check if the chosen seed tracks are really from tau
	int ttre[2]={-1,-1};
	for (int jseed=0; jseed<2; jseed++) { // the two seeds
	  ReconstructedParticle* pfo = highestPtChargedPFO[jseed].second;
	  MCParticle* bestmatch = tauUtils::getBestTrackMatch( pfo, relNavi );
	  if ( !bestmatch ) {
	    cout << "could not find MCP match..." << pfo->getType() << " " << pfo->getEnergy() << " " << relNavi->getRelatedToWeights(pfo).size() << endl;
	  } else {
	    // which tau does this MC particle come from?
	    for (int jtau=0; jtau<2; jtau++) {
	      if ( find( allMCtaudaughters[jtau].begin(), allMCtaudaughters[jtau].end(), bestmatch ) != allMCtaudaughters[jtau].end() ) {
		ttre[jseed]=jtau;
	      }
	    }
	  }
	}
	bool goodMCMatchSeeds = ttre[0]>=0 && ttre[1]>=0 && ttre[0]!=ttre[1];
	if (goodMCMatchSeeds) {
	  _nTwoWellMatchedSeeds[isample]++;
	}
	// then gather particles within cone around seeds

	std::vector < ReconstructedParticle* > coneparticles[2];
	std::vector < ReconstructedParticle* > outsideconeparticles;
	std::vector < ReconstructedParticle* > outsideconephotons;
	for (int j=0; j<pfocol->getNumberOfElements(); j++) {
	  ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (pfocol->getElementAt(j));
	  TLorentzVector ttr = tauUtils::getTLV(pfo);
	  float angle0=ttr.Vect().Angle( coneDir[0] );
	  float angle1=ttr.Vect().Angle( coneDir[1] );
	  if ( angle0<_conesize ) {
	    coneparticles[0].push_back( pfo );
	  } else if ( angle1<_conesize ) {
	    coneparticles[1].push_back( pfo );
	  } else {
	    // cout << "outdside cone: " << pfo->getType() << " " << pfo->getEnergy() << " " << ttr.CosTheta() << " " << angle0 << " " << angle1 << endl;
	    outsideconeparticles.push_back( pfo );
	    if ( pfo->getType()==22 ) 
	      outsideconephotons.push_back( pfo );
	  }
	}

	TLorentzVector jet4mom[2];
	int cone_chg[2]={0,0};
	int cone_ngam[2]={0,0};
	int cone_nchg[2]={0,0};
	int cone_nnhad[2]={0,0};

	for (int j=0; j<2; j++) {

	  if ( isample==0 ) { // extra checks for multi-prongs

	    int decay=-1;
	    //cout << " looking at cone for tau " << j << " matched to mctau " << matchedTauByDir[j];
	    if ( matchedTauByDir[j] ) {
	      decay =  tauUtils::getMCdecayMode(  matchedTauByDir[j]  );
	      //cout << " real tau decay " << decay ;
	    }
	    //cout << endl;

	    if ( decay == tauUtils::decayA1_3p ) {
	      TVector3 tauDir =  tauUtils::getTLV( matchedTauByDir[j] ).Vect();
	      int nmcInCone(0);
	      //cout << "MC a3p decay!" << endl;
	      //cout << " mc stable daughters" << endl;
	      std::vector < EVENT::MCParticle* > stmc = tauUtils::getstablemctauDaughters( matchedTauByDir[j] );
	      for (size_t k=0; k<stmc.size(); k++) {
		float angle=tauUtils::getTLV(  stmc[k] ).Angle( tauDir );
		if ( stmc[k]->getCharge()!=0 && angle<0.1 ) nmcInCone++;
		//cout << stmc[k]->getPDG() << " " << stmc[k]->getEnergy() << " " << angle << endl;
	      }
	      //	      cout << " cone menbers:  " << endl;
	      //for ( size_t i=0; i<coneparticles[j].size(); i++) {
	      //	cout << coneparticles[j][i]->getType() << " " << coneparticles[j][i]->getEnergy() << " " << tauUtils::getTLV( coneparticles[j][i] ).Angle( tauDir ) << endl;
	      // }
	      int npfoInCone(0);
	      //cout << "all charged pfos" << endl;
	      for (int jj=0; jj<pfocol->getNumberOfElements(); jj++) {
		ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*> (pfocol->getElementAt(jj));
		if (pfo->getCharge()!=0) {
		  float angle = tauUtils::getTLV( pfo ).Angle( tauDir );
		  if ( angle<0.1 ) npfoInCone++;
		  //  cout << pfo->getType() << " " << pfo->getEnergy() << " " <<  angle << endl;
		}
	      }
	      int nTrkInCone(0);
	      try {
		LCCollection* trkcol = evt->getCollection( "MarlinTrkTracks" );
		for (int jj=0; jj<trkcol->getNumberOfElements(); jj++) {
		  Track* trk = dynamic_cast<Track*> (trkcol->getElementAt(jj));
		  TLorentzVector trktlv = tauUtils::getFourMomentum( trk->getTrackState( TrackState::AtIP ), 0.139 );
		  float angle=trktlv.Angle( tauDir );
		  if ( angle < 0.1 ) nTrkInCone++;
		  //	  cout << "marlin trk " << jj << " " << trktlv.Energy() << " " << angle << endl;
		}
	      } catch(DataNotAvailableException &e) {};

	      //	      cout << "ch cone contents: mc,pfo,trk = " << nmcInCone << " " << npfoInCone << " " << nTrkInCone << endl;
	      h_a3p_cone_chgMC_pfo[isample]->Fill(nmcInCone,npfoInCone);
	      h_a3p_cone_chgMC_trk[isample]->Fill(nmcInCone,nTrkInCone);
	      h_a3p_cone_trk_pfo  [isample]->Fill(nTrkInCone,npfoInCone);
	    }
	  }

	  jet4mom[j].SetXYZT(0,0,0,0);
	  for ( size_t i=0; i<coneparticles[j].size(); i++) {
	    TLorentzVector thisTLV = tauUtils::getTLV( coneparticles[j][i] );
	    jet4mom[j] += thisTLV;
	    cone_chg[j]+=coneparticles[j][i]->getCharge();
	    if ( coneparticles[j][i]->getType()==22 )       cone_ngam[j]++;
	    else if ( coneparticles[j][i]->getType()==111 ) cone_ngam[j]+=2;
	    else if ( coneparticles[j][i]->getCharge()!=0 ) cone_nchg[j]++;
	    if  ( coneparticles[j][i]->getType()==2112 ) cone_nnhad[j]++;
	  }
	}



	// should look for singleton photons outside cone, which make pi0 mass with single photon inside cone and dont take jet mass above mtau
	for ( size_t i=0; i<outsideconephotons.size(); i++) {
	  TLorentzVector outsidePhotonTlv = tauUtils::getTLV( outsideconephotons[i] );
	  for ( int j=0; j<2; j++) {
	    for ( size_t k=0; k<coneparticles[j].size(); k++) {
	      if ( coneparticles[j][k]->getType()==22 ) {
		TLorentzVector insidePhotonTlv = tauUtils::getTLV( coneparticles[j][k] );
		h_gamgam_mass->Fill( (outsidePhotonTlv+insidePhotonTlv).M() );
	      }
	    }
	  }
	} // photons outside cone
	// this seems not so useful (most pi0 already reconstructed in distilled collection?)



	// check energy outside cone
	float outsideConeEnergy(0);
	float outsideConePt(0);
	for (size_t i=0; i<outsideconeparticles.size(); i++) {
	  TLorentzVector outTlv = tauUtils::getTLV( outsideconeparticles[i] );
	  outsideConeEnergy+=outTlv.E();
	  outsideConePt+=outTlv.Pt();
	}


	// now lets see how well we did by comparing to MC
	float taumatchedEnergy(0);
	float nontaumatchedEnergy(0);
	float unmatchedEnergy(0);
	for ( int j=0; j<2; j++) { // the 2 cones
	  for ( size_t k=0; k<coneparticles[j].size(); k++) {
	    ReconstructedParticle* rp = coneparticles[j][k];
	    int mtau(-1);
	    for (int itau=0; itau<2; itau++) {
	      if ( find( PFOfromtau[itau].begin(), PFOfromtau[itau].end(), rp )!=PFOfromtau[itau].end() ) {
		mtau=itau;
		break;
	      }
	    }
	    if ( mtau>=0 ) {
	      taumatchedEnergy+=rp->getEnergy();
	    } else {
	      if ( find ( noMatchPFOs.begin(), noMatchPFOs.end(), rp ) != noMatchPFOs.end() ) {
		unmatchedEnergy+=rp->getEnergy(); // not matched to a MC particle
	      } else {
		nontaumatchedEnergy+=rp->getEnergy(); // matched to non-tau MC particle
	      }
	    }
	  }
	}
	//	cout << "MC tau-matched / non-tau / unmatched energy = " << taumatchedEnergy << " / " << nontaumatchedEnergy << " / " << unmatchedEnergy << endl;


	
	bool mumu(true);
	for (int icone=0; icone<2; icone++) {
	  int nch(0);
	  int nmu(0);
	  int ngam(0);
	  for ( size_t ip=0; ip<coneparticles[icone].size(); ip++) {
	    ReconstructedParticle* rp = coneparticles[icone][ip];
	    if ( rp->getCharge()!=0 ) {
	      nch++;
	      if ( abs(rp->getType())==13 ) nmu++;
	    } else {
	      if ( rp->getType()==22 ) {
		ngam++;
	      }
	    }
	  }
	  if ( nch>1 || nmu!=1 || ngam>1 ) {
	    mumu=false;
	    break;
	  }
	}
		  
	bool elel(true);
	for (int icone=0; icone<2; icone++) {
	  int nch(0);
	  int nel(0);
	  int ngam(0);
	  for ( size_t ip=0; ip<coneparticles[icone].size(); ip++) {
	    ReconstructedParticle* rp = coneparticles[icone][ip];
	    if ( rp->getCharge()!=0 ) {
	      nch++;
	      if ( abs(rp->getType())==11 ) nel++;
	    } else {
	      if ( rp->getType()==22 ) {
		ngam++;
	      }
	    }
	  }
	  if ( nch>1 || nel!=1 || ngam>1 ) {
	    elel=false;
	    break;
	  }
	}
		  

	float insideConeEnergy = jet4mom[0].E() + jet4mom[1].E() ;
		   
	float maxMass = 2.5;

	bool tauMassSel = jet4mom[0].M()<maxMass && jet4mom[1].M()<maxMass;
	bool prongAngleCut = prongangle < 0.15;
	bool outsideEnergyCut = outsideConeEnergy < 40. && outsideConePt < 20.;
	bool chargeCut = cone_chg[0]*cone_chg[1]==-1;

	bool dilepCut = ! ( mumu && insideConeEnergy>450 ) && ! ( elel  && insideConeEnergy>450 );

	bool select = tauMassSel&&prongAngleCut&&outsideEnergyCut&&chargeCut&&dilepCut;

	



	if (tauMassSel)        _nSel_tmass[isample]++;
	if (outsideEnergyCut)  _nSel_outofcone[isample]++;
	if (prongAngleCut)     _nSel_angle[isample]++;
	if (chargeCut)         _nSel_chg[isample]++;
	if (dilepCut)          _nSel_dilep[isample]++;
	if (select)            _nSel[isample]++;

	
	if ( tautauInvMass>0 && select ) {
	  hSEL_ttmass->Fill(tautauInvMass);
	  hSEL_mcTau_costh[isample]->Fill(taucosth[0], taucosth[1]);
	}

	h_ttmass_prongAngle[isample]->Fill( tautauInvMass, prongangle );
	h_ttmass_outsideEnergy[isample]->Fill( tautauInvMass, outsideConeEnergy );
	h_ttmass_insideEnergy[isample]->Fill( tautauInvMass, insideConeEnergy );
	h_ttmass_outsidePt[isample]->Fill( tautauInvMass, outsideConePt );

	for (int j=0; j<2; j++) {	    // j = cone index

	  if ( finalmctaus.size()==2 && matchedTauByDir[0] && matchedTauByDir[1] ) {

	    int mctauIndex(-1);
	    if ( matchedTauByDir[j] == finalmctaus[0] ) mctauIndex=0;
	    else if ( matchedTauByDir[j] == finalmctaus[1] ) mctauIndex=1;
	    else {
	      cout << "unmatched tau by dir! " << endl;
	      assert(0);
	    }


	    //if (  isample==0 && MCdecayMode[mctauIndex] ==  tauUtils::decayA1_3p ) {
	    //  cout << " FFFFF " << j << " " <<  cone_nchg[j]  << endl;
	    // }

	    h_dec_tmass [isample]-> Fill ( MCdecayMode[mctauIndex], jet4mom[j].M() );
	    h_dec_ngam  [isample]-> Fill ( MCdecayMode[mctauIndex], cone_ngam[j] );
	    h_dec_nchg  [isample]-> Fill ( MCdecayMode[mctauIndex], cone_nchg[j] );
	    h_dec       [isample]-> Fill ( MCdecayMode[mctauIndex] );
	  }

	  h_tchg_tmass[isample]->Fill( cone_chg[j], jet4mom[j].M() );
	  h_ngam_tmass[isample]->Fill( cone_ngam[j], jet4mom[j].M() );
	  h_nchg_tmass[isample]->Fill( cone_nchg[j], jet4mom[j].M() );
	  h_nnhad_tmass[isample]->Fill( cone_nnhad[j], jet4mom[j].M() );
	  h_ngam_nchg[isample]->Fill ( cone_ngam[j], cone_nchg[j] );
	}
	if ( select ) {
	  hSEL_ttmass_prongAngle[isample]->Fill( tautauInvMass, prongangle );
	  hSEL_ttmass_outsideEnergy[isample]->Fill( tautauInvMass, outsideConeEnergy );
	  hSEL_ttmass_insideEnergy[isample]->Fill( tautauInvMass, jet4mom[0].E() + jet4mom[1].E() );
	  hSEL_ttmass_outsidePt[isample]->Fill( tautauInvMass, outsideConePt );
	  for (int j=0; j<2; j++) {
	    hSEL_dec[isample]    -> Fill ( MCdecayMode[j] );
	    hSEL_tchg_tmass[isample]->Fill( cone_chg[j], jet4mom[j].M() );
	    hSEL_ngam_tmass[isample]->Fill( cone_ngam[j], jet4mom[j].M() );
	    hSEL_nchg_tmass[isample]->Fill( cone_nchg[j], jet4mom[j].M() );
	    hSEL_nnhad_tmass[isample]->Fill( cone_nnhad[j], jet4mom[j].M() );
	    hSEL_ngam_nchg[isample]->Fill(  cone_ngam[j], cone_nchg[j] );
	  }
	}


	if ( select ) {

	  // write out the taus
	  LCCollectionVec* newcol = new LCCollectionVec();
	  for ( int j=0; j<2; j++) {

	    ReconstructedParticleImpl* tau = new ReconstructedParticleImpl();
	    for ( size_t k=0; k<coneparticles[j].size(); k++) {
	      tau->addParticle( coneparticles[j][k] );
	    }
	    tau->setType(15);
	    TLorentzVector mom = tauUtils::getTLV( coneparticles[j] );
	    float threemom[3];
	    threemom[0]=mom.X();
	    threemom[1]=mom.Y();
	    threemom[2]=mom.Z();
	    tau->setMomentum(threemom);
	    tau->setEnergy( mom.E() );
	    tau->setMass( mom.M() );
	    newcol->addElement( tau );

	  }
	  
	  evt->addCollection( newcol, "danielKeitaRecoTaus" );

	}

      } // second seeds found
    } // first seed

  } catch(DataNotAvailableException &e) {};


  return;
}


void danielKeitaTauFinderProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void danielKeitaTauFinderProcessor::end(){

  cout << "end" << endl;

  _fout->cd();

  TH1F* eff_ttmass = (TH1F*) h_ttmass -> Clone("seleff_ttmass");
  eff_ttmass->Sumw2();
  eff_ttmass->Divide( hSEL_ttmass, h_ttmass, 1, 1, "B" );
  
  for (int i=0; i<4; i++) {
    TH1D* hh1 = h_mcTau_costh[i]->ProjectionX();
    TH1D* hh2 = hSEL_mcTau_costh[i]->ProjectionX();
    TH1D* heff = (TH1D*) hh1->Clone( "seleff_"+TString(hh1->GetName()) );
    heff->Sumw2();
    heff->Divide( hh2, hh1, 1, 1, "B" );

    TH1F* eff_dec = (TH1F*) h_dec[i] -> Clone("seleff_"+TString(h_dec[i]->GetName()) );
    eff_dec->Sumw2();
    eff_dec->Divide( hSEL_dec[i], h_dec[i], 1, 1, "B" );

  }

  _fout->Write(0);
  _fout->Close();



  cout << " _nOrig          "; for (int i=0; i<4; i++) cout << _nOrig[i] << " "; cout << endl;
  cout << " _nTwoSeeds      "; for (int i=0; i<4; i++) cout << _nTwoSeeds[i] << " "; cout << endl;
  cout << " _nTwoSeedsEff   "; for (int i=0; i<4; i++) cout << 1.*_nTwoSeeds[i]/_nOrig[i] << " "; cout << endl;
  cout << " nGoodSeedDir    "; for (int i=0; i<4; i++) cout << _nGoodSeedDir[i] << " "; cout << endl;
  cout << " nGoodSeedDirEff "; for (int i=0; i<4; i++) cout << 1.*_nGoodSeedDir[i]/_nOrig[i] << " "; cout << endl;
  cout << " _nSel_tmass     "; for (int i=0; i<4; i++) cout << _nSel_tmass[i] << " "; cout << endl;
  cout << " _nSel_tmassEff  "; for (int i=0; i<4; i++) cout << 1.*_nSel_tmass[i]/_nOrig[i] << " "; cout << endl;
  cout << " _nSel_oocone    "; for (int i=0; i<4; i++) cout << _nSel_outofcone[i] << " "; cout << endl;
  cout << " _nSel_ooconeEff "; for (int i=0; i<4; i++) cout << 1.*_nSel_outofcone[i]/_nOrig[i] << " "; cout << endl;
  cout << " _nSel_angle     "; for (int i=0; i<4; i++) cout << _nSel_angle[i] << " "; cout << endl;
  cout << " _nSel_angleEff  "; for (int i=0; i<4; i++) cout << 1.*_nSel_angle[i]/_nOrig[i] << " "; cout << endl;
  cout << " _nSel_chg       "; for (int i=0; i<4; i++) cout << _nSel_chg[i] << " "; cout << endl;
  cout << " _nSel_chgEff    "; for (int i=0; i<4; i++) cout << 1.*_nSel_chg[i]/_nOrig[i] << " "; cout << endl;
  cout << " _nSel_dilep     "; for (int i=0; i<4; i++) cout << _nSel_dilep[i] << " "; cout << endl;
  cout << " _nSel_dilepEff  "; for (int i=0; i<4; i++) cout << 1.*_nSel_dilep[i]/_nOrig[i] << " "; cout << endl;
  cout << " _nSel           "; for (int i=0; i<4; i++) cout << _nSel[i] << " "; cout << endl;
  cout << " _nSelEff        "; for (int i=0; i<4; i++) cout << 1.*_nSel[i]/_nOrig[i] << " "; cout << endl;

//  cout << " n2mctau            " << _nTwoMcTau << endl;
//  cout << " nHighMass          " << _nTwoMcTauHighMass << endl;
//  cout << " n2seeds            " << _nTwoSeeds[0] << " " << 1.*_nTwoSeeds[0]/_nTwoMcTau << endl;
//  cout << " n2seeds HM         " << _nTwoSeeds[1] << " " << 1.*_nTwoSeeds[1]/_nTwoMcTauHighMass << endl;
//  cout << " n2mcgoodseeds      " << _nTwoWellMatchedSeeds[0] << " " << 1.*_nTwoWellMatchedSeeds[0]/_nTwoMcTau << endl;
//  cout << " n2mcgoodseeds  HM  " << _nTwoWellMatchedSeeds[1] << " " << 1.*_nTwoWellMatchedSeeds[1]/_nTwoMcTauHighMass << endl;
//
//  cout << " nGoodSeedDir    " << _nGoodSeedDir[0] << " " << 1.*_nGoodSeedDir[0]/_nTwoMcTau << endl;
//  cout << " nGoodSeedDir HM " << _nGoodSeedDir[1] << " " << 1.*_nGoodSeedDir[1]/_nTwoMcTauHighMass << endl;
//
//  cout << " nseltmass        " <<   _nSel_tmass    [0] << " " << 1.*_nSel_tmass    [0]/_nTwoMcTau << endl;
//  cout << " nseltmass     HM " <<   _nSel_tmass    [1] << " " << 1.*_nSel_tmass    [1]/_nTwoMcTauHighMass << endl;
//  cout << " nseloutofcone    " <<   _nSel_outofcone[0] << " " << 1.*_nSel_outofcone[0]/_nTwoMcTau << endl;
//  cout << " nseloutofcone HM " <<   _nSel_outofcone[1] << " " << 1.*_nSel_outofcone[1]/_nTwoMcTauHighMass << endl;
//  cout << " nselangle        " <<   _nSel_angle    [0] << " " << 1.*_nSel_angle    [0]/_nTwoMcTau << endl;
//  cout << " nselangle     HM " <<   _nSel_angle    [1] << " " << 1.*_nSel_angle    [1]/_nTwoMcTauHighMass << endl;
//  cout << " nselchg          " <<   _nSel_chg      [0] << " " << 1.*_nSel_chg      [0]/_nTwoMcTau << endl;     
//  cout << " nselchg       HM " <<   _nSel_chg      [1] << " " << 1.*_nSel_chg      [1]/_nTwoMcTauHighMass << endl;     
//
//  cout << " nsel               " << _nSel[0] << " " << 1.*_nSel[0]/_nTwoMcTau << endl;
//  cout << " nsel HM            " << _nSel[1] << " " << 1.*_nSel[1]/_nTwoMcTauHighMass << endl;


  std::cout << "danielKeitaTauFinderProcessor::end()  " << name() 
 	    << std::endl ;
}

