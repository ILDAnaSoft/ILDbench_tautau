#include "analyseTau.h"
#include <iostream>
#include <map>

#include <marlin/Global.h>
#include "lcio.h"
#include "tauUtils.h"

#include "EVENT/ReconstructedParticle.h"
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"

using namespace lcio ;
using namespace marlin ;
using std::cout;
using std::endl;

analyseTauProcessor aanalyseTauProcessor ;


analyseTauProcessor::analyseTauProcessor() : Processor("analyseTauProcessor") {
  
  // modify processor description
  _description = "analyseTauProcessor does whatever it does ..." ;



  registerProcessorParameter("outputFilename",
                             "name of output file",
                             _outfile,
                             std::string( "analyseTau.root") );



}


void analyseTauProcessor::init() { 
  cout << "hello from analyseTauProcessor::init" << endl;


  _fout = new TFile(_outfile.c_str(),"recreate");

  for (int i=0; i<3; i++) {
    TString hn = "eonp_ecalFrac_";
    switch (i) {
    case 0:
      hn+="typeE"; break;
    case 1:
      hn+="typeM"; break;
    case 2:
      hn+="typeP"; break;
    default:
      break;
    }
    
    _h_eonp_ecalFrac[i] = new TH2F(hn, hn, 100, 0, 1.1, 100, 0, 1.1);
  }


  for (int i=0; i<4; i++) {
    TString hn = "samp"; hn+=i; hn+="_mcRecDecay";
    _h_mc_rec_decMode[i] = new TH2F( hn, hn, 12, -1.5, 10.5, 12, -1.5, 10.5 );
    _h_mc_rec_decMode[i]->GetXaxis()->SetTitle("MC decay mode");
    _h_mc_rec_decMode[i]->GetYaxis()->SetTitle("RECO decay mode");

    _h_mc_rec_decMode[i]->GetXaxis()->SetBinLabel( 1, "unknown" );
    _h_mc_rec_decMode[i]->GetYaxis()->SetBinLabel( 1, "unknown" );
    for (int j=0; j<tauUtils::NDECAYS; j++) {
      _h_mc_rec_decMode[i]->GetXaxis()->SetBinLabel( j+2, tauUtils::getTauDecLab(j) );
      _h_mc_rec_decMode[i]->GetYaxis()->SetBinLabel( j+2, tauUtils::getTauDecLab(j) );
    }

    hn = "samp"; hn+=i; hn+="_mcRecDecayEff";
    _h_mc_rec_decModeEff[i] = new TH2F( hn, hn, 12, -1.5, 10.5, 12, -1.5, 10.5 );
    _h_mc_rec_decModeEff[i]->GetXaxis()->SetTitle("MC decay mode");
    _h_mc_rec_decModeEff[i]->GetYaxis()->SetTitle("RECO decay mode");

    _h_mc_rec_decModeEff[i]->GetXaxis()->SetBinLabel( 1, "unknown" );
    _h_mc_rec_decModeEff[i]->GetYaxis()->SetBinLabel( 1, "unknown" );
    for (int j=0; j<tauUtils::NDECAYS; j++) {
      _h_mc_rec_decModeEff[i]->GetXaxis()->SetBinLabel( j+2, tauUtils::getTauDecLab(j) );
      _h_mc_rec_decModeEff[i]->GetYaxis()->SetBinLabel( j+2, tauUtils::getTauDecLab(j) );
    }

    hn = "samp"; hn+=i; hn+="_mcdec_jetMass";
    _h_mcdec_jetMass[i] = new TH2F( hn, hn, 12, -1.5, 10.5, 100,0,3 );
    hn = "samp"; hn+=i; hn+="_mcdec_nChg";
    _h_mcdec_nChg[i] = new TH2F( hn, hn, 12, -1.5, 10.5, 10,-0.5,9.5 );
    hn = "samp"; hn+=i; hn+="_mcdec_nTrk";
    _h_mcdec_nTrk[i] = new TH2F( hn, hn, 12, -1.5, 10.5, 10,-0.5,9.5 );
    hn = "samp"; hn+=i; hn+="_mcdec_nGam";
    _h_mcdec_nGam[i] = new TH2F( hn, hn, 12, -1.5, 10.5, 10,-0.5,9.5 );
    hn = "samp"; hn+=i; hn+="_mcdec_nPi0";
    _h_mcdec_nPi0[i] = new TH2F( hn, hn, 12, -1.5, 10.5, 10,-0.5,9.5 );
    hn = "samp"; hn+=i; hn+="_mcdec_nPi0Gam";
    _h_mcdec_nPi0Gam[i] = new TH2F( hn, hn, 12, -1.5, 10.5, 10,-0.5,9.5 );

    _h_mcdec_jetMass[i]->GetXaxis()->SetTitle("MC decay mode");
    _h_mcdec_nChg[i]   ->GetXaxis()->SetTitle("MC decay mode");
    _h_mcdec_nTrk[i]   ->GetXaxis()->SetTitle("MC decay mode");
    _h_mcdec_nGam[i]   ->GetXaxis()->SetTitle("MC decay mode");
    _h_mcdec_nPi0[i]   ->GetXaxis()->SetTitle("MC decay mode");
    _h_mcdec_nPi0Gam[i]->GetXaxis()->SetTitle("MC decay mode");

    _h_mcdec_jetMass[i]->GetXaxis()->SetBinLabel( 1, "unknown" );
    _h_mcdec_nChg[i]   ->GetXaxis()->SetBinLabel( 1, "unknown" );
    _h_mcdec_nTrk[i]   ->GetXaxis()->SetBinLabel( 1, "unknown" );
    _h_mcdec_nGam[i]   ->GetXaxis()->SetBinLabel( 1, "unknown" );
    _h_mcdec_nPi0[i]   ->GetXaxis()->SetBinLabel( 1, "unknown" );
    _h_mcdec_nPi0Gam[i]->GetXaxis()->SetBinLabel( 1, "unknown" );

    for (int j=0; j<tauUtils::NDECAYS; j++) {
      _h_mcdec_jetMass[i]->GetXaxis()->SetBinLabel( j+2, tauUtils::getTauDecLab(j) );
      _h_mcdec_nChg[i]   ->GetXaxis()->SetBinLabel( j+2, tauUtils::getTauDecLab(j) );
      _h_mcdec_nTrk[i]   ->GetXaxis()->SetBinLabel( j+2, tauUtils::getTauDecLab(j) );
      _h_mcdec_nGam[i]   ->GetXaxis()->SetBinLabel( j+2, tauUtils::getTauDecLab(j) );
      _h_mcdec_nPi0[i]   ->GetXaxis()->SetBinLabel( j+2, tauUtils::getTauDecLab(j) );
      _h_mcdec_nPi0Gam[i]->GetXaxis()->SetBinLabel( j+2, tauUtils::getTauDecLab(j) );
    }

  }


  return;
}

void analyseTauProcessor::processRunHeader( LCRunHeader* run) { 
  cout << "hello from analyseTauProcessor::processRunHeader" << endl;
}

void analyseTauProcessor::processEvent( LCEvent * evt ) { 

  std::vector < MCParticle* > mctaus;
  try {
    LCCollection* mccol = evt->getCollection( "MCParticle" );
    mctaus = tauUtils::findFinalTaus( mccol );
  } catch(DataNotAvailableException &e) {};

  //  cout << mctaus.size() ;
  //  cout << " mc taus; decays: ";
  TLorentzVector tlv_alltausmc(0,0,0,0);
  TLorentzVector tlv_tausmc[2];
  int mcdecmode[2]={-1,-1};
  for (size_t i=0; i<mctaus.size(); i++) {
    mcdecmode[i] = tauUtils::getMCdecayMode( mctaus[i] );
    tlv_tausmc[i]=tauUtils::getTLV( mctaus[i] );
    tlv_alltausmc+=tlv_tausmc[i];
  }
  //  cout << " inv mass " << tlv_alltausmc.M();
  //  cout << endl;

  int iSample(0);
  if (mctaus.size()==2 ) {
    if      ( tlv_alltausmc.M() > 480 ) iSample=0;
    else if ( tlv_alltausmc.M() > 80  ) iSample=1;
    else                                iSample=2;
  } else {
    iSample=3;
  }

  int reco_imcMatch[2]={-1,-1};

//  PIDHandler* pidHandler;
//  try {
//    LCCollection* pfocol = evt->getCollection( "PandoraPFOs");
//    pidHandler=new PIDHandler(pfocol);
//  } catch(DataNotAvailableException &e) {};



  try {
    LCCollection* pfocol = evt->getCollection( "danielKeitaRecoTaus" );

    //    cout << pfocol->getNumberOfElements() << " reco tau jets selected " << endl;

    if ( pfocol->getNumberOfElements()==2 ) {

      int reco_dec_mode[2]={-1,-1};

      TLorentzVector visibleTLV[2];
      
      for ( int itau=0; itau<2; itau++ ) {
	visibleTLV[itau].SetXYZT(0,0,0,0);

	ReconstructedParticle* rptau = dynamic_cast <ReconstructedParticle*> (pfocol->getElementAt(itau));
	int nch(0);
	int nneu(0);
	int nchHad(0);
	int npi0(0);
	int nga(0);
	int nel(0);
	int nmu(0);
	for ( size_t ipart = 0; ipart<rptau->getParticles().size(); ipart++ ) {
	  ReconstructedParticle* rp = rptau->getParticles()[ipart];

	  visibleTLV[itau]+=tauUtils::getTLV( rp );

	  int type = rp->getType();
	  if ( rp->getCharge()==0 ) {
	    nneu++;
	    if ( type==22 ) nga++;
	    if ( type==111 ) npi0++;
	  } else {
	    nch++;
	    // check if type is reasonable
	    float mom(0);
	    for (int i=0; i<3; i++) mom+=pow( rp->getMomentum()[i], 2 );
	    mom=sqrt(mom);
	    float subens[6]={0};
	    float toten(0);
	    for (size_t icl = 0; icl<rp->getClusters().size(); icl++) {
	      Cluster* cl = rp->getClusters()[icl];
	      // parameter ClusterSubdetectorNames [string]: ecal, hcal, yoke, lcal, lhcal, bcal, 
	      for (int isub=0; isub<6; isub++) {
		subens[isub]+=cl->getSubdetectorEnergies()[isub];
		toten+=cl->getSubdetectorEnergies()[isub];
	      }
	    }
	    float ecalFrac = subens[0]/toten;
	    float EonP = toten/mom;
	    int itype(-1);
	    if (abs(type)==11 ) itype=0;
	    else if (abs(type)==13 ) itype=1;
	    else if (abs(type)==211 ) itype=2;
	    if ( itype>=0 ) _h_eonp_ecalFrac[itype]->Fill( EonP, ecalFrac );

	    if (abs(type)==11 ) nel++;
	    else if (abs(type)==13 ) nmu++;
	    else if (abs(type)==211 ) nchHad++;

	  }
	}



	int nTrkInCone(0);
	try {
	  LCCollection* trkcol = evt->getCollection( "MarlinTrkTracks" );
	  for (int jj=0; jj<trkcol->getNumberOfElements(); jj++) {
	    Track* trk = dynamic_cast<Track*> (trkcol->getElementAt(jj));

	    std::vector < int > nhits = trk->getSubdetectorHitNumbers();
	    //cout << "trk hits ";
	    //for (size_t k=0; k<nhits.size(); k++) cout << nhits[k] << " " ;
	    //cout << endl;

	    bool isGood = nhits[0]>4 &&  // VXD hits?
	      nhits[6]>20;                // TPC hits?

	    if (isGood) {
	      TLorentzVector trktlv = tauUtils::getFourMomentum( trk->getTrackState( TrackState::AtIP ), 0.139 );
	      float angle=trktlv.Angle(  visibleTLV[itau].Vect() );
	      if ( angle < 0.1 ) nTrkInCone++;
	    }
	  }
	} catch(DataNotAvailableException &e) {};
	


	//	enum { decayChPi=0, decayRho, decayA1_3p, decayA1_1p , decayEl, decayMu , decayW , decayK , decayMultiprong , decayOthersingleprong, decayUnrecognised, NDECAYS};

	// find a consistent MC tau
	for (int imc=0; imc<2; imc++) {
	  if ( tlv_tausmc[imc].Angle( visibleTLV[itau].Vect() ) < 0.1 )
	    reco_imcMatch[itau] = imc;
	}
	int imcmode(-1);
	if ( reco_imcMatch[itau] >=0 ) imcmode =  mcdecmode[reco_imcMatch[itau]] ;

	float jetmass = visibleTLV[itau].M();

	_h_mcdec_jetMass[iSample]->Fill(imcmode, jetmass);
	_h_mcdec_nChg   [iSample]->Fill(imcmode, nch);
	_h_mcdec_nTrk   [iSample]->Fill(imcmode, nTrkInCone);
	_h_mcdec_nGam   [iSample]->Fill(imcmode, nga);
	_h_mcdec_nPi0   [iSample]->Fill(imcmode, npi0);
	_h_mcdec_nPi0Gam[iSample]->Fill(imcmode, npi0+2*nga);

	if ( nch==1 ) {
	  if      ( nel==1 && jetmass<1 )   reco_dec_mode[itau]=tauUtils::decayEl;
	  else if ( nmu==1 && jetmass<0.4 ) reco_dec_mode[itau]=tauUtils::decayMu;
	  else if ( npi0==1 && nga==0 )     reco_dec_mode[itau]=tauUtils::decayRho;
	  else if ( npi0==0 && nga==1 ) {
	    if ( jetmass > 0.5 )
	      reco_dec_mode[itau]=tauUtils::decayRho;
	    else 
	      reco_dec_mode[itau]=tauUtils::decayChPi;
	  }
	  else if ( npi0==0 && nga==0 ) reco_dec_mode[itau]=tauUtils::decayChPi;
	  else if ( npi0+2*nga>2 ) {
	    if ( jetmass > 1.0 )
	      reco_dec_mode[itau]=tauUtils::decayA1_1p;
	    else reco_dec_mode[itau]=tauUtils::decayRho;
	  }
	  else                          reco_dec_mode[itau]=tauUtils::decayOthersingleprong;
	} else if ( nch==3 ) {
	  if    (  npi0==0 && nga<2  ) reco_dec_mode[itau]=tauUtils::decayA1_3p;
	  else                         reco_dec_mode[itau]=tauUtils::decayMultiprong;
	} else {
	  reco_dec_mode[itau]=tauUtils::decayMultiprong;
	}

	_h_mc_rec_decMode[iSample]->Fill( imcmode,  reco_dec_mode[itau] );
	_h_mc_rec_decModeEff[iSample]->Fill( imcmode,  reco_dec_mode[itau] );


      } // the 2 reco taus




      //cout << "reco tau decays: " << reco_dec_mode[0] << " " << reco_dec_mode[1] << endl;
      //cout << " masses " << visibleTLV[0].M() << " " << visibleTLV[1].M() << endl;

      if ( ( reco_dec_mode[0] == tauUtils::decayRho || reco_dec_mode[0] == tauUtils::decayChPi ) && 
	   ( reco_dec_mode[1] == tauUtils::decayRho || reco_dec_mode[1] == tauUtils::decayChPi ) ) {
	
	//cout << " both pi or rho reco " << endl;
	
	
      }

    }
  } catch(DataNotAvailableException &e) {};

  return;
}



void analyseTauProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void analyseTauProcessor::end(){

  // rescale the eff histo
  for (int i=0; i<4; i++) {
    for (int ib=1; ib<=_h_mc_rec_decModeEff[i]->GetNbinsX(); ib++) {
      float totcol(0);
      for (int jb=1; jb<=_h_mc_rec_decModeEff[i]->GetNbinsY(); jb++) {
	totcol+=_h_mc_rec_decModeEff[i]->GetBinContent(ib,jb);
      }
      if ( totcol>0 ) {
	for (int jb=1; jb<=_h_mc_rec_decModeEff[i]->GetNbinsY(); jb++) {
	
	  cout << ib << " " << jb << " " << _h_mc_rec_decModeEff[i]->GetBinContent(ib,jb) << " " << totcol << endl;

	  _h_mc_rec_decModeEff[i]->SetBinContent( ib,jb, 100.*_h_mc_rec_decModeEff[i]->GetBinContent(ib,jb)/totcol );
	}
      }
    }
  }

  _fout->Write(0);
  _fout->Close(0);

  std::cout << "analyseTauProcessor::end()  " << name() 
 	    << std::endl ;
}

