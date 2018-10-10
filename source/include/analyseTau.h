#ifndef analyseTauProcessor_h
#define analyseTauProcessor_h 1
#include "marlin/Processor.h"

using namespace marlin ;

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

class analyseTauProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new analyseTauProcessor ; }
  
  analyseTauProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
    
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 

  virtual void check( LCEvent * evt ) ; 
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;

 protected:

  std::string _outfile;
  TFile* _fout;
  TH2F* _h_eonp_ecalFrac[3];

  TH2F* _h_mc_rec_decMode[4];
  TH2F* _h_mc_rec_decModeEff[4];

  TH2F* _h_mcdec_jetMass[4];
  TH2F* _h_mcdec_nChg[4];
  TH2F* _h_mcdec_nTrk[4];
  TH2F* _h_mcdec_nGam[4];
  TH2F* _h_mcdec_nPi0[4];
  TH2F* _h_mcdec_nPi0Gam[4];

  //  enum {tdec_e=0, tdec_m, tdec_pi, tdec_rho, tdec_a1p, tdec_had1p, tdec_a3p, tdec_had3p, tdec_had5p};


};


#endif



