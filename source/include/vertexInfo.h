#ifndef _MYVTX_H_
#define _MYVTX_H_

#include "EVENT/Track.h"

#include "vertex_lcfi/inc/event.h"
#include "vertex_lcfi/inc/track.h"

#include "vertex_lcfi/util/inc/memorymanager.h"

#include "IMPL/ReconstructedParticleImpl.h"

#include "TVector3.h"

class vertexInfo{
 public:
  vertexInfo() {
    _vtxValid=false; 
    _vtxPos.SetXYZ(999,999,999); 
    _chisq=9999;
    _useIPcons =false;
    _trimTracks=false;
    // bfield
    _bField=3.5;
    // IP region size in mm, used for IP constraint
    _ipSize[0] = 10.e-3;
    _ipSize[1] = 10.e-6;
    _ipSize[2] = 0.3;
    for (int i=0; i<3; i++) {
      _eigenValues[i]=-999;
      _eigenVectors[i].SetXYZ(999,999,999);
    }
  }
  ~vertexInfo() {cleanup();}

  void setBeamSpotSize( float* xyz ) { 
    for ( int i=0; i<3; i++) _ipSize[i]=xyz[i];
  }

  void setBField(float field) {_bField=field;}
  float getBField() {return _bField;}

  void addTrack( EVENT::Track* trk );

  TVector3 getVertexPosition() {calculateVertexPosition(); return _vtxPos;}
  void calculateVertexPosition( );
  int getNtrack() { return lcfi_tracks.size() ; }
  float getVertexChisq() {calculateVertexPosition(); return _chisq;}

  TVector3 getEigenVector(int i) { assert (i>=0 && i<3); calculateVertexPosition(); return _eigenVectors[i];}
  float getEigenValue(int i) {assert (i>=0 && i<3); calculateVertexPosition(); return _eigenValues[i]; }

  float getVertexZ0(float x=0, float y=0);

  bool isValid() {return _vtxValid;}

  void trimTracks( bool b=true ) { _trimTracks=b; _vtxValid=false; }
  void useIPcon( bool b=true ) { _useIPcons=b; _vtxValid=false; }

 private:

  bool _useIPcons;
  bool _trimTracks;
  float _bField;

  std::vector < vertex_lcfi::Track* > lcfi_tracks;
  std::vector < vertex_lcfi::TrackState* > lcfi_trackstates;

  TVector3 _vtxPos;
  bool _vtxValid;
  float _chisq;

  TVector3 _eigenVectors[3];
  float _eigenValues[3];

  float _ipSize[3];

  std::vector < IMPL::ReconstructedParticleImpl* > myrecoparts;

  void cleanup();

};

#endif
