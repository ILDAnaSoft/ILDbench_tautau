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
    for (int i=0; i<3; i++) {
      _eigenValues[i]=-999;
      _eigenVectors[i].SetXYZ(999,999,999);
    }
    //_lcfi_evt = new vertex_lcfi::Event();
    //vertex_lcfi::MemoryManager<vertex_lcfi::Event>::Event()->registerObject( _lcfi_evt );
  }
  ~vertexInfo() {cleanup();}

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

  std::vector < vertex_lcfi::Track* > lcfi_tracks;
  std::vector < vertex_lcfi::TrackState* > lcfi_trackstates;
  //  vertex_lcfi::Event* _lcfi_evt;

  TVector3 _vtxPos;
  bool _vtxValid;
  float _chisq;

  TVector3 _eigenVectors[3];
  float _eigenValues[3];

  static const float _ipSize[3];

  std::vector < IMPL::ReconstructedParticleImpl* > myrecoparts;

  void cleanup();

};

#endif
