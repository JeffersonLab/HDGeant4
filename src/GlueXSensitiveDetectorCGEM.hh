//
// GlueXSensitiveDetectorCGEM - class header
//
// author: richard.t.jones at uconn.edu
// version: october 11, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorCGEM_h
#define GlueXSensitiveDetectorCGEM_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitCGEMlayer.hh"
#include "GlueXHitCGEMpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorCGEM : public G4VSensitiveDetector {
public:
  GlueXSensitiveDetectorCGEM(const G4String& name);
  GlueXSensitiveDetectorCGEM(const GlueXSensitiveDetectorCGEM &right);
  GlueXSensitiveDetectorCGEM &operator=(const GlueXSensitiveDetectorCGEM &right);
  virtual ~GlueXSensitiveDetectorCGEM();
  
  virtual void Initialize(G4HCofThisEvent* hitCollection);
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* unused);
  virtual void EndOfEvent(G4HCofThisEvent* hitCollection);
  
  int GetIdent(std::string div, const G4VTouchable *touch);
  
private:
  double Ei(double x);
  
 private:
  GlueXHitsMapCGEMlayer* fLayersMap;
  GlueXHitsMapCGEMpoint* fPointsMap;
  
  std::map<G4LogicalVolume*, int> fVolumeTable;
  
  static int MAX_HITS;
  static double TWO_HIT_TIME_RESOL;  
  static double THRESH_MEV;

  static int instanceCount;
  static G4Mutex fMutex;
};

#endif
