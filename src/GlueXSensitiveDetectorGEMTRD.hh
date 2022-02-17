//
// GlueXSensitiveDetectorGEMTRD - class header
//
// author: staylor at jlab.org
// version: february 7, 2022
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorGEMTRD_h
#define GlueXSensitiveDetectorGEMTRD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitGEMTRDhit.hh"
#include "GlueXHitGEMTRDpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorGEMTRD : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorGEMTRD(const G4String& name);
   GlueXSensitiveDetectorGEMTRD(const GlueXSensitiveDetectorGEMTRD &right);
   GlueXSensitiveDetectorGEMTRD &operator=(const GlueXSensitiveDetectorGEMTRD &right);
   virtual ~GlueXSensitiveDetectorGEMTRD();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* unused);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
  GlueXHitsMapGEMTRDhit* fHitsMap;
  GlueXHitsMapGEMTRDpoint* fPointsMap;
  
  std::map<G4LogicalVolume*, int> fVolumeTable;
  
  static int MAX_HITS;
  static const double ELECTRON_CHARGE;
  static double GAS_GAIN;
  static double W_EFF_PER_ION;
  static double N_SECOND_PER_PRIMARY;

  static int instanceCount;
  static G4Mutex fMutex;
};

#endif
