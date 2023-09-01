//
// GlueXSensitiveDetectorITOF - class header
//
// author: staylor at jlab.org
// version: february 7, 2022
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorITOF_h
#define GlueXSensitiveDetectorITOF_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitITOFhit.hh"
#include "GlueXHitITOFpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorITOF : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorITOF(const G4String& name);
   GlueXSensitiveDetectorITOF(const GlueXSensitiveDetectorITOF &right);
   GlueXSensitiveDetectorITOF &operator=(const GlueXSensitiveDetectorITOF &right);
   virtual ~GlueXSensitiveDetectorITOF();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* unused);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
  GlueXHitsMapITOFhit* fHitsMap;
  GlueXHitsMapITOFpoint* fPointsMap;
  
  std::map<G4LogicalVolume*, int> fVolumeTable;
  
  static int MAX_HITS;

  static int instanceCount;
  static G4Mutex fMutex;
};

#endif
