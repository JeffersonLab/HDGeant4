//
// GlueXSensitiveDetectorCTOF - class header
//
// author: staylor at jlab.org
// version: october 25, 2021
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorCTOF_h
#define GlueXSensitiveDetectorCTOF_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitCTOFbar.hh"
#include "GlueXHitCTOFpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorCTOF : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorCTOF(const G4String& name);
   GlueXSensitiveDetectorCTOF(const GlueXSensitiveDetectorCTOF &right);
   GlueXSensitiveDetectorCTOF &operator=(const GlueXSensitiveDetectorCTOF &right);
   virtual ~GlueXSensitiveDetectorCTOF();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* unused);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   GlueXHitsMapCTOFbar* fBarHitsMap;
   GlueXHitsMapCTOFpoint* fPointsMap;

   std::map<G4LogicalVolume*, int> fVolumeTable;

   static int MAX_HITS;
   static int MAX_HITS_PER_BAR;
   static double ATTENUATION_LENGTH;
   static double C_EFFECTIVE;
   static double FULL_BAR_LENGTH;
   static double TWO_HIT_TIME_RESOL;
   static double THRESH_MEV;
   static double MAX_TOF;

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif
