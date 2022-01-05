//
// GlueXSensitiveDetectorFMWPC - class header
//
// author: richard.t.jones at uconn.edu
// version: november 28, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorFMWPC_h
#define GlueXSensitiveDetectorFMWPC_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitFMWPCwire.hh"
#include "GlueXHitFMWPCpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorFMWPC : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorFMWPC(const G4String& name);
   GlueXSensitiveDetectorFMWPC(const GlueXSensitiveDetectorFMWPC &right);
   GlueXSensitiveDetectorFMWPC &operator=(const GlueXSensitiveDetectorFMWPC &right);
   virtual ~GlueXSensitiveDetectorFMWPC();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* unused);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   GlueXHitsMapFMWPCwire* fWireHitsMap;
   GlueXHitsMapFMWPCpoint* fPointsMap;

   std::map<G4LogicalVolume*, int> fVolumeTable;

   static int MAX_HITS;
   static double TWO_HIT_TIME_RESOL;
   static double THRESH_KEV;
   static double WIRE_OFFSET;
   static double WIRE_PITCH;

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif
