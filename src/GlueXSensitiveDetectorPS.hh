//
// GlueXSensitiveDetectorPS - class header
//
// author: richard.t.jones at uconn.edu
// version: october 28, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorPS_h
#define GlueXSensitiveDetectorPS_h 1

#include "G4VSensitiveDetector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include "GlueXHitPStile.hh"
#include "GlueXHitPSpoint.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorPS : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorPS(const G4String& name);
   GlueXSensitiveDetectorPS(const GlueXSensitiveDetectorPS &right);
   GlueXSensitiveDetectorPS &operator=(const GlueXSensitiveDetectorPS &right);
   virtual ~GlueXSensitiveDetectorPS();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* unused);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   GlueXHitsMapPStile* fTileHitsMap;
   GlueXHitsMapPSpoint* fPointsMap;

   static std::map<G4LogicalVolume*, int> fVolumeTable;

   static int MAX_HITS;
   static int NUM_COLUMNS_PER_ARM;
   static double TWO_HIT_TIME_RESOL;
   static double THRESH_MEV;

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif
