//
// GlueXSensitiveDetectorFCAL - class header
//
// author: richard.t.jones at uconn.edu
// version: october 25, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorFCAL_h
#define GlueXSensitiveDetectorFCAL_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitFCALblock.hh"
#include "GlueXHitFCALpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorFCAL : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorFCAL(const G4String& name);
   GlueXSensitiveDetectorFCAL(const GlueXSensitiveDetectorFCAL &right);
   GlueXSensitiveDetectorFCAL &operator=(const GlueXSensitiveDetectorFCAL &right);
   virtual ~GlueXSensitiveDetectorFCAL();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* ROhist);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   GlueXHitsMapFCALblock* fBlocksMap;
   GlueXHitsMapFCALpoint* fPointsMap;

   std::map<G4LogicalVolume*, int> fVolumeTable;

   static int MAX_HITS;
   static int CENTRAL_COLUMN;
   static int CENTRAL_ROW;
   static double WIDTH_OF_BLOCK;
   static double LENGTH_OF_BLOCK;
   static double ACTIVE_RADIUS;
   static double ATTENUATION_LENGTH;
   static double C_EFFECTIVE;
   static double TWO_HIT_TIME_RESOL;
   static double THRESH_MEV;
   static double SHOWER_ENERGY_SCALE_FACTOR;
   static double MIP_ENERGY_SCALE_FACTOR;

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif
