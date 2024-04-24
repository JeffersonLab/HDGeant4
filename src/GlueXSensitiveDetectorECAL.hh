//
// GlueXSensitiveDetectorECAL - class header
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorECAL_h
#define GlueXSensitiveDetectorECAL_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitECALblock.hh"
#include "GlueXHitECALpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorECAL : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorECAL(const G4String& name);
   GlueXSensitiveDetectorECAL(const GlueXSensitiveDetectorECAL &right);
   GlueXSensitiveDetectorECAL &operator=(const GlueXSensitiveDetectorECAL &right);
   virtual ~GlueXSensitiveDetectorECAL();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* ROhist);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   GlueXHitsMapECALblock* fBlocksMap;
   GlueXHitsMapECALpoint* fPointsMap;

   std::map<G4LogicalVolume*, int> fVolumeTable;

   static int MAX_HITS;
   static int CENTRAL_COLUMN;
   static int CENTRAL_ROW;
   static double WIDTH_OF_BLOCK;
   static double LENGTH_OF_BLOCK;
   static double ATTENUATION_LENGTH;
   static double C_EFFECTIVE;
   static double TWO_HIT_TIME_RESOL;
   static double THRESH_MEV;

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif
