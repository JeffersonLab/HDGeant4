//
// GlueXSensitiveDetectorBCAL - class header
//
// author: richard.t.jones at uconn.edu
// version: october 25, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorBCAL_h
#define GlueXSensitiveDetectorBCAL_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitBCALcell.hh"
#include "GlueXHitBCALpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorBCAL : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorBCAL(const G4String& name);
   GlueXSensitiveDetectorBCAL(const GlueXSensitiveDetectorBCAL &right);
   GlueXSensitiveDetectorBCAL &operator=(const GlueXSensitiveDetectorBCAL &right);
   virtual ~GlueXSensitiveDetectorBCAL();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* ROhist);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   GlueXHitsMapBCALcell* fCellsMap;
   GlueXHitsMapBCALpoint* fPointsMap;

   std::map<G4LogicalVolume*, int> fVolumeTable;

   static int MAX_HITS;
   static double THRESH_MEV;
   static double TWO_HIT_TIME_RESOL;
   static double ATTENUATION_LENGTH;
   static double C_EFFECTIVE;
   static double SIPM_TIME_BIN_WIDTH;
   static double MODULE_FULL_LENGTH;
   static double ATTENUATION_FULL_LENGTH;
   static double THRESH_ATTENUATED_GEV;
   static double SHOWER_ENERGY_SCALE_FACTOR;

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif
