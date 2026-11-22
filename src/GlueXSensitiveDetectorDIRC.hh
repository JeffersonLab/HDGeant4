//
// GlueXSensitiveDetectorDIRC - class header
//
// author: richard.t.jones at uconn.edu
// version: november 28, 2016
//          may 7, 2017 R. Dzhygadlo - added adjustments for
//                      GlueXHitDIRCBar GlueXHitDIRCPmt    
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorDIRC_h
#define GlueXSensitiveDetectorDIRC_h 1

#include <TGraph.h>

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitDIRCflash.hh"
#include "GlueXHitDIRCPmt.hh"
#include "GlueXHitDIRCBar.hh"
#include "GlueXHitDIRCWob.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorDIRC : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorDIRC(const G4String& name);
   GlueXSensitiveDetectorDIRC(const GlueXSensitiveDetectorDIRC &right);
   GlueXSensitiveDetectorDIRC &operator=(const GlueXSensitiveDetectorDIRC &right);
   virtual ~GlueXSensitiveDetectorDIRC();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* ROhist);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

   static double GetDetectionEfficiency(double energy_GeV);
  
 private:
   std::vector<GlueXHitDIRCBar> fHitsBar;
   std::vector<GlueXHitDIRCWob> fHitsWob;
   std::vector<GlueXHitDIRCPmt> fHitsPmt;
   int fLutId;
   bool fLED;
   
   std::map<G4LogicalVolume*, int> fVolumeTable;
   
   static int MAX_HITS;
   static int MAX_PIXELS;
   // put all other detector response parameters here
   static double TWO_HIT_TIME_RESOL;
   
   static int instanceCount;
   static G4Mutex fMutex;
 
   static TGraph *fDetEff;
   static void InitializeDetEff();
};

#endif
