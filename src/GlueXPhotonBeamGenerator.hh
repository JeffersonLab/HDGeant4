//
// GlueXPhotonBeamGenerator class header
//
// author: richard.t.jones at uconn.edu
// version: december 24, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread.
// Resources are created once when the first object is instantiated,
// and destroyed once when the last object is destroyed.

#ifndef GlueXPhotonBeamGenerator_H
#define GlueXPhotonBeamGenerator_H

#include <G4VPrimaryGenerator.hh>
#include <G4GenericMessenger.hh>
#include <CobremsGeneration.hh>
#include <ImportanceSampler.hh>
#include <GlueXPseudoDetectorTAG.hh>
#include <G4Event.hh>

class GlueXPhotonBeamGenerator: public G4VPrimaryGenerator
{
 public:
   GlueXPhotonBeamGenerator(CobremsGeneration *gen);
   virtual ~GlueXPhotonBeamGenerator();

   virtual void GeneratePrimaryVertex(G4Event *event);
   virtual void GenerateBeamPhoton(G4Event *event, double t0);

   double GenerateTriggerTime(const G4Event *event);
   int GenerateTaggerHit(const G4Event *event, 
                         double energy, double time, int bg=0);
   void GenerateRFsync(const G4Event *event);

 protected:
   CobremsGeneration *fCobrems;
   GlueXPseudoDetectorTAG *fTagger;

   static int fGenerateNotSimulate;
   static int fBeamBackgroundTagOnly;

   static double fBeamBucketPeriod;
   static double fBeamStartZ;
   static double fBeamDiameter;
   static double fBeamVelocity;
   static double fBeamOffset[2];

   ImportanceSampler fCoherentPDFx; 
   ImportanceSampler fIncoherentPDFlogx;
   ImportanceSampler fIncoherentPDFy;
   double fIncoherentPDFtheta02;

   void prepareImportanceSamplingPDFs();

   static int fForceFixedPolarization;
   static double fFixedPolarization;
   static double fFixedPolarization_phi;

   G4GenericMessenger *fMessenger;

 public:
   static void setBeamDiameter(double D) {
      fBeamDiameter = D;
   }
   static double getBeamDiameter() {
      return fBeamDiameter;
   }
   static double getBeamVelocity() {
      return fBeamVelocity;
   }
   static void setBeamOffset(double x, double y) {
      fBeamOffset[0] = x;
      fBeamOffset[1] = y;
   }
   static double getBeamOffset(int i) {
      return fBeamOffset[i];
   }
   static double getRFreferencePlaneZ(int runno=0);
   static double getBeamBucketPeriod(int runno=0);
   static void setBeamBucketPeriod(double period) {
      fBeamBucketPeriod = period;
   }
   static void setBeamStartZ(double z) {
      fBeamStartZ = z;
   }
   static double getBeamStartZ() {
      return fBeamStartZ;
   }
   void enableFixedPolarization(double polar, double phi_deg) {
      fFixedPolarization = polar;
      fFixedPolarization_phi = phi_deg * M_PI/180;
      fForceFixedPolarization = true;
   }
   void disableFixedPolarization() {
      fForceFixedPolarization = false;
   }

 private:
   GlueXPhotonBeamGenerator();
   GlueXPhotonBeamGenerator(const GlueXPhotonBeamGenerator &src) = delete;
   GlueXPhotonBeamGenerator &operator=(const GlueXPhotonBeamGenerator &src) = delete;
};

#endif
