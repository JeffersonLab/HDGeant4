//
// GlueXBeamConversionProcess class header
//
// author: richard.t.jones at uconn.edu
// version: december 24, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread,
// but virtually all of its functions need to be serialized, so
// it maintains its own interlocks for this purpose. Resources
// are created once when the first object is instantiated, and
// destroyed once when the last object is destroyed.

#ifndef _GLUEXBEAMCONVERSIONPROCESS_H_
#define _GLUEXBEAMCONVERSIONPROCESS_H_

#include <G4VDiscreteProcess.hh>
#include <G4VParticleChange.hh>
#include <ImportanceSampler.hh>
#include <PairConversionGeneration.hh>
#include <G4Step.hh>

class GlueXBeamConversionProcess: public G4VDiscreteProcess
{
 public:
   GlueXBeamConversionProcess(const G4String &name, 
                              G4ProcessType aType=fGeneral);
   GlueXBeamConversionProcess(GlueXBeamConversionProcess &src);
   virtual ~GlueXBeamConversionProcess();

   virtual G4double PostStepGetPhysicalInteractionLength(const G4Track &track,
                                                  G4double previousStepSize,
                                                  G4ForceCondition *condition);
   virtual G4VParticleChange *PostStepDoIt(const G4Track &track, 
                                           const G4Step &step);

 protected:
   void GenerateBeamPairConversion(const G4Step &step);
   virtual G4double GetMeanFreePath(const G4Track &track, 
                                    G4double previousStepSize,
                                    G4ForceCondition *condition);
   void GenerateBetheHeitlerProcess(const G4Step &step);

#if USING_DIRACXX
   static PairConversionGeneration *fPairsGeneration;
#endif

   int fStopBeamBeforeConverter;
   int fStopBeamAfterConverter;
   int fStopBeamAfterTarget;

   void prepareImportanceSamplingPDFs();

   static ImportanceSampler fPaircohPDF;
   static ImportanceSampler fTripletPDF;
   static G4double fBHpair_mass_min;

 private:
   GlueXBeamConversionProcess operator=(GlueXBeamConversionProcess &src);
};

#endif
