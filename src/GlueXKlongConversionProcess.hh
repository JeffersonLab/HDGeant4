//
// GlueXKlongConversionProcess class header
//
// author: richard.t.jones at uconn.edu
// version: may 24, 2021
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread.
// Resources are created once when the first object is instantiated,
// and destroyed once when the last object is destroyed.

#ifndef _GLUEXKLONGCONVERSIONPROCESS_H_
#define _GLUEXKLONGCONVERSIONPROCESS_H_

#include "globals.hh"
#include <G4VEmProcess.hh>
#include <G4VParticleChange.hh>
#include <G4Gamma.hh>
#include <G4AutoLock.hh>

class GlueXKlongConversionProcess: public G4VEmProcess
{
 public:
   explicit GlueXKlongConversionProcess(const G4String &name = "Klong_conversion", 
                                          G4ProcessType aType=fElectromagnetic);
   virtual ~GlueXKlongConversionProcess();

   virtual G4double PostStepGetPhysicalInteractionLength(const G4Track &track,
                                       G4double previousStepSize,
                                       G4ForceCondition *condition) override;
   virtual G4VParticleChange *PostStepDoIt(const G4Track &track, 
                                           const G4Step &step) override;

   virtual G4bool IsApplicable(const G4ParticleDefinition&) final;

   virtual G4double MinPrimaryEnergy(const G4ParticleDefinition*,
				    const G4Material*) override;

   // Print few lines of informations about the process: validity range,
   virtual void PrintInfo() override;

   // print documentation in html format
   virtual void ProcessDescription(std::ostream&) const override;

 protected:
   virtual void InitialiseProcess(const G4ParticleDefinition*) override;

   virtual void GenerateConversionVertex(const G4Track &track,
                                         const G4Step &step);
   virtual G4double GetMeanFreePath(const G4Track &track, 
                                    G4double previousStepSize,
                                    G4ForceCondition *condition) override;

   G4double fPIL;

   static int fStopBeamAfterConverter;

 private:
   GlueXKlongConversionProcess operator=(GlueXKlongConversionProcess &src) = delete;
   GlueXKlongConversionProcess(GlueXKlongConversionProcess &src) = delete;

   G4bool  isInitialised;

   static G4Mutex fMutex;
   static int fConfigured;
};

#endif
