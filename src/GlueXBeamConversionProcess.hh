//
// GlueXBeamConversionProcess class header
//
// author: richard.t.jones at uconn.edu
// version: december 24, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread.
// Resources are created once when the first object is instantiated,
// and destroyed once when the last object is destroyed.

#ifndef _GLUEXBEAMCONVERSIONPROCESS_H_
#define _GLUEXBEAMCONVERSIONPROCESS_H_

#include "globals.hh"
#include <G4VEmProcess.hh>
#include <G4VParticleChange.hh>
#include <ImportanceSampler.hh>
#include <PairConversionGeneration.hh>
#include <AdaptiveSampler.hh>
#include <G4AutoLock.hh>
#include <G4Step.hh>

class GlueXBeamConversionProcess: public G4VEmProcess
{
 public:
   explicit GlueXBeamConversionProcess(const G4String &name = "beam_conversion", 
                                       G4ProcessType aType=fElectromagnetic);
   virtual ~GlueXBeamConversionProcess();

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

   enum FormFactorType {kGEp, kGMp, kF1p, kF2p,
                        kGEn, kGMn, kF1n, kF2n};
   virtual double nucleonFormFactor(double Q2_GeV, FormFactorType t);
   virtual double nuclearFormFactor(double Q2_GeV);
   virtual double nuclearFermiMomentum_GeV();
   virtual double nuclearBindingEnergy_GeV();
   virtual double nuclearMass_GeV();

 protected:
   virtual void InitialiseProcess(const G4ParticleDefinition*) override;
   virtual G4double GetMeanFreePath(const G4Track &track, 
                                    G4double previousStepSize,
                                    G4ForceCondition *condition) override;

   virtual void GenerateBeamPairConversion(const G4Step &step);
   virtual void GenerateBetheHeitlerProcess(const G4Step &step);
   virtual void GenerateTripletProcess(const G4Step &step);

   G4double fPIL;

#ifdef USING_DIRACXX
   PairConversionGeneration *fPairsGeneration;
#endif

   static int fStopBeamBeforeConverter;
   static int fStopBeamAfterConverter;
   static int fStopBeamAfterTarget;
   static int fLeptonPairFamily;
   static G4double fBHpair_mass_min;

   void prepareImportanceSamplingPDFs();

   ImportanceSampler *fPaircohPDF;
   ImportanceSampler *fTripletPDF;

   AdaptiveSampler *fAdaptiveSampler;
   void prepareAdaptiveSampler();

   void setConverterMaterial(double Z, double A);

 private:
   GlueXBeamConversionProcess() = delete;
   GlueXBeamConversionProcess(const GlueXBeamConversionProcess &src) = delete;
   GlueXBeamConversionProcess operator=(const GlueXBeamConversionProcess &src) = delete;

   G4bool  isInitialised;

   static G4Mutex fMutex;
   static int fConfigured;
   static std::vector<AdaptiveSampler*> fAdaptiveSamplerRegistry;

   static int fFormFactorChoice;

   double fTargetZ;
   double fTargetA;
};

inline void GlueXBeamConversionProcess::setConverterMaterial(double Z, double A)
{
   fTargetZ = Z;
   fTargetA = A;
}

#endif
