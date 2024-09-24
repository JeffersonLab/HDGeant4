//
// GlueXKlongConversionProcess class header
//
// author: richard.t.jones at uconn.edu
// version: may 24, 2024
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread.
// Resources are created once when the first object is instantiated,
// and destroyed once when the last object is destroyed.

#ifndef _GLUEXKLONGCONVERSIONPROCESS_H_
#define _GLUEXKLONGCONVERSIONPROCESS_H_

#include <G4VEmProcess.hh>
#include <G4VEmModel.hh>
#include <G4ParticleChangeForGamma.hh>
#include <G4ProductionCutsTable.hh>

class GlueXKlongConversionProcess: public G4VEmProcess
{
 public:
   explicit GlueXKlongConversionProcess(const G4String &name = "Klong_conversion", 
                                        G4ProcessType aType=fElectromagnetic);
   virtual ~GlueXKlongConversionProcess();

   virtual G4double PostStepGetPhysicalInteractionLength(const G4Track &track,
                                       G4double previousStepSize,
                                       G4ForceCondition *condition) override;

   virtual G4bool IsApplicable(const G4ParticleDefinition&) final;

   // Print few lines of informations about the process: validity range,
   virtual void PrintInfo() override;

 protected:
   virtual void InitialiseProcess(const G4ParticleDefinition*) override;

 private:
   GlueXKlongConversionProcess operator=(GlueXKlongConversionProcess &src) = delete;
   GlueXKlongConversionProcess(GlueXKlongConversionProcess &src) = delete;

   G4bool  isInitialised;
};


class GlueXKlongConversionModel : public G4VEmModel
{
 public:
   explicit GlueXKlongConversionModel();
   virtual ~GlueXKlongConversionModel();

   virtual void Initialise(const G4ParticleDefinition*,
                           const G4DataVector&);

   virtual G4double ComputeCrossSectionPerAtom(
                           const G4ParticleDefinition* photon,
                                 G4double kinEnergy, 
                                 G4double Z, 
                                 G4double A=0, 
                                 G4double cut=0,
                                 G4double emax=DBL_MAX);

   virtual void SampleSecondaries(std::vector<G4DynamicParticle*>* secondaries,
                                  const G4MaterialCutsCouple* material_couple,
                                  const G4DynamicParticle* photon,
                                  G4double tmin,
                                  G4double maxEnergy);

 private:
   GlueXKlongConversionModel &operator=(const GlueXKlongConversionModel &right) = delete;
   GlueXKlongConversionModel(const GlueXKlongConversionModel&) = delete;

   G4bool isInitialised;
   G4int verboseLevel;

   G4ParticleChangeForGamma* fParticleChange;
 
 public:
   static double gammaPhiXS_dE;
   static double gammaPhiXS_ub[220];
   static double gammaPhiXS_scale_factor;
   static double gammaPhiXS_tslope;
   static double gammaPhiXS_Emin;
   static double gammaPhiXS_Emax;
};

#endif
