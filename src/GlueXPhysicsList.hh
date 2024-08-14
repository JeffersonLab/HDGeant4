//
// GlueXPhysicsList class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. The
// Geant4 runtime makes sure that its constructor runs in the master
// thread, and then its Construct methods are run again each time
// a new worker thread is initialized.

#ifndef GlueXPhysicsList_h
#define GlueXPhysicsList_h 1

#include <GlueXDetectorConstruction.hh>
#include <GlueXUserOptions.hh>
#include <GlueXSpecialCuts.hh>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4Version.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VEmProcess.hh"
#include "G4OpticalPhysics.hh"
#include "CompileTimeConstraints.hh"

class GlueXPhysicsList: public G4VModularPhysicsList
{
 public:
   GlueXPhysicsList(const GlueXDetectorConstruction *geometry=0, 
                    G4int verbosity=0);
   virtual ~GlueXPhysicsList();

   virtual void ConstructParticle();
   virtual void ConstructProcess();
   virtual void SetCuts();

   virtual void ListActiveProcesses();
   virtual void SelectActiveProcesses(G4int verbosity=0);

   virtual void DoMultipleScattering(G4int flag);
   virtual void DoBremsstrahlung(G4int flag);
   virtual void DoComptonScattering(G4int flag);
   virtual void DoIonizationEnergyLoss(G4int flag);
   virtual void DoPairConversion(G4int flag);
   virtual void DoParticleDecay(G4int flag);
   virtual void DoDeltaRayProduction(G4int flag);
   virtual void DoHadronicInteractions(G4int flag);
   virtual void DoCerenkovRadiation(G4int flag);
   virtual void DoOpticalAbsorption(G4int flag);

   virtual void DoProcessReordering();
   virtual void CheckProcessOrdering();

   G4VEmProcess *getBeamConversionProcess() const;
   G4VEmProcess *getBernardConversionProcess() const;
   G4VEmProcess *getKlongConversionProcess() const;
   G4OpticalPhysics *getOpticalPhysicsProcess() const;

 protected:
   GlueXUserOptions *fOptions;

   G4VEmProcess *fBeamConversion;
   G4VEmProcess *fBernardConversion;
   G4VEmProcess *fKlongConversion;
   G4OpticalPhysics *fOpticalPhysics;

#if G4VERSION_NUMBER < 1030
   // This member function gets introduced into base class
   // G4VUserPhysicsList in release Geant4.10.03, but until
   // we abandon ability to build under previous releases,
   // we need to have this member function defined.
   G4ParticleTable::G4PTblDicIterator* GetParticleIterator() const
   {
      return theParticleIterator;
   }
#endif
};

inline G4VEmProcess *GlueXPhysicsList::getBeamConversionProcess() const
{
   return fBeamConversion;
}

inline G4VEmProcess *GlueXPhysicsList::getBernardConversionProcess() const
{
   return fBernardConversion;
}

inline G4VEmProcess *GlueXPhysicsList::getKlongConversionProcess() const
{
   return fKlongConversion;
}

inline G4OpticalPhysics *GlueXPhysicsList::getOpticalPhysicsProcess() const
{
   return fOpticalPhysics;
}

#endif
