//
// GlueXPhysicsList class header
//
// author: richard.t.jones at uconn.edu
// version: january 22, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. The
// Geant4 runtime makes sure that its constructor runs in the master
// thread, and then its Construct methods are run again each time
// a new worker thread is initialized.

#ifndef GlueXSpecialCuts_h
#define GlueXSpecialCuts_h 1

#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4UserLimits.hh"
#include "G4AutoLock.hh"

class G4LossTableManager;

class GlueXSpecialCuts : public G4VProcess 
{
 public:
   GlueXSpecialCuts(const G4String& processName ="GlueXSpecialCut" );
   virtual ~GlueXSpecialCuts();

   virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            );

   virtual G4VParticleChange* PostStepDoIt(
                             const G4Track& ,
                             const G4Step& 
                            );

   //  no operation in  AtRestGPIL
   virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& ,
                             G4ForceCondition* 
                            ){ return -1.0; };

   //  no operation in  AtRestDoIt      
   virtual G4VParticleChange* AtRestDoIt(
                             const G4Track& ,
                             const G4Step&
                            ){return 0;};

   //  no operation in  AlongStepGPIL
   virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track&,
                             G4double  ,
                             G4double  ,
                             G4double& ,
                             G4GPILSelection*
                            ){ return -1.0; };

   //  no operation in  AlongStepDoIt
   virtual G4VParticleChange* AlongStepDoIt(
                             const G4Track& ,
                             const G4Step& 
                            ) {return 0;};

   const G4UserLimits *GetUserLimits() const;
   void SetUserLimits(const G4UserLimits *userlimits);

 private:
   // hide assignment operator as private 
   GlueXSpecialCuts(GlueXSpecialCuts&);
   GlueXSpecialCuts& operator=(const GlueXSpecialCuts& right);

   G4LossTableManager* theLossTableManager;
   static G4UserLimits *fUserLimits;
   static G4Mutex fMutex;
};

#endif
