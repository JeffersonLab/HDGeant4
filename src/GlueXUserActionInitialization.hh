//
// GlueXUserActionInitalization - class header
//
// author: richard.t.jones at uconn.edu
// version: august 24, 2015
//
// In the context of the Geant4 event-level multithreading model,
// this class is "shared", ie. has no thread-local state.

#ifndef _GlueXUserActionInitialization_
#define _GlueXUserActionInitialization_

#include <G4VUserActionInitialization.hh>

#include <GlueXPrimaryGeneratorAction.hh>
#include <GlueXRunAction.hh>
#include <GlueXEventAction.hh>
#include <GlueXStackingAction.hh>
#include <GlueXTrackingAction.hh>
#include <GlueXSteppingAction.hh>
#include <GlueXSteppingVerbose.hh>
#include <GlueXPhysicsList.hh>

class GlueXUserActionInitialization : public G4VUserActionInitialization
{
 public:
   GlueXUserActionInitialization(GlueXPhysicsList *plist)
    : fPhysicsList(plist) {}
   ~GlueXUserActionInitialization() {}
   
   virtual void Build() const {
      SetUserAction(new GlueXRunAction(fPhysicsList));
      SetUserAction(new GlueXEventAction());
      SetUserAction(new GlueXStackingAction());
      SetUserAction(new GlueXTrackingAction());
      SetUserAction(new GlueXSteppingAction());
      SetUserAction(new GlueXPrimaryGeneratorAction());
   }

   virtual void BuildForMaster() const {
      SetUserAction(new GlueXRunAction(fPhysicsList));
   }

   virtual G4VSteppingVerbose* InitializeSteppingVerbose() const {
      return new GlueXSteppingVerbose();
   }

 protected:
   GlueXPhysicsList *fPhysicsList;
};

#endif // _GlueXUserActionInitialization_
