//
// GlueXSteppingAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread.

#ifndef GlueXSteppingAction_h
#define GlueXSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

class GlueXSteppingAction : public G4UserSteppingAction
{
 public:
   GlueXSteppingAction();
   GlueXSteppingAction(const GlueXSteppingAction &src);
   ~GlueXSteppingAction();
   void UserSteppingAction(const G4Step*);

 private:
   GlueXSteppingAction &operator=(const GlueXSteppingAction &src);

 protected:
   static int fStopTracksInCollimator;
   static int fSaveTrajectories;
   static G4Mutex fMutex;
};

#endif
