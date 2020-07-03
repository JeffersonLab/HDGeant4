//
// GlueXEventAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread.
 
#ifndef GlueXEventAction_h
#define GlueXEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4Event.hh"

class G4Event;

class GlueXEventAction : public G4UserEventAction
{
 public:
   GlueXEventAction();
   void BeginOfEventAction(const G4Event*);
   void EndOfEventAction(const G4Event*);

 protected:
   int fProgressStep; 
};

#endif
