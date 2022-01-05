//
// GlueXUserActionInitalization - class header
//
// author: richard.t.jones at uconn.edu
// version: january 29, 2017
//
// In the context of the Geant4 event-level multithreading model,
// this class is "shared", ie. has no thread-local state.

#ifndef GlueXStackingAction_h
#define GlueXStackingAction_h 1

#include "G4UserStackingAction.hh"
#include "G4ClassificationOfNewTrack.hh"

class G4StackManager;
class G4Track;

class GlueXStackingAction : public G4UserStackingAction
{
 public:
   GlueXStackingAction();
   virtual ~GlueXStackingAction();

   virtual G4ClassificationOfNewTrack 
           ClassifyNewTrack(const G4Track* aTrack);
   virtual void NewStage();
   virtual void PrepareNewEvent();

 protected:
   int nosecondaries;

 private:
  double fBarEnd;
};

#endif

