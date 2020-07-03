//
// GlueXEventAction class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
 
#include "GlueXEventAction.hh"
#include "GlueXUserOptions.hh"

#include "G4ios.hh"

GlueXEventAction::GlueXEventAction()
 : G4UserEventAction(),
   fProgressStep(0)
{}

void GlueXEventAction::BeginOfEventAction(const G4Event*)
{}

void GlueXEventAction::EndOfEventAction(const G4Event* evt)
{
   G4int event_id = evt->GetEventID();

   // periodic printing

   if (fProgressStep == 0) {
      long tens=1;
      while (tens * 10 <= event_id)
         tens *= 10;
      if (event_id % tens == 0)
         G4cout << event_id << " events simulated" << G4endl;
   }
   else if (event_id % fProgressStep == 0) {
      G4cout << event_id << " events simulated" << G4endl;
   }
}
