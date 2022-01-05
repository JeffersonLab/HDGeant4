//
// GlueXRunAction class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
 
#include "GlueXRunAction.hh"
#include "GlueXPhysicsList.hh"
#include "GlueXUserEventInformation.hh"

#include "G4VisManager.hh"
#include "G4ViewParameters.hh"
#include "G4VViewer.hh"
#include "G4ios.hh"

GlueXRunAction::GlueXRunAction(GlueXPhysicsList *plist)
 : fPhysicsList(plist)
{}

void GlueXRunAction::BeginOfRunAction(const G4Run*)
{
   G4VisManager *vmgr = dynamic_cast<G4VisManager*>
                        (G4VVisManager::GetConcreteInstance());
   if (vmgr) {
      G4VViewer *viewer = vmgr->GetCurrentViewer();
      if (viewer) {
         G4ViewParameters pars(viewer->GetViewParameters());
         pars.SetPicking(false);
         viewer->SetViewParameters(pars);
      }
   }

   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXRunAction constructor - "
             << "GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   std::map<int,int> save;
   if (user_opts->Find("SAVEHITS", save)) {
      if (save[1] != 0)
         GlueXUserEventInformation::setWriteNoHitEvents(1);
   }

   fPhysicsList->SelectActiveProcesses();
   //fPhysicsList->ListActiveProcesses();

   // Reorder processes to user-defined order
   fPhysicsList->DoProcessReordering();
}

void GlueXRunAction::EndOfRunAction(const G4Run* evt)
{}
