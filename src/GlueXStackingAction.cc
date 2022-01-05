//
// GlueXStackingAction class implementation
//
// author: richard.t.jones at uconn.edu
// version: january 29, 2017

#include "GlueXStackingAction.hh"
#include "GlueXUserOptions.hh"
#include "GlueXUserTrackInformation.hh"
#include "GlueXSensitiveDetectorDIRC.hh"

#include "G4Track.hh"
#include "G4ios.hh"
#include "G4ParticleTable.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalVolumeStore.hh"

#ifdef DIRC_MONITORING_HISTOS
#include "TH1.h"
#include "TCanvas.h"
TH1F *hbouncez = new TH1F("hbouncez",";bounces along z [#];entries [#]",1000,0,2000);
TH1F *hbouncey = new TH1F("hbouncey",";bounces along y [#];entries [#]",1000,0,2000);
#endif

GlueXStackingAction::GlueXStackingAction()
{
   if (!(G4ParticleTable::GetParticleTable()->GetReadiness())) {
      G4String msg;
      msg =  " You are instantiating GlueXStackingAction BEFORE your\n";
      msg += "G4VUserPhysicsList is instantiated and assigned to ";
      msg += "G4RunManager.\nSuch an instantiation is prohibited by ";
      msg += "Geant4 version 8.0. To fix this problem,\n";
      msg += "please make sure that your main() instantiates ";
      msg += "G4VUserPhysicsList AND\nset it to G4RunManager before ";
      msg += "instantiating other user action classes\n";
      msg += "such as GlueXStackingAction.";
      G4Exception("GlueXStackingAction::GlueXStackingAction()",
                  "Event0031",FatalException,msg);
   }

   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXStackingAction constructor - "
             << "GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   nosecondaries = 0;
   std::map<int,int> opt;
   if (user_opts->Find("NOSECONDARIES", opt)) {
      nosecondaries = (opt[1] != 0);
   }

   fBarEnd = 0;
   G4PhysicalVolumeStore* pvStore = G4PhysicalVolumeStore::GetInstance();
   for (size_t i=0; i<pvStore->size(); i++) {
      if ((*pvStore)[i]->GetName()=="QZWN") {
         fBarEnd = fabs((*pvStore)[i]->GetTranslation().x());
      }
   }
}

GlueXStackingAction::~GlueXStackingAction()
{
#ifdef DIRC_MONITORING_HISTOS
   TCanvas *c = new TCanvas("c","c",800,400);
   hbouncez->Draw();
   c->Print("cbounces_z.png");
   c->Print("cbounces_z.C");
   hbouncey->Draw();
   c->Print("cbounces_y.png");
   c->Print("cbounces_y.C");
#endif
}

G4ClassificationOfNewTrack GlueXStackingAction::ClassifyNewTrack(
                                                const G4Track *aTrack)
{
   // Determine what to do with this secondary, which has just been
   // made into a G4Track object and needs to be classified as follows.
   //
   //    enum G4ClassificationOfNewTrack
   //    {
   //      fUrgent,    // put into the urgent stack
   //      fWaiting,   // put into the waiting stack
   //      fPostpone,  // postpone to the next event
   //      fKill       // kill without stacking
   //    };
   //
   //    The parent_ID of the track indicates the origin of it.
   //                
   //    G4int parent_ID = aTrack->get_parentID();
   //   
   //      parent_ID = 0 : primary particle
   //                > 0 : secondary particle
   //                < 0 : postponed from the previous event

   if (nosecondaries && aTrack->GetParentID() != 0)
      return fKill;

   // apply detection efficiency for the DIRC at production stage:
   G4String ParticleName = aTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
   if (aTrack->GetParentID() != 0) { // particle is secondary
      if (ParticleName == "opticalphoton") {
         Double_t Ephoton = aTrack->GetMomentum().mag();
         Double_t ra = G4UniformRand();
         if (ra > GlueXSensitiveDetectorDIRC::GetDetectionEfficiency(Ephoton))
            return fKill;

         G4ThreeVector v = aTrack->GetPosition();
         G4ThreeVector n = aTrack->GetMomentumDirection().unit();

         Double_t bary = 35;     // bar width
         Double_t barz = 17.25;  // bar height
         Double_t barx = 4*1225; // bar length
 
         Double_t lenx;
         if (v.y()<0) {
           lenx = fabs(fBarEnd+v.x());
           if (n.x()>0)
              lenx = 2*barx - lenx;
          }
          else {
             lenx = fabs(v.x()-fBarEnd);
             if (n.x()<0)
                lenx = 2*barx - lenx;
          }

          Double_t lenz = lenx*n.z()/fabs(n.x());
          Double_t leny = lenx*n.y()/fabs(n.x());

          int bouncesz = fabs(lenz/barz);
          int bouncesy = fabs(leny/bary);

#ifdef DIRC_MONITORING_HISTOS
          hbouncez->Fill(bouncesz);
          hbouncey->Fill(bouncesy);
#endif

          Double_t anglez = fabs(n.getTheta()-CLHEP::pi/2.);
          Double_t angley = fabs(n.angle(G4ThreeVector(0,1,0))-CLHEP::pi/2.);

          Double_t lambda = 197.0*2.0*CLHEP::pi/(aTrack->GetMomentum().mag()*1.0E6);
          Double_t lambda2 = lambda*lambda; 

          // calculate bounce probability
          Double_t n_quartz = sqrt(1 + (0.696*lambda2/(lambda2-pow(0.068,2))) + (0.407*lambda2/(lambda2-pow(0.116,2))) + 0.897*lambda2/(lambda2-pow(9.896,2)));
          Double_t bounce_probz = 1 - pow(4*CLHEP::pi*cos(anglez)*0.5*n_quartz/lambda,2);// 0.5 [nm] - roughness
          Double_t bounce_proby = 1 - pow(4*CLHEP::pi*cos(angley)*0.5*n_quartz/lambda,2);    
          Double_t prob = pow(bounce_probz,bouncesz) * pow(bounce_proby,bouncesy);

          // transport efficiency
          if (G4UniformRand() > prob)
             return fKill;       
      } //else {return fKill;} // remove this condition!!!
   } // if particle is secondary
   return fUrgent;
}

void GlueXStackingAction::NewStage()
{
   // This method is called by G4StackManager when the urgentStack
   // becomes empty and contents in the waitingStack are transtered
   // to the urgentStack. Note that this method is not called at the
   // begining of each event, see "PrepareNewEvent" for that.
   //
   // In case re-classification of the stacked tracks is needed,
   // use the following method to request to G4StackManager.
   //
   //    stackManager->ReClassify();
   //
   // All of the stacked tracks in the waitingStack will be re-classified 
   // by "ClassifyNewTrack" method. To abort current event, do
   //
   //    stackManager->clear();
   //
   // Note that this way is valid and safe only for the case it is called
   // from this user class. The more global way of event abortion is
   //
   //    G4UImanager * UImanager = G4UImanager::GetUIpointer();
   //    UImanager->ApplyCommand("/event/abort");
}

void GlueXStackingAction::PrepareNewEvent()
{
   // This method is called by G4StackManager at the begining of each event.
   // Be careful that the urgentStack and the waitingStack of G4StackManager
   // are empty at this moment, because this method is called before 
   // accepting primary particles. Also, note that the postponeStack of
   // G4StackManager may have some postponed tracks.
}
