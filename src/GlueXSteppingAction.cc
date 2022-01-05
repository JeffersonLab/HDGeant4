//
// GlueXSteppingAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXSteppingAction.hh"
#include "GlueXUserOptions.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Threading.hh"
#include "G4Event.hh"
#include "G4Decay.hh"

#include <stdio.h>
#include <exception>
#include <map>

//#define BACKGROUND_PROFILING 1
#if BACKGROUND_PROFILING
#include <TFile.h>
#include <TTree.h>
TFile *bgprofiles_file = 0;
std::map<int, TTree*> bgprofiles;
#endif

G4Mutex GlueXSteppingAction::fMutex = G4MUTEX_INITIALIZER;
int GlueXSteppingAction::fStopTracksInCollimator = 0;
int GlueXSteppingAction::fSaveTrajectories = 0;

GlueXSteppingAction::GlueXSteppingAction()
{
   G4AutoLock barrier(&fMutex);
   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXSteppingAction constructor - "
             << "GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   std::map<int,double> showersInCol;
   if (user_opts->Find("SHOWERSINCOL", showersInCol)) {
      fStopTracksInCollimator = (showersInCol[1] == 0);
   }
   else {
      fStopTracksInCollimator = 0;
   }

   std::map<int,int> trajectories;
   if (user_opts->Find("TRAJECTORIES", trajectories)) {
      fSaveTrajectories = trajectories[1];
   }
   else {
      fSaveTrajectories = 0;
  }

#if BACKGROUND_PROFILING
   if (bgprofiles_file == 0) {
      bgprofiles_file = new TFile("bgprofiles.root", "recreate");
   }
#endif

}

GlueXSteppingAction::GlueXSteppingAction(const GlueXSteppingAction &src)
{}

GlueXSteppingAction::~GlueXSteppingAction()
{
#if BACKGROUND_PROFILING
   G4AutoLock barrier(&fMutex);
   if (bgprofiles_file) {
      std::map<int, TTree*>::iterator iter;
      for (iter = bgprofiles.begin(); iter != bgprofiles.end(); ++iter) {
         iter->second->Write();
      }
      delete bgprofiles_file;
      bgprofiles_file = 0;
   }
#endif
}

void GlueXSteppingAction::UserSteppingAction(const G4Step* step)
{
   G4Track *track = (G4Track*)step->GetTrack();
   const G4Event *event = G4RunManager::GetRunManager()
                                      ->GetCurrentEvent();
   GlueXUserEventInformation *eventinfo = (GlueXUserEventInformation*)
                                          event->GetUserInformation();
   GlueXUserTrackInformation *trackinfo;
   trackinfo = (GlueXUserTrackInformation*)track->GetUserInformation();
   const G4Step *hyperstep = G4ParallelWorldProcess::GetHyperStep();
   G4VPhysicalVolume *pvol = hyperstep->GetPostStepPoint()->GetPhysicalVolume();

   // Kill tracks at entry to primary collimator or active collimator
   // if this was asked for in the control file.
   if (fStopTracksInCollimator) {
      if (pvol && (pvol->GetName() == "INSU" ||
                   pvol->GetName() == "PCTT" ||
                   pvol->GetName() == "PCPB") )
      {
         track->SetTrackStatus(fStopAndKill);
      }
   }

   // Kill reflections inside DIRC's PMT
   if (pvol && (pvol->GetName() == "BWPV")) {
      track->SetTrackStatus(fStopAndKill);
   }

   // Kill tracks when they enter the walls / ceiling / floor,
   // otherwise a lot of time is spent showering in the walls.
   if (pvol && (pvol->GetName() == "World")) {
      track->SetTrackStatus(fStopAndKill);
   }

   // Post new vertices to the MC record for primary particle decays
   if (trackinfo) {
      int primeID = trackinfo->GetGlueXTrackID();
      if (primeID > 0) {
         const G4VProcess* process;
         process = step->GetPostStepPoint()->GetProcessDefinedStep();
         if (process && process->GetProcessType() == fDecay) {
            G4TrackVector &secondary = *(G4TrackVector*)
                                        step->GetSecondaryInCurrentStep();
            G4TrackVector::iterator iter;
            for (iter = secondary.begin(); iter != secondary.end(); ++iter) {
               int newID = eventinfo->AssignNextGlueXTrackID();
               trackinfo = new GlueXUserTrackInformation();
               if (eventinfo) {
                  trackinfo->SetGlueXTrackID(newID);
               }
               (*iter)->SetUserInformation(trackinfo);
            }
            if (eventinfo) {
               int mech[2];
               char *cmech = (char*)mech;
               snprintf(cmech, 5, "%c%c%c%c", 'D', 'C', 'A', 'Y');
               eventinfo->AddSecondaryVertex(secondary, primeID, mech[0]);
            }
         }
      }
   }
 
   // Save mc trajectory information if requested
   if (fSaveTrajectories) {
      eventinfo->AddMCtrajectoryPoint(*step, fSaveTrajectories);
   }

#if BACKGROUND_PROFILING

//  The following defines a ROOT tree containing information on particle 
//  type, energy, position, polarization, and at what virtual detector the 
//  particle passes through. These virtual detectors are filled with air
//  and are called "DETx" where x is an index currently between 1 and 8,
//  stored in the tree as integer element "det".  The xint[i][3] array
//  records the records the vertices of the interaction sequence leading
//  to the detected particle.

   struct profiler_row_t {
      float totE;
      float x[7];
      float ppol;
      float xspot[2];
      int ptype;
      int det;
      int mint;
      float xint[9][3];
   };
   const int mint_max = 9;
   static struct profiler_row_t prow[256];
   int id = G4Threading::G4GetThreadId() + 1;
   assert(id < 256);

   if (track->GetCurrentStepNumber() == 1) {
      if (track->GetParentID() == 0) {
         G4StepPoint *point = step->GetPostStepPoint();
         const G4ThreeVector &pos = point->GetPosition();
         prow[id].xspot[0] = pos[0]/cm;
         prow[id].xspot[1] = pos[1]/cm;
         prow[id].mint = 0;
      }
      else {
         G4StepPoint *point = step->GetPreStepPoint();
         const G4ThreeVector &pos = point->GetPosition();
         prow[id].xint[prow[id].mint][0] = pos[0]/cm;
         prow[id].xint[prow[id].mint][1] = pos[1]/cm;
         prow[id].xint[prow[id].mint][2] = pos[2]/cm;
         if (prow[id].mint < mint_max - 1)
            ++prow[id].mint;
      }
   }

   int det=0;
   if (pvol == 0)
      det = 0;
   else if (pvol->GetName() == "DET1")
      det = 1;
   else if (pvol->GetName() == "DET2")
      det = 2;
   else if (pvol->GetName() == "DET3")
      det = 3;
   else if (pvol->GetName() == "DET4")
      det = 4;
   else if (pvol->GetName() == "DET5")
      det = 5;
   else if (pvol->GetName() == "DET6")
      det = 6;
   else if (pvol->GetName() == "DET7")
      det = 7;
   else if (pvol->GetName() == "DET8")
      det = 8;
   else
      det = 0;

   TTree *proftree;
   if (det > 0) {
      G4AutoLock barrier(&fMutex);
      try {
         proftree = bgprofiles.at(det);
      }
      catch(std::out_of_range &err) {
         std::stringstream names;
         names << "det" << det;
         std::stringstream titles;
         titles << "hits in virtual detector " << det;
         proftree = new TTree(names.str().c_str(), titles.str().c_str());
         proftree->Branch("totE", &prow[0].totE, "totE/F");
         proftree->Branch("x", &prow[0].x[0], "x[7]/F");
         proftree->Branch("ppol", &prow[0].ppol, "ppol/F");
         proftree->Branch("xspot", &prow[0].xspot[0], "xspot[2]/F");
         proftree->Branch("ptype", &prow[0].ptype, "ptype/I");
         proftree->Branch("det", &prow[0].det, "det/I");
         proftree->Branch("mint", &prow[0].mint, "mint/I");
         proftree->Branch("xint", &prow[0].xint[0][0], "xint[mint][3]/F");
         bgprofiles[det] = proftree;
      }
      G4StepPoint *point = step->GetPostStepPoint();
      const G4ThreeVector &pos = point->GetPosition();
      const G4ThreeVector &mom = point->GetMomentum();
      const G4ThreeVector &pol = point->GetPolarization();
      double pmag = mom.mag();
      double Etot = point->GetTotalEnergy();
      int pdgcode = track->GetDynamicParticle()->GetPDGcode();
      int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgcode);
      prow[0].totE = Etot/GeV;
      prow[0].x[0] = pos[0]/cm;
      prow[0].x[1] = pos[1]/cm;
      prow[0].x[2] = pos[2]/cm;
      prow[0].x[3] = mom[0]/pmag;
      prow[0].x[4] = mom[1]/pmag;
      prow[0].x[5] = mom[2]/pmag;
      prow[0].x[6] = pmag/GeV;
      prow[0].ppol= pol.mag();
      prow[0].xspot[0] = prow[id].xspot[0];
      prow[0].xspot[1] = prow[id].xspot[1];
      prow[0].ptype = g3type;
      prow[0].det = det;
      prow[0].mint = prow[id].mint;
      for (int i=0; i < prow[0].mint; ++i) {
         prow[0].xint[i][0] = prow[id].xint[i][0];
         prow[0].xint[i][1] = prow[id].xint[i][1];
         prow[0].xint[i][2] = prow[id].xint[i][2];
      }
      proftree->Fill();
   }

#endif

}
