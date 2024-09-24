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

//#define VERBOSE_PRINTOUT_OF_GAMMA_PROCESSES 1
#include "G4Gamma.hh"

#define BACKGROUND_PROFILING 1
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
#ifdef VERBOSE_PRINTOUT_OF_GAMMA_PROCESSES
    // Get the current track
    G4Track* track = step->GetTrack();

    // Get the particle type
    G4ParticleDefinition* particle = track->GetDefinition();

    // Only print information for gamma particles
    if (particle == G4Gamma::GammaDefinition()) {
        // Get the process manager of the gamma particle
        G4ProcessManager* processManager = particle->GetProcessManager();
        if (processManager) {
            G4ProcessVector* processVector = processManager->GetProcessList();
            G4int numProcesses = processVector->size();

            // Print the current step number and particle information
            G4cout << "Track ID: " << track->GetTrackID()
                   << " | Step Number: " << track->GetCurrentStepNumber()
                   << " | Gamma particle" << G4endl;

            // Loop over all processes
            for (G4int i = 0; i < numProcesses; i++) {
                G4VProcess* process = (*processVector)[i];

                // Get the physical interaction length (GPIL) for each process
				G4ForceCondition forced;
                G4double gpil = process->PostStepGetPhysicalInteractionLength(*track, step->GetStepLength(), &forced);
                
                G4cout << "  Process: " << process->GetProcessName()
                       << " | GPIL: " << gpil / mm << " mm" 
                       << " | forced: " << ((forced == InActivated)? "InActivated" : (forced == Forced)? "Forced" : (forced == NotForced)? "NotForced" :
                                            (forced == Conditionally)? "Conditionally" : (forced == ExclusivelyForced)? "ExclusivelyForced" :
                                            (forced == StronglyForced)? "StronglyForced" : "(unknown code)") << G4endl;
            }
        }
    }
#else
   G4Track *track = (G4Track*)step->GetTrack();
#endif
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
            G4TrackVector decaySecondaries;
            for (iter = secondary.begin(); iter != secondary.end(); ++iter) {
               if ((*iter)->GetCreatorProcess()->GetProcessType() == fDecay) {
                  int newID = eventinfo->AssignNextGlueXTrackID();
                  trackinfo = new GlueXUserTrackInformation();
                  if (eventinfo) {
                     trackinfo->SetGlueXTrackID(newID);
                  }
                  (*iter)->SetUserInformation(trackinfo);
                  decaySecondaries.push_back(*iter);
               }
            }
            if (eventinfo and decaySecondaries.size() > 0) {
               int mech[2];
               char *cmech = (char*)mech;
               snprintf(cmech, 5, "%c%c%c%c", 'D', 'C', 'A', 'Y');
               eventinfo->AddSecondaryVertex(decaySecondaries, primeID, mech[0]);
            }
         }
      }
   }
 
   // Post new vertices to the MC record for Klong conversions
   if (trackinfo) {
      const G4VProcess* process;
      process = step->GetPostStepPoint()->GetProcessDefinedStep();
      if (process && process->GetProcessName() == "KlongConversion") {
         G4TrackVector &secondary = *(G4TrackVector*)
                                     step->GetSecondaryInCurrentStep();
         G4TrackVector::iterator iter;
         G4TrackVector decaySecondaries;
         for (iter = secondary.begin(); iter != secondary.end(); ++iter) {
            int newID = eventinfo->AssignNextGlueXTrackID();
            trackinfo = new GlueXUserTrackInformation();
            if (eventinfo) {
               trackinfo->SetGlueXTrackID(newID);
            }
            (*iter)->SetUserInformation(trackinfo);
            decaySecondaries.push_back(*iter);
         }
         if (eventinfo and decaySecondaries.size() > 0) {
            int mech[2];
            char *cmech = (char*)mech;
            snprintf(cmech, 5, "%c%c%c%c", 'C', 'O', 'N', 'V');
            eventinfo->AddSecondaryVertex(decaySecondaries, 1, mech[0]);
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
//  stored in the tree as integer element "det".  The Xint[i] arrays
//  contain the record of the vertices in the interaction sequence leading
//  to the detected particle, where i=0 refers to the birth vertex of the
//  ptype particle itself, i=1 to its immediate ancestor, and so on.

   const int mint_max = 999;
   struct profiler_row_t {
      float totE;                // total energy (GeV)
      float time;                // global tracking time at det
      float x[7];                // x,y,z,p,px/p,py/y,pz/p (cm, GeV)
      float ppol;                // polarization (if any)
      float xspot[2];            // x,y of primary vertex (cm)
      int ptype;                 // particle type (G3 code)
      int det;                   // index of detected volume
      int mint;                  // interaction history depth
      float xint[mint_max][3];   // x,y,z of ancestor birthplace (cm)
      float tint[mint_max];      // global time of ancestor birth (ns)
      float pint[mint_max];      // particle type of ancestor (G3 code)
   };
   static struct profiler_row_t prow[256];

   struct parent_history_t {
      int parent_id;
      float xint[3];
      float tint;
      float pint;
   };
   static std::map<int, struct parent_history_t> parent_history[256];

   int id = G4Threading::G4GetThreadId() + 1;
   assert(id < 256);

   if (track->GetCurrentStepNumber() == 1) {
      int trackId = track->GetTrackID();
      int parentId = track->GetParentID();
      G4StepPoint *point = step->GetPreStepPoint();
      const G4ThreeVector &pos = point->GetPosition();
      if (parentId == 0) {
         prow[id].xspot[0] = pos[0]/cm;
         prow[id].xspot[1] = pos[1]/cm;
         parent_history[id].clear();
      }
      parent_history[id][trackId].parent_id = parentId;
      parent_history[id][trackId].xint[0] = pos[0]/cm;
      parent_history[id][trackId].xint[1] = pos[1]/cm;
      parent_history[id][trackId].xint[2] = pos[2]/cm;
      parent_history[id][trackId].tint = point->GetGlobalTime()/ns;
      int pdgcode = track->GetDynamicParticle()->GetPDGcode();
      int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgcode);
      parent_history[id][trackId].pint = g3type;
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
         proftree->Branch("time", &prow[0].time, "time/F");
         proftree->Branch("x", &prow[0].x[0], "x[7]/F");
         proftree->Branch("ppol", &prow[0].ppol, "ppol/F");
         proftree->Branch("xspot", &prow[0].xspot[0], "xspot[2]/F");
         proftree->Branch("ptype", &prow[0].ptype, "ptype/I");
         proftree->Branch("det", &prow[0].det, "det/I");
         proftree->Branch("mint", &prow[0].mint, "mint/I");
         proftree->Branch("xint", &prow[0].xint[0][0], "xint[mint][3]/F");
         proftree->Branch("tint", &prow[0].tint[0], "tint[mint]/F");
         proftree->Branch("pint", &prow[0].pint[0], "pint[mint]/I");
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
      prow[0].time = point->GetGlobalTime()/ns;
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
      int trackId = track->GetTrackID();
      for (int i=0; i < mint_max; ++i) {
         prow[0].xint[i][0] = parent_history[id][trackId].xint[0];
         prow[0].xint[i][1] = parent_history[id][trackId].xint[1];
         prow[0].xint[i][2] = parent_history[id][trackId].xint[2];
         prow[0].tint[i] = parent_history[id][trackId].tint;
         prow[0].pint[i] = parent_history[id][trackId].pint;
         prow[0].mint = i + 1;
         trackId = parent_history[id][trackId].parent_id;
         if (trackId == 0)
            break;
      }
      proftree->Fill();
   }

#endif

}
