//
// class implementation for GlueXBernardConversionProcess
//
// author: richard.t.jones at uconn.edu
// version: december 24, 2016
//

#include "GlueXBernardConversionProcess.hh"
#include "GlueXPhotonBeamGenerator.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"
#include "GlueXUserOptions.hh"
#include "G4ParallelWorldProcess.hh"

#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

// If you set this flag to 1 then all beam photons that reach
// the TPOL converter target will convert to e+e- pairs inside,
// otherwise the standard pair conversion probabilities apply.
#define FORCED_PTAR_PAIR_CONVERSION 0

// If you set this flag to 1 then all beam photons that reach
// the LIH2 target will undergo Bethe-Heitler conversion to e+e-
// pairs inside, otherwise the standard pair conversion
// probabilities apply.
#define FORCED_LIH2_PAIR_CONVERSION 0
#define FORCED_LIHE_PAIR_CONVERSION 0
#define FORCED_BETG_PAIR_CONVERSION 0
#define FORCED_TGT0_PAIR_CONVERSION 0

#include <G4SystemOfUnits.hh>
#include "G4PhysicalConstants.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4RunManager.hh"
#include "G4TrackVector.hh"
#include "G4GammaConversion.hh"
#include "G4PairProductionRelModel.hh"
#include "G4EmParameters.hh"

#include <stdio.h>
#include <iomanip>

G4Mutex GlueXBernardConversionProcess::fMutex = G4MUTEX_INITIALIZER;
int GlueXBernardConversionProcess::fStopBeamBeforeConverter = 0;
int GlueXBernardConversionProcess::fStopBeamAfterConverter = 0;
int GlueXBernardConversionProcess::fStopBeamAfterTarget = 0;
int GlueXBernardConversionProcess::fLeptonPairFamily = 0;
int GlueXBernardConversionProcess::fConfigured = 0;

GlueXBernardConversionProcess::GlueXBernardConversionProcess(
                                           const G4String &name, 
                                           G4ProcessType aType)
 : G4VEmProcess(name, aType),
   isInitialised(false)
{
   SetMinKinEnergy(2.0*electron_mass_c2);
   SetProcessSubType(fGammaConversion);
   SetStartFromNullFlag(true);
   SetBuildTableFlag(true);
   SetSecondaryParticle(G4Electron::Electron());
   SetLambdaBinning(220);

   verboseLevel = 0;

   G4AutoLock barrier(&fMutex);
   if (fConfigured)
      return;

   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXBernardConversionProcess constructor - "
             << "GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   fStopBeamBeforeConverter = 0;
   fStopBeamAfterConverter = 0;
   fStopBeamAfterTarget = 0;

   std::map<int,std::string> infile;
   std::map<int,double> beampars;
   std::map<int,std::string> genbeampars;
   if (!user_opts->Find("INFI", infile) &&
       user_opts->Find("GENBEAM", genbeampars))
   {
      if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "POSTCONV_BERNARD" ||
               genbeampars[1] == "postconv_bernard" ||
               genbeampars[1] == "Postconv_Bernard" ||
               genbeampars[1] == "PostConv_Bernard" ))
      {
         fStopBeamAfterConverter = 1;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHgen_Bernard" ||
               genbeampars[1] == "BHGEN_BERNARD" ||
               genbeampars[1] == "bhgen_bernard" ))
      {
         fStopBeamAfterTarget = 1;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHmuons_Bernard" ||
               genbeampars[1] == "BHmuons_BERNARD" ||
               genbeampars[1] == "bhmuons_bernard" ))
      {
         fStopBeamAfterTarget = 1;
         fLeptonPairFamily = 1;
      }
   }

   fConfigured = 1;

   if (verboseLevel > 0) {
       G4cout << GetProcessName() << " is created " << G4endl
              << "    Stop beam before converter? "
              << (fStopBeamBeforeConverter? "yes" : "no") << G4endl
              << "    Stop beam after converter? "
              << (fStopBeamAfterConverter? "yes" : "no") << G4endl
              << "    Stop beam after target? "
              << (fStopBeamAfterTarget? "yes" : "no") << G4endl;
   }
}

GlueXBernardConversionProcess::~GlueXBernardConversionProcess()
{}

G4bool GlueXBernardConversionProcess::IsApplicable(const G4ParticleDefinition& p)
{
   return (&p == G4Gamma::Gamma());
}

void GlueXBernardConversionProcess::InitialiseProcess(const G4ParticleDefinition*)
{
   if (!isInitialised) {
      isInitialised = true;
      G4EmParameters* param = G4EmParameters::Instance();
      G4double emin = std::max(param->MinKinEnergy(), 2*electron_mass_c2);
      G4double emax = param->MaxKinEnergy();

      SetMinKinEnergy(emin);

      if (!EmModel(0)) {
         SetEmModel(new G4PairProductionRelModel(), 0);
      }
      EmModel(0)->SetLowEnergyLimit(emin);
      EmModel(0)->SetHighEnergyLimit(emax);
      AddEmModel(1, EmModel(0));
   } 
}

G4double GlueXBernardConversionProcess::MinPrimaryEnergy(const G4ParticleDefinition*,
                                                         const G4Material*)
{
  return 2*electron_mass_c2;
}

void GlueXBernardConversionProcess::PrintInfo()
{}         


void GlueXBernardConversionProcess::ProcessDescription(std::ostream& out) const
{
  out << "  Bernard Gamma conversion (forced)";
  G4VEmProcess::ProcessDescription(out);
}

G4double GlueXBernardConversionProcess::GetMeanFreePath(
                                        const G4Track& /* track */, 
                                        G4double /* previousStepSize */,
                                        G4ForceCondition* /* condition */)
{
   return 100*cm;
}

G4double GlueXBernardConversionProcess::PostStepGetPhysicalInteractionLength(
                                        const G4Track &track,
                                        G4double previousStepSize,
                                        G4ForceCondition *condition)
{
   const G4Step *step = G4ParallelWorldProcess::GetHyperStep();
   G4VPhysicalVolume *pvol = step->GetPostStepPoint()->GetPhysicalVolume();
   if (track.GetTrackID() == 1 && pvol && pvol->GetName() == "PTAR" &&
       (FORCED_PTAR_PAIR_CONVERSION || 
        fStopBeamBeforeConverter ||
        fStopBeamAfterConverter ))
   {
      fPIL = G4VEmProcess::PostStepGetPhysicalInteractionLength(
                             track, previousStepSize, condition);
      *condition = Forced;
      return 100*cm;
   }
   else if (track.GetTrackID() == 1 && pvol && pvol->GetName() == "LIH2" &&
       (FORCED_LIH2_PAIR_CONVERSION || 
        fStopBeamAfterTarget ))
   {
      fPIL = G4VEmProcess::PostStepGetPhysicalInteractionLength(
                             track, previousStepSize, condition);
      *condition = Forced;
      return 100*cm;
   }
   else if (track.GetTrackID() == 1 && pvol && pvol->GetName() == "TGT0" &&
       (FORCED_TGT0_PAIR_CONVERSION || 
        fStopBeamAfterTarget ))
   {
      fPIL = G4VEmProcess::PostStepGetPhysicalInteractionLength(
                             track, previousStepSize, condition);
      *condition = Forced;
      return 100*cm;
   }
   else if (track.GetTrackID() == 1 && pvol && pvol->GetName() == "BETG" &&
       (FORCED_BETG_PAIR_CONVERSION || 
        fStopBeamAfterTarget ))
   {
      fPIL = G4VEmProcess::PostStepGetPhysicalInteractionLength(
                             track, previousStepSize, condition);
      *condition = Forced;
      return 100*cm;
   }
   else if (track.GetTrackID() == 1 && pvol && pvol->GetName() == "LIHE" &&
       (FORCED_LIHE_PAIR_CONVERSION || 
        fStopBeamAfterTarget ))
   {
      fPIL = G4VEmProcess::PostStepGetPhysicalInteractionLength(
                             track, previousStepSize, condition);
      *condition = Forced;
      return 100*cm;
   }
   *condition = NotForced;
   return 1e99;
}

G4VParticleChange *GlueXBernardConversionProcess::PostStepDoIt(
                                                  const G4Track &track, 
                                                  const G4Step &step)
{
   pParticleChange->Initialize(track);
   GlueXUserEventInformation *eventinfo;
   const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
   eventinfo = (GlueXUserEventInformation*)event->GetUserInformation();
   if (fStopBeamBeforeConverter) {
      double tvtx = step.GetPreStepPoint()->GetGlobalTime();
      G4ThreeVector vtx = step.GetPreStepPoint()->GetPosition();
      G4ThreeVector mom = step.GetPreStepPoint()->GetMomentum();
      G4ThreeVector pol = step.GetPreStepPoint()->GetPolarization();
      G4PrimaryVertex* vertex = new G4PrimaryVertex(vtx, tvtx);
      G4PrimaryParticle* photon = new G4PrimaryParticle(G4Gamma::Definition(),
                                                        mom[0], mom[1], mom[2]);
      photon->SetPolarization(pol);
      vertex->SetPrimary(photon);
      eventinfo->AddPrimaryVertex(*vertex);
      eventinfo->AddBeamParticle(1, tvtx, vtx, mom, pol);

      if (verboseLevel > 0) {
         G4cout << "GlueXBernardConversionProcess: beam particle stopped"
                << " before converter, stored in ouptut primary vertex."
                << G4endl;
      }
   }
   else if (fStopBeamAfterConverter) {
      double tvtx = step.GetPreStepPoint()->GetGlobalTime();
      G4ThreeVector vtx = step.GetPreStepPoint()->GetPosition();
      G4ThreeVector mom = step.GetPreStepPoint()->GetMomentum();
      G4ThreeVector pol = step.GetPreStepPoint()->GetPolarization();
      eventinfo->AddBeamParticle(1, tvtx, vtx, mom, pol);
      GenerateConversionVertex(track, step);

      if (verboseLevel > 0) {
         G4cout << "GlueXBernardConversionProcess: beam particle stopped"
                << " at converter exit, pair conversion forced."
                << G4endl;
      }
   }
   else if (fStopBeamAfterTarget) {
      double tvtx = step.GetPreStepPoint()->GetGlobalTime();
      G4ThreeVector vtx = step.GetPreStepPoint()->GetPosition();
      G4ThreeVector mom = step.GetPreStepPoint()->GetMomentum();
      G4ThreeVector pol = step.GetPreStepPoint()->GetPolarization();
      eventinfo->AddBeamParticle(1, tvtx, vtx, mom, pol);
      GenerateConversionVertex(track, step);
   }
   else {
      G4VEmProcess::PostStepDoIt(track, step);
      if (verboseLevel > 0) {
         G4cout << "GlueXBernardConversionProcess: beam particle stopped"
                << " unexpectedly, pair conversion forced. But WHY?"
                << G4endl;
      }
   }
   pParticleChange->ProposeTrackStatus(fStopAndKill);
   eventinfo->SetKeepEvent(1);
   return pParticleChange;
}

void GlueXBernardConversionProcess::GenerateConversionVertex(const G4Track &track,
                                                             const G4Step &step)
{
   G4ThreeVector pol = step.GetPreStepPoint()->GetPolarization();
   double pmag = pol.mag();
   if (G4UniformRand() > pmag) {
      G4Track *mtrack = (G4Track*)&track;
      mtrack->SetPolarization(G4ThreeVector(0,0,0));
   }
   else {
      G4Track *mtrack = (G4Track*)&track;
      mtrack->SetPolarization(pol / pmag); 
   }
   G4VParticleChange *pchange = G4VEmProcess::PostStepDoIt(track, step);

   // Generate a new vertex for the pair
   double beamVelocity = GlueXPhotonBeamGenerator::getBeamVelocity();
   double steplength = pParticleChange->GetTrueStepLength();
   G4ThreeVector direction(track.GetMomentumDirection());
   G4ThreeVector x0(track.GetPosition());
   double t0 = track.GetGlobalTime();
   double uvtx = G4UniformRand();
   double lvtx = steplength + fPIL * log(1 - uvtx * (1 - exp(-steplength / fPIL)));
   x0 -= lvtx * direction;
   t0 -= lvtx / beamVelocity;
   G4PrimaryVertex vertex(x0, t0);

   G4TrackVector secondaries;
   GlueXUserEventInformation *event_info;
   const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
   event_info = (GlueXUserEventInformation*)event->GetUserInformation();
   for (int is=0; is < pchange->GetNumberOfSecondaries(); ++is) {
      GlueXUserTrackInformation *trackinfo = new GlueXUserTrackInformation();
      if (event_info) {
         trackinfo->SetGlueXTrackID(event_info->AssignNextGlueXTrackID());
      }
      pchange->GetSecondary(is)->SetPosition(x0);
      pchange->GetSecondary(is)->SetGlobalTime(t0);
      pchange->GetSecondary(is)->SetUserInformation(trackinfo);
      pchange->GetSecondary(is)->SetTrackStatus(fStopAndKill);
      secondaries.push_back(pchange->GetSecondary(is));
   }

   // append secondary vertex to MC record
   if (event_info) {
      int mech[2];
      char *cmech = (char*)mech;
      snprintf(cmech, 5, "%c%c%c%c", 'B', 'E', 'R', 'N');
      event_info->AddSecondaryVertex(secondaries, 1, mech[0]);
      hddm_s::ReactionList rea = event_info->getOutputRecord()->getReactions();
      if (rea.size() > 0) {
         rea(0).setType(221); // Bethe-Heitler process
         rea(0).setWeight(1.0);
      }
   }
}

G4ParticleDefinition *GlueXBernardConversionProcess::GetLepton(int charge) {
   if (fLeptonPairFamily == 0) {
      if (charge > 0)
         return G4Positron::Definition();
      else
         return G4Electron::Definition();
   }
   else if (fLeptonPairFamily == 1) {
      if (charge > 0)
         return G4MuonPlus::Definition();
      else
         return G4MuonMinus::Definition();
   }
   return 0;
}
