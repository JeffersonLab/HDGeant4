//
// class implementation for GlueXKlongConversionProcess
//
// author: richard.t.jones at uconn.edu
// version: august 2, 2024
//

#include "GlueXKlongConversionProcess.hh"
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

// Phi photoproduction total cross section as function of photon energy,
// digitized from Fig. 3 in Wang et al, arXiv:2208.10289 (22 Aug 2022).

double gammaPhiXS_dE_GeV(0.1);
double gammaPhiXS_ub[220] = {
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.079,
0.185, 0.24, 0.291, 0.314, 0.337, 0.36, 0.383, 0.395, 0.403, 0.411,
0.419, 0.427, 0.435, 0.443, 0.451, 0.459, 0.465, 0.468, 0.47, 0.473,
0.475, 0.478, 0.48, 0.483, 0.485, 0.487, 0.49, 0.492, 0.495, 0.497,
0.5, 0.502, 0.505, 0.507, 0.51, 0.512, 0.515, 0.517, 0.519, 0.522,
0.524, 0.526, 0.528, 0.529, 0.531, 0.532, 0.533, 0.535, 0.536, 0.538,
0.539, 0.541, 0.542, 0.543, 0.545, 0.546, 0.548, 0.549, 0.551, 0.552,
0.554, 0.555, 0.556, 0.558, 0.559, 0.561, 0.562, 0.564, 0.565, 0.566,
0.568, 0.569, 0.57, 0.571, 0.572, 0.573, 0.574, 0.575, 0.576, 0.577,
0.578, 0.579, 0.58, 0.581, 0.582, 0.583, 0.584, 0.585, 0.586, 0.587,
0.588, 0.589, 0.59, 0.591, 0.592, 0.593, 0.594, 0.595, 0.596, 0.597,
0.598, 0.599, 0.6, 0.6, 0.601, 0.602, 0.603, 0.604, 0.605, 0.605,
0.606, 0.607, 0.608, 0.609, 0.61, 0.61, 0.611, 0.612, 0.613, 0.614,
0.615, 0.615, 0.616, 0.617, 0.618, 0.619, 0.62, 0.62, 0.621, 0.622,
0.623, 0.624, 0.625, 0.625, 0.626, 0.627, 0.628, 0.629, 0.629, 0.63,
0.631, 0.631, 0.632, 0.633, 0.634, 0.634, 0.635, 0.636, 0.637, 0.637,
0.638, 0.639, 0.639, 0.64, 0.641, 0.642, 0.642, 0.643, 0.644, 0.644,
0.645, 0.646, 0.647, 0.647, 0.648, 0.649, 0.65, 0.65, 0.651, 0.652,
0.652, 0.653, 0.653, 0.654, 0.654, 0.654, 0.655, 0.655, 0.655, 0.656,
0.656, 0.656, 0.656, 0.657, 0.657, 0.657, 0.658, 0.658, 0.658, 0.659,
0.659, 0.659, 0.66, 0.66, 0.66, 0.661, 0.661, 0.661, 0.662, 0.662,
};
double gammaPhiXS_tslope_GeV2(6.2);
double gammaPhiXS_Emin_GeV(1.9);

// conversion target material properties
double beryllium_A(9)
double beryllium_density_gcm3(1.85);
double tungesten_A(0.9 * 183.8 + 0.1 * 63.5);
double tungesten_density_gcm3(1 / ((0.9 / 19.38) + (0.1 / 8.96)));

G4Mutex GlueXKlongConversionProcess::fMutex = G4MUTEX_INITIALIZER;
int GlueXKLongConversionProcess::fStopBeamAfterConverter = 0;
int GlueXKlongConversionProcess::fConfigured = 0;

GlueXKlongConversionProcess::GlueXKlongConversionProcess(
                                       const G4String &name, 
                                       G4ProcessType aType)
 : G4VEmProcess(name, aType),
   isInitialised(false)
{
   SetMinKinEnergy(gammaPhiXS_Emin_GeV*GeV);
   SetProcessSubType(fGammaConversion);
   SetStartFromNullFlag(true);
   SetBuildTableFlag(true);
   SetLambdaBinning(220);

   verboseLevel = 0;

   G4AutoLock barrier(&fMutex);
   if (fConfigured)
      return;

   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXKlongConversionProcess constructor - "
             << "GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   std::map<int,std::string> infile;
   std::map<int,double> beampars;
   std::map<int,std::string> genbeampars;
   if (!user_opts->Find("INFI", infile) &&
       user_opts->Find("GENBEAM", genbeampars))
   {
      if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "KLgen" ||
               genbeampars[1] == "KLGEN" ||
               genbeampars[1] == "KL_gen" ||
               genbeampars[1] == "klong_gen" ))
      {
         fStopBeamAfterConverter = 1;
      }
   }

   fConfigured = 1;

   if (verboseLevel > 0) {
       G4cout << GetProcessName() << " is created " << G4endl
              << "    Stop beam after converter? "
              << (fStopBeamAfterConverter? "yes" : "no") << G4endl
   }
}

GlueXKlongConversionProcess::~GlueXKlongConversionProcess()
{}

G4bool GlueXKlongConversionProcess::IsApplicable(const G4ParticleDefinition& p)
{
   return (&p == G4Gamma::Gamma());
}

void GlueXKlongConversionProcess::InitialiseProcess(const G4ParticleDefinition*)
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

G4double GlueXKlongConversionProcess::MinPrimaryEnergy(const G4ParticleDefinition*,
                                                         const G4Material*)
{
  return 2*electron_mass_c2;
}

void GlueXKlongConversionProcess::PrintInfo()
{}         


void GlueXKlongConversionProcess::ProcessDescription(std::ostream& out) const
{
  out << "  Klong Gamma conversion (forced)";
  G4VEmProcess::ProcessDescription(out);
}

G4double GlueXKlongConversionProcess::GetMeanFreePath(
                                     const G4Track &track, 
                                     G4double previousStepSize,
                                     G4ForceCondition *condition)
{
   return 100*cm;
}

G4double GlueXKlongConversionProcess::PostStepGetPhysicalInteractionLength(
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

G4VParticleChange *GlueXKlongConversionProcess::PostStepDoIt(
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
         G4cout << "GlueXKlongConversionProcess: beam particle stopped"
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
         G4cout << "GlueXKlongConversionProcess: beam particle stopped"
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
         G4cout << "GlueXKlongConversionProcess: beam particle stopped"
                << " unexpectedly, pair conversion forced. But WHY?"
                << G4endl;
      }
   }
   pParticleChange->ProposeTrackStatus(fStopAndKill);
   eventinfo->SetKeepEvent(1);
   return pParticleChange;
}

void GlueXKlongConversionProcess::GenerateConversionVertex(const G4Track &track,
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
   for (int s=0; s < pchange->GetNumberOfSecondaries(); ++s) {
      GlueXUserTrackInformation *trackinfo = new GlueXUserTrackInformation();
      if (event_info) {
         trackinfo->SetGlueXTrackID(event_info->AssignNextGlueXTrackID());
      }
      pchange->GetSecondary(s)->SetPosition(x0);
      pchange->GetSecondary(s)->SetGlobalTime(t0);
      pchange->GetSecondary(s)->SetUserInformation(trackinfo);
      pchange->GetSecondary(s)->SetTrackStatus(fStopAndKill);
      secondaries.push_back(pchange->GetSecondary(s));
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

G4ParticleDefinition *GlueXKlongConversionProcess::GetLepton(int charge) {
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
