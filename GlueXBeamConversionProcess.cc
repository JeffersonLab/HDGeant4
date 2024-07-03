//
// class implementation for GlueXBeamConversionProcess
//
// author: richard.t.jones at uconn.edu
// version: december 24, 2016
//
// programmer's notes
// 1) may 5, 2021 [rtj]
//    Added support for pair conversion off the electrons in the target,
//    also called triplet production.
// 1) july 5, 2021 [rtj]
//    Added the ability to generate muon pairs instead of e+e- pairs,
//    enabled using the control.in flag 'BHmuons' instead of 'BHgen'.
// 2) july 19, 2021 [rtj]
//    I extended the functionality of the GenerateBetheHeitlerConversion
//    method to cover the case of nuclear targets, eg. 208Pb for the CPP
//    experiment. This necessarily raises the question of how to handle
//    recoil nuclear excitations. To keep things simple, I handle just
//    three discrete cases:
//        a. coherent nuclear recoil in the ground state
//        b. single proton knock-out (quasi-elastic)
//        c. single neutron knock-out (quasi-elastic)
//    Each of the three cases is weighted by the appropriate production
//    cross section for that reaction, governed by the nuclear charge
//    form factor. In the case of quasi-elastic scattering, Pauli blocking
//    by the other nucleons at low Q is taken into account in an approximate
//    way by weighting the cross section by |1 - FF(Q)^2|. Other than that,
//    the remaining A-1 nucleons are treated as spectators in the process.
//    I recognize that the quasi-deuteron process is probably just as
//    important at the relevant Q-scale as the quasi-neutron process, but
//    I am not interested in getting into that level of detail in this
//    generator.

#include "GlueXBeamConversionProcess.hh"
#include "GlueXPhotonBeamGenerator.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"
#include "GlueXUserOptions.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4PairProductionRelModel.hh"
#include "G4NistManager.hh"

#define VERBOSE_PAIRS_SPLITTING 1
#define DO_PAIRCOH_IMPORTANCE_SAMPLE 1
#define USE_ADAPTIVE_SAMPLER 1

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

// If the pair conversion target is a nucleus, the generator
// needs to know what fraction of the events to throw as
// elastic-nuclear vs quasi-elastic-nucleon, and in the case
// of quasi-elestic-nucleon, what fraction to throw as
// quasi-elastic-proton vs quasi-elastic-neutron. The output
// events are properly weighted by the respective production
// rates no matter what values are assigned here, provided
// that none of the values are 0 or 1. Adjusting these
// values enables the user to improve the statistics of one
// of these output channels at the expense of others. The
// values below should be good initial values for most targets.
// The three values should be >0 and should sum up to unity.
// If the target is hydogen, these values are ignored.
#define ELASTIC_NUCLEAR_FRACTION 0.75
#define QUASIELASTIC_PROTON_FRACTION 0.20
#define QUASIELASTIC_NEUTRON_FRACTION 0.05

#define AMU_GEV 0.9313680888469

#include <G4SystemOfUnits.hh>
#include "G4PhysicalConstants.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4RunManager.hh"
#include "G4TrackVector.hh"

#include <stdio.h>
#include <iomanip>

#ifdef USING_DIRACXX
#include <TLorentzBoost.h>
#include <TPhoton.h>
#include <TLepton.h>
#include <TCrossSection.h>

#else
const double mElectron=0.51099907e-3;
const double mMuon=0.105658389;
const double mProton=0.9382723128;
const double mNeutron=0.9395656328;
#endif

void unif01(int n, double *u) {
   G4Random::getTheEngine()->flatArray(n,u);
}

G4double GlueXBeamConversionProcess::fBHpair_mass_min = 0;

std::vector<AdaptiveSampler*> GlueXBeamConversionProcess::fAdaptiveSamplerRegistry;

G4Mutex GlueXBeamConversionProcess::fMutex = G4MUTEX_INITIALIZER;
int GlueXBeamConversionProcess::fStopBeamBeforeConverter = 0;
int GlueXBeamConversionProcess::fStopBeamAfterConverter = 0;
int GlueXBeamConversionProcess::fStopBeamAfterTarget = 0;
int GlueXBeamConversionProcess::fLeptonPairFamily = 0;
int GlueXBeamConversionProcess::fConfigured = 0;

int GlueXBeamConversionProcess::fFormFactorChoice = 0;

GlueXBeamConversionProcess::GlueXBeamConversionProcess(const G4String &name, 
                                                       G4ProcessType aType)
 : G4VEmProcess(name, aType),
#  ifdef USING_DIRACXX
   fPairsGeneration(0),
#  endif
   fPaircohPDF(0),
   fTripletPDF(0),
   fAdaptiveSampler(0),
   isInitialised(false),
   fTargetZ(0),
   fTargetA(0)
{
   verboseLevel = 0;

   G4AutoLock barrier(&fMutex);
   if (fConfigured)
      return;

   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXBeamConversionProcess constructor - "
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
         (genbeampars[1] == "POSTCOL" ||
          genbeampars[1] == "postcol" ||
          genbeampars[1] == "Postcol" ||
          genbeampars[1] == "PostCol" ))
      {
         fStopBeamBeforeConverter = 1;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "POSTCONV" ||
               genbeampars[1] == "postconv" ||
               genbeampars[1] == "Postconv" ||
               genbeampars[1] == "PostConv" ))
      {
         fStopBeamAfterConverter = 1;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHgen" ||
               genbeampars[1] == "BHGEN" ||
               genbeampars[1] == "BHGen" ||
               genbeampars[1] == "bhgen" ))
      {
         fStopBeamAfterTarget = 1;
         if (genbeampars.find(2) != genbeampars.end()) {
            fBHpair_mass_min = std::atof(genbeampars[2].c_str());
         }
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHgen_A" ||
               genbeampars[1] == "BHGEN_A" ||
               genbeampars[1] == "BHGen_A" ||
               genbeampars[1] == "bhgen_a" ))
      {
         fStopBeamAfterTarget = 1;
         if (genbeampars.find(2) != genbeampars.end()) {
            fBHpair_mass_min = std::atof(genbeampars[2].c_str());
         }
         fFormFactorChoice = 1;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHgen_C" ||
               genbeampars[1] == "BHGEN_C" ||
               genbeampars[1] == "BHGen_C" ||
               genbeampars[1] == "bhgen_c" ))
      {
         fStopBeamAfterTarget = 1;
         if (genbeampars.find(2) != genbeampars.end()) {
            fBHpair_mass_min = std::atof(genbeampars[2].c_str());
         }
         fFormFactorChoice = 3;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHgen_D" ||
               genbeampars[1] == "BHGEN_D" ||
               genbeampars[1] == "BHGen_D" ||
               genbeampars[1] == "bhgen_d" ))
      {
         fStopBeamAfterTarget = 1;
         if (genbeampars.find(2) != genbeampars.end()) {
            fBHpair_mass_min = std::atof(genbeampars[2].c_str());
         }
         fFormFactorChoice = 4;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHgen_E" ||
               genbeampars[1] == "BHGEN_E" ||
               genbeampars[1] == "BHGen_E" ||
               genbeampars[1] == "bhgen_e" ))
      {
         fStopBeamAfterTarget = 1;
         if (genbeampars.find(2) != genbeampars.end()) {
            fBHpair_mass_min = std::atof(genbeampars[2].c_str());
         }
         fFormFactorChoice = 5;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHgen_F" ||
               genbeampars[1] == "BHGEN_F" ||
               genbeampars[1] == "BHGen_F" ||
               genbeampars[1] == "bhgen_f" ))
      {
         fStopBeamAfterTarget = 1;
         if (genbeampars.find(2) != genbeampars.end()) {
            fBHpair_mass_min = std::atof(genbeampars[2].c_str());
         }
         fFormFactorChoice = 6;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHgen_G" ||
               genbeampars[1] == "BHGEN_G" ||
               genbeampars[1] == "BHGen_G" ||
               genbeampars[1] == "bhgen_g" ))
      {
         fStopBeamAfterTarget = 1;
         if (genbeampars.find(2) != genbeampars.end()) {
            fBHpair_mass_min = std::atof(genbeampars[2].c_str());
         }
         fFormFactorChoice = 7;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHmuons" ||
               genbeampars[1] == "BHmuons" ||
               genbeampars[1] == "BHmuons" ||
               genbeampars[1] == "bhmuons" ))
      {
         fStopBeamAfterTarget = 1;
         fLeptonPairFamily = 1;
         if (genbeampars.find(2) != genbeampars.end()) {
            fBHpair_mass_min = std::atof(genbeampars[2].c_str());
         }
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHmuons_A" ||
               genbeampars[1] == "BHMUONS_A" ||
               genbeampars[1] == "BHMuons_A" ||
               genbeampars[1] == "bhmuons_a" ))
      {
         fStopBeamAfterTarget = 1;
         fLeptonPairFamily = 1;
         if (genbeampars.find(2) != genbeampars.end()) {
            fBHpair_mass_min = std::atof(genbeampars[2].c_str());
         }
         fFormFactorChoice = 1;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHmuons_C" ||
               genbeampars[1] == "BHMUONS_C" ||
               genbeampars[1] == "BHMuons_C" ||
               genbeampars[1] == "bhmuons_c" ))
      {
         fStopBeamAfterTarget = 1;
         fLeptonPairFamily = 1;
         if (genbeampars.find(2) != genbeampars.end()) {
            fBHpair_mass_min = std::atof(genbeampars[2].c_str());
         }
         fFormFactorChoice = 3;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHmuons_D" ||
               genbeampars[1] == "BHMUONS_D" ||
               genbeampars[1] == "BHMuons_D" ||
               genbeampars[1] == "bhmuons_d" ))
      {
         fStopBeamAfterTarget = 1;
         fLeptonPairFamily = 1;
         if (genbeampars.find(2) != genbeampars.end()) {
            fBHpair_mass_min = std::atof(genbeampars[2].c_str());
         }
         fFormFactorChoice = 4;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHmuons_E" ||
               genbeampars[1] == "BHMUONS_E" ||
               genbeampars[1] == "BHMuons_E" ||
               genbeampars[1] == "bhmuons_e" ))
      {
         fStopBeamAfterTarget = 1;
         fLeptonPairFamily = 1;
         if (genbeampars.find(2) != genbeampars.end()) {
            fBHpair_mass_min = std::atof(genbeampars[2].c_str());
         }
         fFormFactorChoice = 5;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHmuons_F" ||
               genbeampars[1] == "BHMUONS_F" ||
               genbeampars[1] == "BHMuons_F" ||
               genbeampars[1] == "bhmuons_f" ))
      {
         fStopBeamAfterTarget = 1;
         fLeptonPairFamily = 1;
         if (genbeampars.find(2) != genbeampars.end()) {
            fBHpair_mass_min = std::atof(genbeampars[2].c_str());
         }
         fFormFactorChoice = 6;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHmuons_G" ||
               genbeampars[1] == "BHMUONS_G" ||
               genbeampars[1] == "BHMuons_G" ||
               genbeampars[1] == "bhmuons_g" ))
      {
         fStopBeamAfterTarget = 1;
         fLeptonPairFamily = 1;
         if (genbeampars.find(2) != genbeampars.end()) {
            fBHpair_mass_min = std::atof(genbeampars[2].c_str());
         }
         fFormFactorChoice = 7;
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

GlueXBeamConversionProcess::~GlueXBeamConversionProcess()
{
   if (fPaircohPDF)
      delete fPaircohPDF;
   if (fTripletPDF)
      delete fTripletPDF;
#ifdef USING_DIRACXX
   if (fPairsGeneration)
      delete fPairsGeneration;
#endif
   G4AutoLock barrier(&fMutex);
   int nsamplers = fAdaptiveSamplerRegistry.size();
   if (nsamplers > 0) {
      fAdaptiveSampler = fAdaptiveSamplerRegistry[0];
      for (int i=nsamplers-1; i > 0; --i) {
         std::stringstream astatefile;
         astatefile << "BHgen_thread_" << i << ".astate";
         fAdaptiveSamplerRegistry[i]->saveState(astatefile.str());
         fAdaptiveSampler->mergeState(astatefile.str());
         delete fAdaptiveSamplerRegistry[i];
         fAdaptiveSamplerRegistry.erase(fAdaptiveSamplerRegistry.begin() + i);
      }
      fAdaptiveSamplerRegistry[0]->saveState("BHgen_stats.astate");
      delete fAdaptiveSamplerRegistry[0];
      fAdaptiveSamplerRegistry.clear();
   }
}

G4bool GlueXBeamConversionProcess::IsApplicable(const G4ParticleDefinition& p)
{
   return (&p == G4Gamma::Gamma());
}

void GlueXBeamConversionProcess::InitialiseProcess(const G4ParticleDefinition*)
{
   if (!isInitialised) {
      isInitialised = true;
      G4EmParameters* param = G4EmParameters::Instance();
      G4double emin = std::max(param->MinKinEnergy(), 2*electron_mass_c2);
      if (fLeptonPairFamily == 1)
         emin = std::max(emin, 2*G4MuonMinus::Definition()->GetPDGMass());
      SetMinKinEnergy(emin);

      if (!EmModel(0)) {
         SetEmModel(new G4PairProductionRelModel(), 0);
      }
      EmModel(0)->SetLowEnergyLimit(emin);
      AddEmModel(1, EmModel(0));

#ifdef USING_DIRACXX
      std::vector<double> Z;
      std::vector<double> A;
      std::vector<double> w;
      Z.push_back(4);    // assumes 4Be converter target for TPOL
      A.push_back(9.012);
      w.push_back(1);
      fPairsGeneration = new PairConversionGeneration(Z, A, w);

#if defined DO_TRIPLET_IMPORTANCE_SAMPLE || defined DO_PAIRCOH_IMPORTANCE_SAMPLE
      if (fTripletPDF == 0 || fPaircohPDF == 0) {
         G4cout << "GlueXBeamConversionProcess::InitialiseProcess:"
                << G4endl
                << "   Setting up cross section tables, please wait... "
                << std::flush;
         prepareImportanceSamplingPDFs();
         G4cout << "finished." << G4endl;
      }
#endif

#if USE_ADAPTIVE_SAMPLER
      if (fAdaptiveSampler == 0) {
         prepareAdaptiveSampler();
      }
#endif

#endif
   }
}

G4double GlueXBeamConversionProcess::MinPrimaryEnergy(const G4ParticleDefinition*,
                                                      const G4Material*)
{
  if (fLeptonPairFamily == 1)
      return 2 * G4MuonMinus::MuonMinus()->GetPDGMass();
  return 2 * electron_mass_c2;
}

void GlueXBeamConversionProcess::PrintInfo()
{}         


void GlueXBeamConversionProcess::ProcessDescription(std::ostream& out) const
{
  out << "  Dirac++ Gamma conversion (forced)";
  G4VEmProcess::ProcessDescription(out);
}


G4double GlueXBeamConversionProcess::GetMeanFreePath(
                                     const G4Track &track, 
                                     G4double previousStepSize,
                                     G4ForceCondition *condition)
{
   return 100*cm;
}

G4double GlueXBeamConversionProcess::PostStepGetPhysicalInteractionLength(
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
      setConverterMaterial(4, 9);
      *condition = Forced;
      return 100*cm;
   }
   else if (track.GetTrackID() == 1 && pvol && pvol->GetName() == "LIH2" &&
       (FORCED_LIH2_PAIR_CONVERSION || 
        fStopBeamAfterTarget ))
   {
      fPIL = G4VEmProcess::PostStepGetPhysicalInteractionLength(
                             track, previousStepSize, condition);
      //setConverterMaterial(1, 1);
      setConverterMaterial(4, 9);
      *condition = Forced;
      return 100*cm;
   }
   else if (track.GetTrackID() == 1 && pvol && pvol->GetName() == "TGT0" &&
       (FORCED_TGT0_PAIR_CONVERSION || 
        fStopBeamAfterTarget ))
   {
      fPIL = G4VEmProcess::PostStepGetPhysicalInteractionLength(
                             track, previousStepSize, condition);
      setConverterMaterial(82, 208);
      *condition = Forced;
      return 100*cm;
   }
   else if (track.GetTrackID() == 1 && pvol && pvol->GetName() == "BETG" &&
       (FORCED_BETG_PAIR_CONVERSION || 
        fStopBeamAfterTarget ))
   {
      fPIL = G4VEmProcess::PostStepGetPhysicalInteractionLength(
                             track, previousStepSize, condition);
      setConverterMaterial(4, 9);
      *condition = Forced;
      return 100*cm;
   }
   else if (track.GetTrackID() == 1 && pvol && pvol->GetName() == "LIHE" &&
       (FORCED_LIHE_PAIR_CONVERSION || 
        fStopBeamAfterTarget ))
   {
      fPIL = G4VEmProcess::PostStepGetPhysicalInteractionLength(
                             track, previousStepSize, condition);
      setConverterMaterial(2, 4);
      *condition = Forced;
      return 100*cm;
   }
   *condition = NotForced;
   return 1e99;
}

G4VParticleChange *GlueXBeamConversionProcess::PostStepDoIt(
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
         G4cout << "GlueXBeamConversionProcess: beam particle stopped"
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
      GenerateBeamPairConversion(step);

      if (verboseLevel > 0) {
         G4cout << "GlueXBeamConversionProcess: beam particle stopped"
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
      double maxTripletMass2 = 2 * 0.511*MeV *
                               step.GetPreStepPoint()->GetKineticEnergy();
      double BetheHeitler_fraction = fTargetZ / (fTargetZ + 1);
      if (fBHpair_mass_min*GeV > sqrt(maxTripletMass2))
         BetheHeitler_fraction = 1;
      double weight_factor;
      if (G4UniformRand() < BetheHeitler_fraction) {
         if (verboseLevel > 0) {
            G4cout << "GlueXBeamConversionProcess: beam particle stopped"
                   << " in the LiH2 target, Bethe-Heitler conversion forced."
                   << G4endl;
         }
         GenerateBetheHeitlerProcess(step);
         weight_factor = 1 / BetheHeitler_fraction;
      }
      else {
         if (verboseLevel > 0) {
            G4cout << "GlueXBeamConversionProcess: beam particle stopped"
                   << " in the LiH2 target, triplet conversion forced."
                   << G4endl;
         }
         GenerateTripletProcess(step);
         weight_factor = 1 / (1. - BetheHeitler_fraction);
      }
      GlueXUserEventInformation *event_info;
      const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
      event_info = (GlueXUserEventInformation*)event->GetUserInformation();
      hddm_s::ReactionList rea = event_info->getOutputRecord()->getReactions();
      rea(0).setWeight(rea(0).getWeight() * weight_factor);
   }
   else {
      GenerateBeamPairConversion(step);

      if (verboseLevel > 0) {
         G4cout << "GlueXBeamConversionProcess: beam particle stopped"
                << " unexpectedly, pair conversion forced. But WHY?"
                << G4endl;
      }
   }
   pParticleChange->ProposeTrackStatus(fStopAndKill);
   eventinfo->SetKeepEvent(1);
   return pParticleChange;
}

void GlueXBeamConversionProcess::prepareImportanceSamplingPDFs()
{
   // Construct lookup tables representing the PDFs used for
   // importance-sampling the gamma pair conversion kinematics.

   // Compute 2D histogram containing rate as a function of u0,u1
   // where u0 is the (originally uniform [0,1]) random number used
   // to generate Mpair and u1 generates qR. The algorithm assumes
   // that the mapping u0->Mpair and u1->qR used here is the same
   // as the mmpaing used in the process generation methods, currently
   // GenerateBeamPairConversion and GenerateTargetConversionProcess.

   if (fPaircohPDF == 0 || fTripletPDF == 0) {
      fPaircohPDF = new ImportanceSampler;
      fTripletPDF = new ImportanceSampler;
   }

   fPaircohPDF->Pcut = 60;
   fTripletPDF->Pcut = 15;

#if USING_DIRACXX

   LDouble_t kin = 9.; // GeV
   LDouble_t mLepton = mElectron;
   if (fLeptonPairFamily == 1)
      mLepton = mMuon;
   TPhoton g0;
   TLepton p1(mLepton);
   TLepton e2(mLepton);
   TLepton e3(mLepton);
   TThreeVectorReal p;
   g0.SetMom(p.SetPolar(kin,0,0));
   g0.SetPol(TThreeVector(0,0,0));
   p1.AllPol();
   e2.AllPol();
   e3.AllPol();

   LDouble_t Epos = kin / 2;
   LDouble_t Mmin = 2 * mLepton;
   LDouble_t Mcut = 10 * mLepton; // 5 MeV cutoff parameter
   LDouble_t um0 = 1 + sqr(Mcut / Mmin);
   LDouble_t qRcut = 2 * mLepton; // 1 MeV/c cutoff parameter

   int Nbins = 50;
   fTripletPDF->Psum = 0;
   fPaircohPDF->Psum = 0;
   for (int i0=0; i0 < Nbins; ++i0) {
      LDouble_t u0 = (i0 + 0.5) / Nbins;
      LDouble_t um = pow(um0, u0);
      LDouble_t Mpair = Mcut / sqrt(um - 1);
      LDouble_t k12star2 = sqr(Mpair / 2) - sqr(mLepton);
      LDouble_t qRmin = sqr(Mpair) / (2 * kin);
      LDouble_t uq0 = qRmin / (qRcut + sqrt(sqr(qRcut) + sqr(qRmin)));
      LDouble_t weight0 = sqr(Mpair) * (sqr(Mcut) + sqr(Mpair));
      for (int i1=0; i1 < Nbins; ++i1) {
         LDouble_t u1 = (i1 + 0.5) / Nbins;
         LDouble_t uq = pow(uq0, u1);
         LDouble_t qR = 2 * qRcut * uq / (1 - sqr(uq));
         LDouble_t weight = weight0 * sqr(qR) * sqrt(sqr(qRcut) + sqr(qR));
         LDouble_t E3 = sqrt(sqr(qR) + sqr(mElectron));
         LDouble_t E12 = kin + mElectron - E3;
         if (k12star2 < 0 || E12 < Mpair) {
            fPaircohPDF->density.push_back(0);
            fTripletPDF->density.push_back(0);
            fPaircohPDF->integral.push_back(fPaircohPDF->Psum);
            fTripletPDF->integral.push_back(fTripletPDF->Psum);
            continue;
         }
         LDouble_t k12star = sqrt(k12star2);
         LDouble_t q12mag = sqrt(sqr(E12) - sqr(Mpair));
         LDouble_t costhetastar = (Epos - E12 / 2) * 
                                  Mpair / (k12star * q12mag);
         LDouble_t costhetaR = (sqr(Mpair) / 2 + (kin + mElectron) *
                                (E3 - mElectron)) / (kin * qR);
         if (fabs(costhetastar) > 1 || fabs(costhetaR) > 1) {
            fPaircohPDF->density.push_back(0);
            fTripletPDF->density.push_back(0);
            fPaircohPDF->integral.push_back(fPaircohPDF->Psum);
            fTripletPDF->integral.push_back(fTripletPDF->Psum);
            continue;
         }
         LDouble_t qRlong = qR * costhetaR;
         LDouble_t qRperp = sqrt(sqr(qR) - sqr(qRlong));
         TThreeVectorReal q3(0, qRperp, qRlong);
         LDouble_t sinthetastar = sqrt(1 - sqr(costhetastar));
         TThreeVectorReal k12(k12star * sinthetastar, 0, 
                              k12star * costhetastar);
         TFourVectorReal q1(Mpair / 2, k12);
         TFourVectorReal q2(Mpair / 2, -k12);
         TLorentzBoost tolab(q3[1] / E12, q3[2] / E12, (q3[3] - kin) / E12);
         q1.Boost(tolab);
         q2.Boost(tolab);
         p1.SetMom(q1);
         e2.SetMom(q2);
         e3.SetMom(q3);
         LDouble_t tripXS = fPairsGeneration->DiffXS_triplet(g0,p1,e2,e3);
         LDouble_t pairXS = fPairsGeneration->DiffXS_pair(g0,p1,e2);
         fTripletPDF->Psum += fTripletPDF->Pmax = tripXS * weight;
         fPaircohPDF->Psum += fPaircohPDF->Pmax = pairXS * weight;
         fTripletPDF->density.push_back(fTripletPDF->Pmax);
         fPaircohPDF->density.push_back(fPaircohPDF->Pmax);
         fPaircohPDF->integral.push_back(fPaircohPDF->Psum);
         fTripletPDF->integral.push_back(fTripletPDF->Psum);
      }
   }

   LDouble_t du2 = 1. / sqr(Nbins);
   for (int i0=0, index=0; i0 < Nbins; ++i0) {
      for (int i1=0; i1 < Nbins; ++i1, ++index) {
         LDouble_t randvar = i0 + (i1 + 0.5) / Nbins;
         fTripletPDF->randvar.push_back(randvar);
         fPaircohPDF->randvar.push_back(randvar);
         fTripletPDF->density[index] /= fTripletPDF->Psum * du2;
         fPaircohPDF->density[index] /= fPaircohPDF->Psum * du2;
         fTripletPDF->integral[index] /= fTripletPDF->Psum;
         fPaircohPDF->integral[index] /= fPaircohPDF->Psum;
      }
   }
#endif

   fTripletPDF->Pmax = 0;
   fTripletPDF->Psum = 0;
   fPaircohPDF->Pmax = 0;
   fPaircohPDF->Psum = 0;
}

void GlueXBeamConversionProcess::GenerateBeamPairConversion(const G4Step &step)
{
   // Unlike the other GenerateXXX methods in this class, this method should
   // be invoked after tracking of an event is already underway. Its purpose
   // is to generate pair conversion / triplet production by gamma rays in a
   // converter target with the full polarization-dependent QED differential
   // cross section. Both incoherent (triplet) and nuclear + atomic coherent
   // scattering is simulated by this method. It should be called from your
   // G4SteppingAction::UserSteppingAction method to force pair conversion
   // in special simulations dedicated to the study of the analyzing power
   // of this reaction. The incident gamma ray is stopped and the pair (plus
   // recoil electron, if any) are added to the event stack as new primary
   // particles.

#ifndef USING_DIRACXX

   G4cerr << "GlueXBeamConversionProcess::GenerateBeamPairConversion error:"
          << G4endl
          << "  You have enabled pair/triplet conversion in the PTAR target,"
          << G4endl
          << "  but your have built HDGeant4 without support for the Dirac++"
		  << G4endl
          << "  library. Either rebuild with Dirac++ support or else jack up"
		  << G4endl
          << "  the value of constant FORCED_PTAR_PAIR_CONVERSION_THRESHOLD"
		  << G4endl
          << "  in GlueXSteppingAction.cc and try again. Aborting this run..."
		  << G4endl;
   exit(1);

#else

   LDouble_t mLepton = mElectron;
   if (fLeptonPairFamily == 1)
      mLepton = mMuon;
   TPhoton gIn;
   TLepton p1(mLepton);
   TLepton e2(mLepton);
   TLepton e3(mElectron);
   p1.AllPol();
   e2.AllPol();
   e3.AllPol();
   const G4Track *track = step.GetTrack();
   LDouble_t kin = track->GetKineticEnergy()/GeV;

   // If we are below pair production threshold, do nothing
   const double mTarget = step.GetPreStepPoint()->GetMaterial()->GetA() * 0.932;
   if (kin < 2 * mLepton * (1 + mLepton / mTarget))
      return;

   G4ThreeVector mom(track->GetMomentum());
   TThreeVectorReal mom0(mom[0]/GeV, mom[1]/GeV, mom[2]/GeV);
   gIn.SetMom(mom0);
   G4ThreeVector pol(track->GetPolarization());
   TThreeVectorReal pol0(pol[0], pol[1], pol[2]);
   if (pol0.Length() > 0)
      gIn.SetPlanePolarization(pol0, pol0.Length());
   else
      gIn.SetPol(pol0);

   // Define an angle and axis that rotates zhat into the direction
   // of the incident gamma, so that the generated kinematics is
   // defined with the incident gamma aligned with zhat, and then
   // rotated at the end into the final spatial direction.
   TThreeVectorReal rockaxis(mom0);
   rockaxis.Cross(TThreeVectorReal(0,0,1));
   double rockangle(0);
   if (mom0.Length() > 0) {
      rockangle = asin(rockaxis.Length() / mom0.Length());
      rockaxis /= rockaxis.Length();
   }

   while (true) {
      LDouble_t weight = 1;

      // Generate uniform in E+, phi12, phiR
      LDouble_t Epos = kin * G4UniformRand();
      while (Epos < mLepton) {
         Epos = kin * G4UniformRand();
      }
      weight *= kin - mLepton;
      LDouble_t phi12 = 2*M_PI * G4UniformRand();
      weight *= 2*M_PI;
      LDouble_t phiR = 2*M_PI * G4UniformRand();
      weight *= 2*M_PI;

      LDouble_t u0 = G4UniformRand();
      LDouble_t u1 = G4UniformRand();

#if DO_TRIPLET_IMPORTANCE_SAMPLE

      int i = fTripletPDF->search(u1);
      LDouble_t fi = fTripletPDF->density[i];
      LDouble_t ui = fTripletPDF->integral[i];
      LDouble_t ri = fTripletPDF->randvar[i];
      LDouble_t xi = ri - floor(ri);
      LDouble_t dx = (xi > 0.5)? ri - fTripletPDF->randvar[i-1]:
                              fTripletPDF->randvar[i+1] - ri;
      u1 = xi + dx / 2 - (ui - u1) / (fi * dx);
      u0 = (u0 + floor(ri)) * dx;
      weight /= fi;

#elif DO_PAIRCOH_IMPORTANCE_SAMPLE

      int i = fPaircohPDF->search(u1);
      LDouble_t fi = fPaircohPDF->density[i];
      LDouble_t ui = fPaircohPDF->integral[i];
      LDouble_t ri = fPaircohPDF->randvar[i];
      LDouble_t xi = ri - floor(ri);
      LDouble_t dx = (xi > 0.5)? ri - fPaircohPDF->randvar[i-1]:
                              fPaircohPDF->randvar[i+1] - ri;
      u1 = xi + dx / 2 - (ui - u1) / (fi * dx);
      u0 = (u0 + floor(ri)) * dx;
      weight /= fi;

#endif
   
      // Generate Mpair as 1 / (M [M^2 + Mcut^2])
      LDouble_t Mmin = 2 * mLepton;
      LDouble_t Mcut = 10 * mLepton;  // 5 MeV for electrons
      LDouble_t um0 = 1 + sqr(Mcut / Mmin);
      LDouble_t um = pow(um0, u0);
      LDouble_t Mpair = Mcut / sqrt(um - 1 + 1e-99);
      weight *= Mpair * (sqr(Mcut) + sqr(Mpair)) * log(um0) / (2 * sqr(Mcut));
   
      // Generate qR^2 with weight 1 / [qR^2 sqrt(qRcut^2 + qR^2)]
      LDouble_t qRmin = sqr(Mpair) /(2 * kin);
      LDouble_t qRcut = 2 * mLepton; // 1 MeV for electrons
      LDouble_t uq0 = qRmin / (qRcut + sqrt(sqr(qRcut) + sqr(qRmin)));
      LDouble_t uq = pow(uq0, u1);
      LDouble_t qR = 2 * qRcut * uq / (1 - sqr(uq));
      LDouble_t qR2 = qR * qR;
      weight *= qR2 * sqrt(1 + qR2 / sqr(qRcut)) * (-2 * log(uq0));
   
      // Include overall measure Jacobian factor
      weight *= Mpair / (2 * kin);
   
      // Generate with importance sampling
      LDouble_t Striplet = fTripletPDF->Npassed * 
                        (fTripletPDF->Ntested / (fTripletPDF->Psum + 1e-99));
      LDouble_t Spaircoh = fPaircohPDF->Npassed *
                        (fPaircohPDF->Ntested / (fPaircohPDF->Psum + 1e-99));
      if (Striplet < Spaircoh) {                     // try incoherent generation
         ++fTripletPDF->Ntested;
   
         // Solve for the c.m. momentum of e+ in the pair 1,2 rest frame
         LDouble_t k12star2 = sqr(Mpair / 2) - sqr(mLepton);
         if (k12star2 < 0) {
            // try again (this should never happen!)
            continue;
         }
         LDouble_t k12star = sqrt(k12star2);
         LDouble_t E3 = sqrt(qR2 + sqr(mElectron));
         LDouble_t E12 = kin + mElectron - E3;
         if (E12 < Mpair) {
            // no kinematic solution because E12 < Mpair, try again
            continue;
         }
         LDouble_t q12mag = sqrt(sqr(E12) - sqr(Mpair));
         LDouble_t costhetastar = (Epos - E12 / 2) * Mpair / (k12star * q12mag);
         if (Epos > E12 - mLepton) {
            // no kinematic solution because Epos > E12 - mLepton, try again
            continue;
         }
         else if (fabs(costhetastar) > 1) {
            // no kinematic solution because |costhetastar| > 1, try again
            continue;
         }
   
         // Solve for the recoil electron kinematics
         LDouble_t costhetaR = (sqr(Mpair) / 2 + (kin + mElectron) *
                                (E3 - mElectron)) / (kin * qR);
         if (fabs(costhetaR) > 1) {
            // no kinematic solution because |costhetaR| > 1, try again
            continue;
         }
         LDouble_t sinthetaR = sqrt(1 - sqr(costhetaR));
         TFourVectorReal q3(E3, qR * sinthetaR * cos(phiR),
                                qR * sinthetaR * sin(phiR),
                                qR * costhetaR);
   
         // Boost the pair momenta into the lab
         LDouble_t sinthetastar = sqrt(1 - sqr(costhetastar));
         TThreeVectorReal k12(k12star * sinthetastar * cos(phi12),
                              k12star * sinthetastar * sin(phi12),
                              k12star * costhetastar);
         TFourVectorReal q1(Mpair / 2, k12);
         TFourVectorReal q2(Mpair / 2, -k12);
         TLorentzBoost toLab(q3[1] / E12, q3[2] / E12, (q3[3] - kin) / E12);
         q1.Boost(toLab);
         q2.Boost(toLab);
   
         // To avoid double-counting, return zero if recoil electron
         // momentum is greater than the momentum of the pair electron.
         if (fLeptonPairFamily == 0 && q2.Length() < qR) {
            // recoil/pair electrons switched, try again
            continue;
         }

         // Rotate final-state momenta into the lab frame
         if (rockangle != 0) {
            q1.Rotate(rockaxis, rockangle);
            q2.Rotate(rockaxis, rockangle);
            q3.Rotate(rockaxis, rockangle);
         }
         p1.SetMom(q1);
         e2.SetMom(q2);
         e3.SetMom(q3);

         // Check 4-momentum conservation
         TFourVectorReal pIn, pFi;
         TFourVectorReal::SetResolution(1e-10);
         pIn = gIn.Mom() + TFourVectorReal(mElectron,0,0,0);
         pFi = p1.Mom() + e2.Mom() + e3.Mom();
         if (pIn != pFi) {
            std::cout << "Warning in GenerateBeamPairConversion - "
                      << "momentum conservation violated." << std::endl
                      << "   pIn = "; 
            pIn.Print();
            std::cout << "   pFi = ";
            pFi.Print();
            std::cout << "   pIn - pFi = ";
            (pIn-pFi).Print();
         }

         // Compute the differential cross section (barns/GeV^4)
         // returned as d(sigma)/(dE+ dphi+ d^3qR)
         LDouble_t diffXS = fPairsGeneration->DiffXS_triplet(gIn, p1, e2, e3);

         // Use keep/discard algorithm
         LDouble_t Pfactor = diffXS * weight;
         if (Pfactor > fTripletPDF->Pmax)
            fTripletPDF->Pmax = Pfactor;
         if (Pfactor > fTripletPDF->Pcut) {
            G4cout << "Warning in GenerateBeamPairConversion - Pfactor " 
                   << Pfactor << " exceeds fTripletPDF->Pcut = " 
                   << fTripletPDF->Pcut << G4endl
                   << "  present qR = " << qR << G4endl
                   << "  present Mpair = " << Mpair << G4endl
                   << "  present maximum Pfactor = "
                   << fTripletPDF->Pmax << G4endl
                   << "  current generator efficiency = "
                   << fTripletPDF->Npassed /
                      (fTripletPDF->Ntested + 1e-99)
                   << G4endl;
         }
         fTripletPDF->Psum += Pfactor;
         if (G4UniformRand() * fTripletPDF->Pcut > Pfactor) {
            continue;
         }
         ++fTripletPDF->Npassed;
         break;
      }

      else {                          // try coherent generation
         ++fPaircohPDF->Ntested;
   
         // Solve for the c.m. momentum of e+ in the pair 1,2 rest frame
         // assuming that the atomic target absorbs zero energy
         LDouble_t k12star2 = sqr(Mpair / 2) - sqr(mLepton);
         if (k12star2 < 0) {
            // try again (this should never happen!)
            continue;
         }
         LDouble_t k12star = sqrt(k12star2);
         LDouble_t Eneg = kin - Epos;
         if (kin < Mpair) {
            // no kinematic solution because kin < Mpair, try again
            continue;
         }
         else if (Eneg < mLepton) {
            // no kinematic solution because Eneg < mLepton, try again
            continue;
         }
         LDouble_t q12mag = sqrt(sqr(kin) - sqr(Mpair));
         LDouble_t costhetastar = (Epos - kin / 2) * Mpair / (k12star * q12mag);
         if (fabs(costhetastar) > 1) {
            // no kinematic solution because |costhetastar| > 1, try again
            continue;
         }
   
         // Solve for the recoil kinematics kinematics
         LDouble_t costhetaR = (sqr(Mpair) + qR2) / (2 * kin * qR);
         if (fabs(costhetaR) > 1) {
            // no kinematic solution because |costhetaR| > 1, try again
            continue;
         }
         LDouble_t sinthetaR = sqrt(1 - sqr(costhetaR));
         TThreeVectorReal q3(qR * sinthetaR * cos(phiR),
                             qR * sinthetaR * sin(phiR),
                             qR * costhetaR);
   
         // Boost the pair momenta into the lab
         LDouble_t sinthetastar = sqrt(1 - sqr(costhetastar));
         TThreeVectorReal k12(k12star * sinthetastar * cos(phi12),
                              k12star * sinthetastar * sin(phi12),
                              k12star * costhetastar);
         TFourVectorReal q1(Mpair / 2, k12);
         TFourVectorReal q2(Mpair / 2, -k12);
         TLorentzBoost toLab(q3[1] / kin, q3[2] / kin, (q3[3] - kin) / kin);
         q1.Boost(toLab);
         q2.Boost(toLab);
   
         // Compute the differential cross section (barnes/GeV^4)
         // returned as d(sigma)/(dE+ dphi+ d^3qR)
         if (rockangle != 0) {
            q1.Rotate(rockaxis, rockangle);
            q2.Rotate(rockaxis, rockangle);
         }
         p1.SetMom(q1);
         e2.SetMom(q2);
         e3.SetMom(TThreeVectorReal(0,0,0));
         LDouble_t diffXS = fPairsGeneration->DiffXS_pair(gIn, p1, e2);
   
         // Use keep/discard algorithm
         LDouble_t Pfactor = diffXS * weight;
         if (Pfactor > fPaircohPDF->Pmax)
            fPaircohPDF->Pmax = Pfactor;
         if (Pfactor > fPaircohPDF->Pcut) {
            G4cout << "Warning in GenerateBeamPairConversion - Pfactor " 
                   << Pfactor << " exceeds fPaircohPDF->Pcut = " 
                   << fPaircohPDF->Pcut << G4endl
                   << "  present qR = " << qR << G4endl
                   << "  present Mpair = " << Mpair << G4endl
                   << "  present maximum Pfactor = "
                   << fPaircohPDF->Pmax << G4endl
                   << "  current generator efficiency = "
                   << fPaircohPDF->Npassed /
                      (fPaircohPDF->Ntested + 1e-99)
                   << G4endl;
         }
         fPaircohPDF->Psum += Pfactor;
         if (G4UniformRand() * fPaircohPDF->Pcut > Pfactor) {
            continue;
         }
         ++fPaircohPDF->Npassed;
         break;
      }
   }

   // Generate new primaries for the pair
   double beamVelocity = GlueXPhotonBeamGenerator::getBeamVelocity();
   double steplength = pParticleChange->GetTrueStepLength();
   G4ThreeVector direction(track->GetMomentumDirection());
   G4ThreeVector x0(track->GetPosition());
   double t0 = track->GetGlobalTime();
   double uvtx = G4UniformRand();
   double lvtx = steplength + fPIL * log(1 - uvtx * (1 - exp(-steplength / fPIL)));
   x0 -= lvtx * direction;
   t0 -= lvtx / beamVelocity;
   G4PrimaryVertex vertex(x0, t0);
   G4ParticleDefinition *positron = G4Positron::Definition();
   G4ParticleDefinition *electron = G4Electron::Definition();
   G4ParticleDefinition *leptonplus = positron;
   G4ParticleDefinition *leptonminus = electron;
   if (fLeptonPairFamily == 1) {
      leptonplus = G4MuonPlus::Definition();
      leptonminus = G4MuonMinus::Definition();
   }
   G4ThreeVector psec1(p1.Mom()[1]*GeV, p1.Mom()[2]*GeV, p1.Mom()[3]*GeV);
   G4ThreeVector psec2(e2.Mom()[1]*GeV, e2.Mom()[2]*GeV, e2.Mom()[3]*GeV);
   G4ThreeVector psec3(e3.Mom()[1]*GeV, e3.Mom()[2]*GeV, e3.Mom()[3]*GeV);
   G4DynamicParticle *sec1 = new G4DynamicParticle(leptonplus, psec1);
   G4DynamicParticle *sec2 = new G4DynamicParticle(leptonminus, psec2);
   G4DynamicParticle *sec3 = new G4DynamicParticle(electron, psec3);
   G4TrackVector secondaries;
   secondaries.push_back(new G4Track(sec1, t0, x0));
   secondaries.push_back(new G4Track(sec2, t0, x0));
   if (e3.Mom().Length() > 0) {
      secondaries.push_back(new G4Track(sec3, t0, x0));
   }

   GlueXUserEventInformation *event_info;
   const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
   event_info = (GlueXUserEventInformation*)event->GetUserInformation();
   G4TrackVector::iterator iter;
   for (iter = secondaries.begin(); iter != secondaries.end(); ++iter) {
      GlueXUserTrackInformation *trackinfo = new GlueXUserTrackInformation();
      if (event_info) {
         trackinfo->SetGlueXTrackID(event_info->AssignNextGlueXTrackID());
      }
      (*iter)->SetUserInformation(trackinfo);
      if (fStopBeamAfterConverter == 0) {
         pParticleChange->AddSecondary(*iter);
      }
   }

   // append secondary vertex to MC record
   if (event_info) {
      int mech[2];
      char *cmech = (char*)mech;
      snprintf(cmech, 5, "%c%c%c%c", 'B', 'E', 'A', 'M');
      event_info->AddSecondaryVertex(secondaries, 1, mech[0]);
   }

#if VERBOSE_PAIRS_SPLITTING
   if ((fTripletPDF->Npassed + fPaircohPDF->Npassed) % 500 == 0) {
      G4cout << std::setprecision(5)
             << "triplet cross section is " << fTripletPDF->Pcut *
             fTripletPDF->Npassed / (fTripletPDF->Ntested + 1e-99) 
             << " +/- " << fTripletPDF->Pcut *
             sqrt(fTripletPDF->Npassed) / (fTripletPDF->Ntested + 1e-99) 
             << " barns, efficiency is " 
             << fTripletPDF->Npassed / (fTripletPDF->Ntested + 1e-99)
             << G4endl
             << "pair cross section is " << fPaircohPDF->Pcut *
             fPaircohPDF->Npassed / (fPaircohPDF->Ntested + 1e-99) 
             << " +/- " << fPaircohPDF->Pcut *
             sqrt(fPaircohPDF->Npassed) / (fPaircohPDF->Ntested + 1e-99) 
             << " barns, efficiency is " 
             << fPaircohPDF->Npassed / (fPaircohPDF->Ntested + 1e-99) 
             << G4endl
             << "counts are "
             << fTripletPDF->Npassed << " / " << fPaircohPDF->Npassed
             << " = "
             << fTripletPDF->Npassed / (fPaircohPDF->Npassed + 1e-99)
             << G4endl;
   }
#endif

#endif
}

void GlueXBeamConversionProcess::GenerateBetheHeitlerProcess(const G4Step &step)
{
   // Unlike the other GenerateXXX methods in this class, this method should
   // be invoked after tracking of an event is already underway. Its purpose
   // is to generate Bethe-Heitler pair conversion of the beam photon in the
   // primary GlueX target with the full polarization-dependent QED differential
   // cross section. A minimum invariant mass of the pair may be specified in
   // the member variable fBHpair_mass_min. It should be called from your
   // G4SteppingAction::UserSteppingAction method to force pair conversion
   // in special simulations dedicated to the study of the Bethe-Heitler
   // process. The incident gamma ray is stopped and the e+/e- pair vertex,
   // plus the recoil proton, is added to the output event record.

#ifndef USING_DIRACXX

   G4cerr << "GlueXBeamConversionProcess::GenerateBetheHeitlerProcess error:"
          << G4endl
          << "  You have enabled beam Bethe-Heitler conversion in the target,"
          << G4endl
          << "  but your have built HDGeant4 without support for the Dirac++"
		  << G4endl
          << "  library. Either rebuild with Dirac++ support or else turn off"
		  << G4endl
          << "  this process. Aborting this run..."
		  << G4endl;
   exit(1);

#else

   LDouble_t mLepton = mElectron;
   if (fLeptonPairFamily == 1)
      mLepton = mMuon;
   TPhoton gIn;
   TLepton eOut(mLepton);
   TLepton pOut(mLepton);
   TLepton nIn(mProton);
   TLepton nOut(mProton);
   LDouble_t diffXS(0);
   eOut.AllPol();
   pOut.AllPol();
   nOut.AllPol();
   const TThreeVectorReal zero3vector(0,0,0);
   nIn.SetPol(TThreeVectorReal(zero3vector));
   nIn.SetMom(TThreeVectorReal(zero3vector));
   const G4Track *track = step.GetTrack();
   LDouble_t kin = track->GetKineticEnergy()/GeV;

   // Decide which nuclear recoil state to generate
   enum recoilType { kNucleus, kProton, kNeutron };
   double mRecoil;
   if (fTargetA == 1) {
      mRecoil = (fTargetZ == 1)? mProton : mNeutron;
   }
   else {
      double uchan;
      unif01(1, &uchan);
      if (uchan < ELASTIC_NUCLEAR_FRACTION)
         mRecoil = nuclearMass_GeV();
      else if (uchan - ELASTIC_NUCLEAR_FRACTION < QUASIELASTIC_PROTON_FRACTION)
         mRecoil = mProton;
      else
         mRecoil = mNeutron;
   }
   nIn.SetMass(mRecoil);
   nOut.SetMass(mRecoil);

   // If we are below pair production threshold, do nothing
   if (kin < 2 * mLepton * (1 + mLepton / mRecoil))
      return;

   G4ThreeVector mom(track->GetMomentum());
   TThreeVectorReal mom0(mom[0]/GeV, mom[1]/GeV, mom[2]/GeV);
   gIn.SetMom(mom0);
   G4ThreeVector pol(track->GetPolarization());
   TThreeVectorReal pol0(pol[0], pol[1], pol[2]);
   if (pol0.Length() > 0)
      gIn.SetPlanePolarization(pol0, pol0.Length());
   else
      gIn.SetPol(pol0);

   // Define an angle and axis that rotates zhat into the direction
   // of the incident gamma, so that the generated kinematics is
   // defined with the incident gamma aligned with zhat, and then
   // rotated at the end into the final spatial direction.
   TThreeVectorReal rockaxis(mom0);
   rockaxis.Cross(TThreeVectorReal(0,0,1));
   double rockangle(0);
   if (mom0.Length() > 0) {
      rockangle = asin(rockaxis.Length() / mom0.Length());
      rockaxis /= rockaxis.Length();
   }

   // The beam photon energy is already generated by the cobrems
   // generator. Supply a normalized beam energy as the leading
   // elemnt of the u-vector in the call to sample(), and let it
   // return the remaining 5 random numbers needed to fix the
   // kinematics of the pair-production reaction.
   double E0_GeV = GlueXPrimaryGeneratorAction::getBeamEndpointEnergy()/GeV;
   double u[6];
   u[0] = kin / E0_GeV;
#if USE_ADAPTIVE_SAMPLER
   double weight = fAdaptiveSampler->sample(u);
#else
   unif01(5, u+1);
   double weight=1;
#endif

   // Generate uniform in E+, phi12, phiR
   LDouble_t Epos = kin * u[1];
   weight *= kin;
   LDouble_t phi12 = 2*M_PI * u[2];
   weight *= 2*M_PI;
   LDouble_t phiR = 2*M_PI * u[3];
   weight *= 2*M_PI;

   // Generate Mpair as 1 / (M [M^2 + Mcut^2])
   LDouble_t Mthresh = 2 * mLepton;
   LDouble_t Mcut = 10 * mLepton;  // 5 MeV for electron pairs
   LDouble_t um0 = 1 + sqr(Mcut / Mthresh);
   LDouble_t umax = 1;
   if (fBHpair_mass_min > Mthresh)
      umax = log(1 + sqr(Mcut / fBHpair_mass_min)) / log(um0);
   LDouble_t um = pow(um0, u[4] * umax);
   LDouble_t Mpair = Mcut / sqrt(um - 1 + 1e-99);
   weight *= Mpair * (sqr(Mcut) + sqr(Mpair)) * log(um0) / (2 * sqr(Mcut));
   weight *= umax;

   // Generate qR^2 with weight 1 / [qR^2 sqrt(qRcut^2 + qR^2)]
   LDouble_t qRmin = sqr(Mpair) /(2 * kin);
   LDouble_t qRcut = 2 * mLepton; // 1 MeV for electron pairs
   LDouble_t uq0 = qRmin / (qRcut + sqrt(sqr(qRcut) + sqr(qRmin)));
   LDouble_t uq = pow(uq0, u[5]);
   LDouble_t qR = 2 * qRcut * uq / (1 - sqr(uq));
   LDouble_t qR2 = qR * qR;
   weight *= qR2 * sqrt(1 + qR2 / sqr(qRcut)) * (-2 * log(uq0));

   // Include overall measure Jacobian factor
   weight *= Mpair / (2 * kin);

   try {
      if (Epos < mLepton) {
         throw std::runtime_error("positron energy less than its rest mass.");
      }

      // Solve for the initial state target kinematics
      LDouble_t Etarget(mRecoil);
      TFourVectorReal q0(Etarget, 0, 0, 0);
      if (Etarget < nuclearMass_GeV()) {
         double fermiMomentum = nuclearFermiMomentum_GeV() / sqrt(3);
         if (fermiMomentum > 0) {
            q0[1] = G4RandGauss::shoot(0, fermiMomentum);
            q0[2] = G4RandGauss::shoot(0, fermiMomentum);
            q0[3] = G4RandGauss::shoot(0, fermiMomentum);
         }
         Etarget -= nuclearBindingEnergy_GeV();
         q0[0] = Etarget;
      }
      LDouble_t q02 = q0.LengthSqr();

      // Solve for the recoil particle kinematics
      LDouble_t nu = sqrt(sqr(mRecoil) + qR2 + q02) - Etarget;
      LDouble_t costhetaR = (2 * kin * nu + qR2 + sqr(Mpair) - sqr(nu))
                            / (2 * qR * kin);
      for (int i=0; i < 999999999; ++i) {
         LDouble_t q0dotqR;
         if (fabs(costhetaR) <= 1) {
            LDouble_t sinthetaR = sqrt(1 - sqr(costhetaR));
            q0dotqR = q0[1] * qR * sinthetaR * cos(phiR) +
                      q0[2] * qR * sinthetaR * sin(phiR) +
                      q0[3] * qR * costhetaR;
         }
         else if (costhetaR > 1) {
            q0dotqR = -q0.Length() * qR;
         }
         else {
std::cerr << "no kinematic solution because costhetaR < -1" << std::endl;
            throw std::runtime_error("no kinematic solution because costhetaR < -1");
         }
         nu = sqrt(sqr(mRecoil) + qR2 + q02 + 2 * q0dotqR) - Etarget;
std::cerr << "nu=" << nu << ", qR=" << qR << ", q0=" << sqr(q02) << ", Etarget=" << Etarget << std::endl;
         LDouble_t last_costhetaR = costhetaR;
         costhetaR = (2 * kin * nu + qR2 + sqr(Mpair) - sqr(nu))
                     / (2 * qR * kin);
         if (costhetaR > 1) {
std::cerr << "no kinematic solution because costhetaR > 1" << std::endl;
            throw std::runtime_error("no kinematic solution because costhetaR > 1");
         }
         else if (fabs(costhetaR - last_costhetaR) < 1e-9) {
std::cerr << "costhetaR iterative solver converged at i=" << i << std::endl;
            break;
         }
         else if (i > 99) {
std::cerr << "no kinematic solution because costhetaR iterative solver failed to converge" << std::endl;
            throw std::runtime_error("no kinematic solution because costhetaR"
                                     " iterative solver failed to converge");
         }
      }
      TFourVectorReal q3;
      LDouble_t sinthetaR = sqrt(1 - sqr(costhetaR));
      q3[1] = q0[1] + qR * sinthetaR * cos(phiR);
      q3[2] = q0[2] + qR * sinthetaR * sin(phiR);
      q3[3] = q0[3] + qR * costhetaR;
      q3[0] = sqrt(q3.LengthSqr() + sqr(mRecoil));
      LDouble_t Erecoil(Etarget + nu);
      if (fabs(Erecoil - q3[0]) > 1e-6) {
std::cerr << "kinematics consistency check failed, cannot continue!!!" << std::endl;
std::cerr << "mRecoil=" << mRecoil << ", Etarget=" << Etarget << ", Erecoil=" << Erecoil << ", q3[0]=" << q3[0] << std::endl;
         throw std::runtime_error("no kinematic solution because Erecoil != q3[0]");
      }

      // Solve for the c.m. momentum of e+ in the pair 1,2 rest frame
      LDouble_t k12star2 = sqr(Mpair / 2) - sqr(mLepton);
      if (k12star2 < 0) {
         throw std::runtime_error("no kinematic solution because k12star2 < 0");
      }
      LDouble_t k12star = sqrt(k12star2);
      LDouble_t E12 = kin + Etarget - Erecoil;
      if (E12 < Mpair) {
std::cerr << "throw up #1 with mRecoil=" << mRecoil << std::endl;
         throw std::runtime_error("no kinematic solution because E12 < Mpair");
      }
      LDouble_t q12mag = sqrt(sqr(E12) - sqr(Mpair));
      LDouble_t costhetastar = (Epos - E12 / 2) * Mpair / (k12star * q12mag);
      if (Epos > E12 - mLepton) {
std::cerr << "throw up #2 with mRecoil=" << mRecoil << std::endl;
         throw std::runtime_error("no kinematic solution because Epos > E12 - mLepton");
      }
      else if (fabs(costhetastar) > 1) {
std::cerr << "throw up #3 with mRecoil=" << mRecoil << std::endl;
         throw std::runtime_error("no kinematic solution because |costhetastar| > 1");
      }

      // Boost the pair momenta into the lab
      TThreeVectorReal beta((q3[1] - q0[1]) / E12,
                            (q3[2] - q0[2]) / E12,
                            (q3[3] - q0[3] - kin) / E12);
      TLorentzBoost toLab(beta);
      toLab.SetGamma(E12 / Mpair);
      LDouble_t sinthetastar = sqrt(1 - sqr(costhetastar));
      TThreeVectorReal k12(k12star * sinthetastar * cos(phi12),
                           k12star * sinthetastar * sin(phi12),
                           k12star * costhetastar);
      TFourVectorReal q1(Mpair / 2, k12);
      TFourVectorReal q2(Mpair / 2, -k12);
      q1.Boost(toLab);
      q2.Boost(toLab);

      // Rotate all momenta to the lab frame
      if (rockangle != 0) {
         q0.Rotate(rockaxis, rockangle);
         q1.Rotate(rockaxis, rockangle);
         q2.Rotate(rockaxis, rockangle);
         q3.Rotate(rockaxis, rockangle);
      }
      nIn.SetMom(q0);
      pOut.SetMom(q1);
      eOut.SetMom(q2);
      nOut.SetMom(q3);
 
      // Check 4-momentum conservation
      TFourVectorReal pIn(gIn.Mom() + nIn.Mom());
      TFourVectorReal pFi(pOut.Mom() + eOut.Mom() + nOut.Mom());
      TFourVectorReal::SetResolution(1e-10);
      if (pIn != pFi) {
         std::streamsize defprec = std::cout.precision();
         std::cout << std::setprecision(12);
         std::cout << "Warning in GenerateBetheHeitlerConversion - "
                   << "momentum conservation violated." << std::endl
                   << "   pIn = "; 
         pIn.Print();
         std::cout << "       = ";
         gIn.Mom().Print();
         std::cout << "       + ";
         nIn.Mom().Print();
         std::cout << "   pFi = ";
         pFi.Print();
         std::cout << "       = ";
         pOut.Mom().Print();
         std::cout << "       + ";
         eOut.Mom().Print();
         std::cout << "       + ";
         nOut.Mom().Print();
         std::cout << "   pIn - pFi = ";
         (pIn-pFi).Print();
         std::cout << std::setprecision(defprec);
      }

      // Compute the polarized differential cross section (barnes/GeV^4)
      // returned as d(sigma)/(dE+ dphi+ d^3qR)
      LDouble_t t = sqr(Erecoil - Etarget) - qR2;
      if (fTargetA > 1 && mRecoil == nuclearMass_GeV()) {
         diffXS = TCrossSection::PairProduction(gIn, eOut, pOut);
         fPairsGeneration->SetConverterZ(fTargetZ);
         diffXS *= sqr(1 - fPairsGeneration->FFatomic(qR));
         diffXS *= sqr(fTargetZ * nuclearFormFactor(-t));
         weight /= ELASTIC_NUCLEAR_FRACTION;
      }
      else if (mRecoil == mProton) {
         LDouble_t F1_timelike;
         LDouble_t F2_timelike;
         LDouble_t F1_spacelike;
         LDouble_t F2_spacelike;
         if (fFormFactorChoice == 0) {
            F1_timelike = 1;
            F2_timelike = 0;
            F1_spacelike = nucleonFormFactor(-t, kF1p);
            F2_spacelike = nucleonFormFactor(-t, kF2p);
         }
         else if (fFormFactorChoice == 1) {
            F1_timelike = 1;
            F2_timelike = 0;
            F1_spacelike = nucleonFormFactor(-t, kF1p);
            F2_spacelike = nucleonFormFactor(-t, kF2p);
         }
         else if (fFormFactorChoice == 3) {
            F1_timelike = 0;
            F2_timelike = 0;
            F1_spacelike = nucleonFormFactor(-t, kF1p);
            F2_spacelike = nucleonFormFactor(-t, kF2p);
         }
         else if (fFormFactorChoice == 4) {
            F1_timelike = 0;
            F2_timelike = 0;
            F1_spacelike = nucleonFormFactor(-t, kF1p);
            F2_spacelike = 0;
         }
         else if (fFormFactorChoice == 5) {
            F1_timelike = 1;
            F2_timelike = 1;
            F1_spacelike = 0;
            F2_spacelike = 0;
         }
         else if (fFormFactorChoice == 6) {
            F1_timelike = 1;
            F2_timelike = 0;
            F1_spacelike = 0;
            F2_spacelike = 0;
         }
         else if (fFormFactorChoice == 7) {
            F1_timelike = 0;
            F2_timelike = 1;
            F1_spacelike = 0;
            F2_spacelike = 0;
         }
         else {
            G4cerr << "GlueXBeamConversionProcess::GenerateBetheHeitlerProcess"
                   << " error: Unsupported nucleon form factor option" 
                   << ", cannot continue!" << G4endl;
		    exit(89);
         }
         diffXS = TCrossSection::BetheHeitlerNucleon(gIn, nIn,
                                                     pOut, eOut, nOut,
                                                     F1_spacelike,
                                                     F2_spacelike,
                                                     F1_timelike,
                                                     F2_timelike);
         if (fTargetZ == 1) {
            fPairsGeneration->SetConverterZ(fTargetZ);
            diffXS *= sqr(1 - fPairsGeneration->FFatomic(qR));
         }
         else {
            diffXS *= fTargetZ * exp(-fTargetA * sqr(nuclearFormFactor(-t)));
            weight /= QUASIELASTIC_PROTON_FRACTION;
         }
      }
      else if (mRecoil == mNeutron) {
         LDouble_t F1_timelike = 1;
         LDouble_t F2_timelike = 0;
         LDouble_t F1_spacelike = nucleonFormFactor(-t, kF1n);
         LDouble_t F2_spacelike = nucleonFormFactor(-t, kF2n);
         diffXS = TCrossSection::BetheHeitlerNucleon(gIn, nIn,
                                                     pOut, eOut, nOut,
                                                     F1_spacelike,
                                                     F2_spacelike,
                                                     F1_timelike,
                                                     F2_timelike);
         if (fTargetA > 1) {
            diffXS *= (fTargetA - fTargetZ) * exp(-fTargetA * sqr(nuclearFormFactor(-t)));
            weight /= QUASIELASTIC_NEUTRON_FRACTION;
         }
      }
      else {
         G4cerr << "GlueXBeamConversionProcess::GenerateBetheHeitlerProcess error:"
                << "  Unknown recoil particle of mass " << mRecoil
                << ", cannot continue!" << G4endl;
		 exit(82);
      }
#if USE_ADAPTIVE_SAMPLER
      fAdaptiveSampler->feedback(u, weight * diffXS);
#endif
   }
   catch (const std::exception &e) {

      // These events have no cross section, but do not discard
      // them because they are needed to get the right MC integral.
 
std::cerr << "catching throw up with mRecoil=" << mRecoil << std::endl;
      nOut.SetMom(TThreeVectorReal(0,0,1e-12));
      eOut.SetMom(TThreeVectorReal(0,0,1e-12));
      pOut.SetMom(TThreeVectorReal(0,0,1e-12));
#if USE_ADAPTIVE_SAMPLER
      fAdaptiveSampler->feedback(u, 0);
#endif
      diffXS = 0;
   }

   // Generate the new vertex
   double beamVelocity = GlueXPhotonBeamGenerator::getBeamVelocity();
   double steplength = pParticleChange->GetTrueStepLength();
   G4ThreeVector direction(track->GetMomentumDirection());
   G4ThreeVector x0(track->GetPosition());
   double t0 = track->GetGlobalTime();
   double uvtx = G4UniformRand();
   x0 -= uvtx * steplength * direction;
   t0 -= uvtx * steplength / beamVelocity;
   G4PrimaryVertex vertex(x0, t0);
   G4ParticleDefinition *leptonplus = G4Positron::Definition();
   G4ParticleDefinition *leptonminus = G4Electron::Definition();
   if (fLeptonPairFamily == 1) {
      leptonplus = G4MuonPlus::Definition();
      leptonminus = G4MuonMinus::Definition();
   }
   G4ThreeVector psec1(pOut.Mom()[1]*GeV, pOut.Mom()[2]*GeV, pOut.Mom()[3]*GeV);
   G4ThreeVector psec2(eOut.Mom()[1]*GeV, eOut.Mom()[2]*GeV, eOut.Mom()[3]*GeV);
   G4ThreeVector psec3(nOut.Mom()[1]*GeV, nOut.Mom()[2]*GeV, nOut.Mom()[3]*GeV);
   G4DynamicParticle *sec1 = new G4DynamicParticle(leptonplus, psec1);
   G4DynamicParticle *sec2 = new G4DynamicParticle(leptonminus, psec2);
   G4DynamicParticle *sec3;
   if (mRecoil == mProton) {
      G4ParticleDefinition *proton = G4Proton::Definition();
      sec3 = new G4DynamicParticle(proton, psec3);
   }
   else if (mRecoil == mNeutron) {
      G4ParticleDefinition *neutron = G4Neutron::Definition();
      sec3 = new G4DynamicParticle(neutron, psec3);
   }
   else {
      G4ParticleDefinition *ion = G4ParticleTable::GetParticleTable()
                           ->GetIonTable()->GetIon(fTargetZ, fTargetA);
      sec3 = new G4DynamicParticle(ion, psec3);
   }
   G4TrackVector secondaries;
   secondaries.push_back(new G4Track(sec1, t0, x0));
   secondaries.push_back(new G4Track(sec2, t0, x0));
   if (nOut.Mom().Length() > 0) {
      secondaries.push_back(new G4Track(sec3, t0, x0));
   }

   GlueXUserEventInformation *event_info;
   const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
   event_info = (GlueXUserEventInformation*)event->GetUserInformation();
   G4TrackVector::iterator iter;
   for (iter = secondaries.begin(); iter != secondaries.end(); ++iter) {
      GlueXUserTrackInformation *trackinfo = new GlueXUserTrackInformation();
      if (event_info) {
         trackinfo->SetGlueXTrackID(event_info->AssignNextGlueXTrackID());
      }
      (*iter)->SetUserInformation(trackinfo);
      if (fStopBeamAfterTarget == 0) {
         pParticleChange->AddSecondary(*iter);
      }
   }

   // append secondary vertex to MC record
   if (event_info) {
      int mech[2];
      char *cmech = (char*)mech;
      snprintf(cmech, 5, "%c%c%c%c", 'C', 'O', 'N', 'V');
      event_info->AddSecondaryVertex(secondaries, 1, mech[0]);
      hddm_s::ReactionList rea = event_info->getOutputRecord()->getReactions();
      if (rea.size() > 0) {
         rea(0).setType(221); // Bethe-Heitler process
         rea(0).setWeight(weight * diffXS);
      }
   }

#if USE_ADAPTIVE_SAMPLER
   static int passes = 0;
   ++passes; // thread safety not important here
   int tens = 1000;
   while (passes >= tens*10)
      tens *= 10;
   if (passes % tens == 0) {
      std::cout << "AdaptiveSampler reports efficiency " 
                << fAdaptiveSampler->getEfficiency()
                << std::endl;
      std::stringstream astatefile;
      astatefile << "BHgen_thread_" << G4Threading::G4GetThreadId() << ".astate";
      fAdaptiveSampler->saveState(astatefile.str());
   }
#endif

#endif
}

void GlueXBeamConversionProcess::GenerateTripletProcess(const G4Step &step)
{
   // This method is completemary to GenerateBetheHeitlerProcess, for the case
   // of pair production off a target electron rather than coherently from the
   // atom as a whole. It generates the triplet process in the GlueX liquid
   // hydrogen target with the full polarization-dependent QED differential
   // cross section. A minimum invariant mass of the pair may be specified in
   // the member variable fBHpair_mass_min. It should be called from your
   // G4SteppingAction::UserSteppingAction method to force pair conversion
   // in special simulations dedicated to the study of the triplet proces.
   // The incident gamma ray is stopped and the e+/e- pair vertex plus the
   // recoil electron is added to the output event record.

#ifndef USING_DIRACXX

   G4cerr << "GlueXBeamConversionProcess::GenerateTripletProcess error:"
          << G4endl
          << "  You have enabled beam triplet conversion in the target,"
          << G4endl
          << "  but your have built HDGeant4 without support for the Dirac++"
		  << G4endl
          << "  library. Either rebuild with Dirac++ support or else turn off"
		  << G4endl
          << "  this process. Aborting this run..."
		  << G4endl;
   exit(1);

#else

   LDouble_t mLepton = mElectron;
   if (fLeptonPairFamily == 1)
      mLepton = mMuon;
   TPhoton gIn;
   TLepton eOut(mLepton);
   TLepton pOut(mLepton);
   TLepton eIn(mElectron);
   TLepton eOut3(mElectron);
   LDouble_t diffXS(0);
   eOut.AllPol();
   pOut.AllPol();
   eOut3.AllPol();
   const TThreeVectorReal zero3vector(0,0,0);
   eIn.SetPol(TThreeVectorReal(zero3vector));
   eIn.SetMom(TThreeVectorReal(zero3vector));
   const G4Track *track = step.GetTrack();
   LDouble_t kin = track->GetKineticEnergy()/GeV;

   // If we are below pair production threshold, do nothing
   if (kin < 2 * mLepton * (1 + mLepton / mElectron))
      return;

   G4ThreeVector mom(track->GetMomentum());
   TThreeVectorReal mom0(mom[0]/GeV, mom[1]/GeV, mom[2]/GeV);
   gIn.SetMom(mom0);
   G4ThreeVector pol(track->GetPolarization());
   TThreeVectorReal pol0(pol[0], pol[1], pol[2]);
   if (pol0.Length() > 0)
      gIn.SetPlanePolarization(pol0, pol0.Length());
   else
      gIn.SetPol(pol0);

   // Define an angle and axis that rotates zhat into the direction
   // of the incident gamma, so that the generated kinematics is
   // defined with the incident gamma aligned with zhat, and then
   // rotated at the end into the final spatial direction.
   TThreeVectorReal rockaxis(mom0);
   rockaxis.Cross(TThreeVectorReal(0,0,1));
   double rockangle(0);
   if (mom0.Length() > 0) {
      rockangle = asin(rockaxis.Length() / mom0.Length());
      rockaxis /= rockaxis.Length();
   }

   // The beam photon energy is already generated by the cobrems
   // generator. Supply a normalized beam energy as the leading
   // elemnt of the u-vector in the call to sample(), and let it
   // return the remaining 5 random numbers needed to fix the
   // kinematics of the pair-production reaction.
   double E0_GeV = GlueXPrimaryGeneratorAction::getBeamEndpointEnergy()/GeV;
   double u[6];
   u[0] = kin / E0_GeV;
#if USE_ADAPTIVE_SAMPLER
   double weight = fAdaptiveSampler->sample(u);
#else
   unif01(5, u+1);
   double weight=1;
#endif

   // Generate uniform in E+, phi12, phiR
   LDouble_t Epos = kin * u[1];
   weight *= kin;
   LDouble_t phi12 = 2*M_PI * u[2];
   weight *= 2*M_PI;
   LDouble_t phiR = 2*M_PI * u[3];
   weight *= 2*M_PI;

   // Generate Mpair as 1 / (M [M^2 + Mcut^2])
   LDouble_t Mthresh = 2 * mLepton;
   LDouble_t Mcut = 10 * mLepton;  // 5 MeV for electrons
   LDouble_t um0 = 1 + sqr(Mcut / Mthresh);
   LDouble_t umax = 1;
   if (fBHpair_mass_min > Mthresh)
      umax = log(1 + sqr(Mcut / fBHpair_mass_min)) / log(um0);
   LDouble_t um = pow(um0, u[4] * umax);
   LDouble_t Mpair = Mcut / sqrt(um - 1 + 1e-99);
   weight *= Mpair * (sqr(Mcut) + sqr(Mpair)) * log(um0) / (2 * sqr(Mcut));
   weight *= umax;

   // Generate qR^2 with weight 1 / [qR^2 sqrt(qRcut^2 + qR^2)]
   LDouble_t qRmin = sqr(Mpair) /(2 * kin);
   LDouble_t qRcut = 2 * mLepton; // 1 MeV for electrons
   LDouble_t uq0 = qRmin / (qRcut + sqrt(sqr(qRcut) + sqr(qRmin)));
   LDouble_t uq = pow(uq0, u[5]);
   LDouble_t qR = 2 * qRcut * uq / (1 - sqr(uq));
   LDouble_t qR2 = qR * qR;
   weight *= qR2 * sqrt(1 + qR2 / sqr(qRcut)) * (-2 * log(uq0));

   // Include overall measure Jacobian factor
   weight *= Mpair / (2 * kin);

   try {
      if (Epos < mLepton) {
         throw std::runtime_error("positive lepton energy less than its rest mass.");
      }

      // Solve for the c.m. momentum of e+ in the pair 1,2 rest frame
      LDouble_t k12star2 = sqr(Mpair / 2) - sqr(mLepton);
      if (k12star2 < 0) {
         throw std::runtime_error("no kinematic solution because k12star2 < 0");
      }
      LDouble_t k12star = sqrt(k12star2);
      LDouble_t Erec = sqrt(qR2 + sqr(mElectron));
      LDouble_t E12 = kin + mElectron - Erec;
      if (E12 < Mpair) {
         throw std::runtime_error("no kinematic solution because E12 < Mpair");
      }
      LDouble_t q12mag = sqrt(sqr(E12) - sqr(Mpair));
      LDouble_t costhetastar = (Epos - E12 / 2) * Mpair / (k12star * q12mag);
      if (Epos > E12 - mLepton) {
         throw std::runtime_error("no kinematic solution because Epos > E12 - mLepton");
      }
      else if (fabs(costhetastar) > 1) {
         throw std::runtime_error("no kinematic solution because |costhetastar| > 1");
      }

      // Solve for the recoil electron kinematics
      LDouble_t costhetaR = (sqr(Mpair) / 2 + (kin + mElectron) *
                            (Erec - mElectron)) / (kin * qR);
      if (fabs(costhetaR) > 1) {
         throw std::runtime_error("no kinematic solution because |costhetaR| > 1");
      }
      LDouble_t sinthetaR = sqrt(1 - sqr(costhetaR));
      TFourVectorReal q3(Erec, qR * sinthetaR * cos(phiR),
                               qR * sinthetaR * sin(phiR),
                               qR * costhetaR);

      // Boost the pair momenta into the lab
      TLorentzBoost toLab(q3[1] / E12, q3[2] / E12, (q3[3] - kin) / E12);
      LDouble_t sinthetastar = sqrt(1 - sqr(costhetastar));
      TThreeVectorReal k12(k12star * sinthetastar * cos(phi12),
                           k12star * sinthetastar * sin(phi12),
                           k12star * costhetastar);
      TFourVectorReal q1(Mpair / 2, k12);
      TFourVectorReal q2(Mpair / 2, -k12);
      q1.Boost(toLab);
      q2.Boost(toLab);

      // Rotate final state momenta to the lab frame
      if (rockangle != 0) {
         q1.Rotate(rockaxis, rockangle);
         q2.Rotate(rockaxis, rockangle);
         q3.Rotate(rockaxis, rockangle);
      }
      pOut.SetMom(q1);
      eOut.SetMom(q2);
      eOut3.SetMom(q3);

      // Check 4-momentum conservation
      TFourVectorReal pIn(gIn.Mom() + eIn.Mom());
      TFourVectorReal pFi(pOut.Mom() + eOut.Mom() + eOut3.Mom());
      TFourVectorReal::SetResolution(1e-10);
      if (pIn != pFi) {
         std::cout << "Warning in GenerateTripletConversion - "
                   << "momentum conservation violated." << std::endl
                   << "   pIn = "; 
         pIn.Print();
         std::cout << "       = ";
         gIn.Mom().Print();
         std::cout << "       + ";
         eIn.Mom().Print();
         std::cout << "   pFi = ";
         pFi.Print();
         std::cout << "       = ";
         pOut.Mom().Print();
         std::cout << "       + ";
         eOut.Mom().Print();
         std::cout << "       + ";
         eOut3.Mom().Print();
         std::cout << "   pIn - pFi = ";
         (pIn-pFi).Print();
      }

      // Compute the polarized differential cross section (barnes/GeV^4)
      // returned as d(sigma)/(dE+ dphi+ d^3qR)
      diffXS = TCrossSection::TripletProduction(gIn, eIn, pOut, eOut, eOut3);
      fPairsGeneration->SetConverterZ(fTargetZ);
      diffXS *= fTargetZ * exp(-fTargetZ * sqr(fPairsGeneration->FFatomic(qR)));
#if USE_ADAPTIVE_SAMPLER
      fAdaptiveSampler->feedback(u, weight * diffXS);
#endif
   }
   catch (const std::exception &e) {

      // These events have no cross section, but do not discard
      // them because they are needed to get the right MC integral.
 
      eOut3.SetMom(TThreeVectorReal(0,0,1e-12));
      eOut.SetMom(TThreeVectorReal(0,0,1e-12));
      pOut.SetMom(TThreeVectorReal(0,0,1e-12));
#if USE_ADAPTIVE_SAMPLER
      fAdaptiveSampler->feedback(u, 0);
#endif
      diffXS = 0;
   }

   // Generate the new vertex
   double beamVelocity = GlueXPhotonBeamGenerator::getBeamVelocity();
   double steplength = pParticleChange->GetTrueStepLength();
   G4ThreeVector direction(track->GetMomentumDirection());
   G4ThreeVector x0(track->GetPosition());
   double t0 = track->GetGlobalTime();
   double uvtx = G4UniformRand();
   x0 -= uvtx * steplength * direction;
   t0 -= uvtx * steplength / beamVelocity;
   G4PrimaryVertex vertex(x0, t0);
   G4ParticleDefinition *positron = G4Positron::Definition();
   G4ParticleDefinition *electron = G4Electron::Definition();
   G4ParticleDefinition *leptonplus = positron;
   G4ParticleDefinition *leptonminus = electron;
   if (fLeptonPairFamily == 1) {
      leptonplus = G4MuonPlus::Definition();
      leptonminus = G4MuonMinus::Definition();
   }
   G4ThreeVector psec1(pOut.Mom()[1]*GeV, pOut.Mom()[2]*GeV, pOut.Mom()[3]*GeV);
   G4ThreeVector psec2(eOut.Mom()[1]*GeV, eOut.Mom()[2]*GeV, eOut.Mom()[3]*GeV);
   G4ThreeVector psec3(eOut3.Mom()[1]*GeV, eOut3.Mom()[2]*GeV, eOut3.Mom()[3]*GeV);
   G4DynamicParticle *sec1 = new G4DynamicParticle(leptonplus, psec1);
   G4DynamicParticle *sec2 = new G4DynamicParticle(leptonminus, psec2);
   G4DynamicParticle *sec3 = new G4DynamicParticle(electron, psec3);
   G4TrackVector secondaries;
   secondaries.push_back(new G4Track(sec1, t0, x0));
   secondaries.push_back(new G4Track(sec2, t0, x0));
   if (eOut3.Mom().Length() > 0) {
      secondaries.push_back(new G4Track(sec3, t0, x0));
   }

   GlueXUserEventInformation *event_info;
   const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
   event_info = (GlueXUserEventInformation*)event->GetUserInformation();
   G4TrackVector::iterator iter;
   for (iter = secondaries.begin(); iter != secondaries.end(); ++iter) {
      GlueXUserTrackInformation *trackinfo = new GlueXUserTrackInformation();
      if (event_info) {
         trackinfo->SetGlueXTrackID(event_info->AssignNextGlueXTrackID());
      }
      (*iter)->SetUserInformation(trackinfo);
      if (fStopBeamAfterTarget == 0) {
         pParticleChange->AddSecondary(*iter);
      }
   }

   // append secondary vertex to MC record
   if (event_info) {
      int mech[2];
      char *cmech = (char*)mech;
      snprintf(cmech, 5, "%c%c%c%c", 'C', 'O', 'N', 'V');
      event_info->AddSecondaryVertex(secondaries, 1, mech[0]);
      hddm_s::ReactionList rea = event_info->getOutputRecord()->getReactions();
      if (rea.size() > 0) {
         rea(0).setType(231); // triplet process
         rea(0).setWeight(weight * diffXS);
      }
   }

#if USE_ADAPTIVE_SAMPLER
   static int passes = 0;
   ++passes; // thread safety not important here
   int tens = 1000;
   while (passes >= tens*10)
      tens *= 10;
   if (passes % tens == 0) {
      std::cout << "AdaptiveSampler reports efficiency " 
                << fAdaptiveSampler->getEfficiency()
                << std::endl;
      std::stringstream astatefile;
      astatefile << "BHgen_thread_" << G4Threading::G4GetThreadId() << ".astate";
      fAdaptiveSampler->saveState(astatefile.str());
   }
#endif

#endif
}

void GlueXBeamConversionProcess::prepareAdaptiveSampler()
{
   fAdaptiveSampler = new AdaptiveSampler(6, &unif01, 1);
   fAdaptiveSampler->restoreState("BHgen.astate");
   fAdaptiveSampler->setVerbosity(2);
   fAdaptiveSampler->reset_stats();

   G4AutoLock barrier(&fMutex);
   fAdaptiveSamplerRegistry.push_back(fAdaptiveSampler);
}

double GlueXBeamConversionProcess::nucleonFormFactor(double Q2_GeV,
                       GlueXBeamConversionProcess::FormFactorType t)
{
   // The following fits were taken from J.J. Kelly, Simple Parameterization
   // of Nucleon Form Factors, PHYSICAL REVIEW C 70, 068202 (2004)

   const double proton_magnetic_moment = 2.793;
   const double neutron_magnetic_moment = -1.913;

   const std::map<GlueXBeamConversionProcess::FormFactorType, std::vector<double> >
      a = {{kGEp, { 0.24 }},
           {kGMp, { 0.12 }},
           {kGMn, { 2.33 }},
           {kGEn, { 1.70 }}
          };
   const std::map<GlueXBeamConversionProcess::FormFactorType, std::vector<double> >
      b = {{kGEp, { 10.98, 12.82, 21.97 }},
           {kGMp, { 10.97, 18.86,  6.55 }},
           {kGMn, { 14.72, 24.20, 84.10 }},
           {kGEn, { 3.30 }}
          };

   double tau = Q2_GeV / sqr(2 * mProton);
   if (t == kGEp || t == kGMp || t == kGMn) {
      double numer = 1;
      double denom = 1;
      for (unsigned n=0; n < a.at(t).size(); ++n)
         numer += a.at(t)[n] * pow(tau, n+1);
      for (unsigned n=0; n < b.at(t).size(); ++n)
         denom += b.at(t)[n] * pow(tau, n+1);
      if (t == kGMp)
         numer *= proton_magnetic_moment;
      else if (t == kGMn)
         numer *= neutron_magnetic_moment;
      return numer / denom;
   }
   else if (t == kGEn) {
      double GDipole = 1 / sqr(1 + Q2_GeV/0.71);
      return GDipole * a.at(t)[0] * tau / (1 + b.at(t)[0] * tau);
   }
   else if (t == kF1p) {
      return (nucleonFormFactor(Q2_GeV, kGEp) + 
              nucleonFormFactor(Q2_GeV, kGMp) * tau) / (1 + tau);
   }
   else if (t == kF2p) {
      return (nucleonFormFactor(Q2_GeV, kGMp) -
              nucleonFormFactor(Q2_GeV, kGEp)) / (1 + tau);
   }
   else if (t == kF1n) {
      return (nucleonFormFactor(Q2_GeV, kGEn) + 
              nucleonFormFactor(Q2_GeV, kGMn) * tau) / (1 + tau);
   }
   else if (t == kF2n) {
      return (nucleonFormFactor(Q2_GeV, kGMn) -
              nucleonFormFactor(Q2_GeV, kGEn)) / (1 + tau);
   }

   std::cerr << "GlueXBeamConversionProcess::nucleonFormFactor error - "
             << "unknown form factor requested, returning zero!" << std::endl;
   return 0;
}

double GlueXBeamConversionProcess::nuclearFormFactor(double Q2_GeV)
{
   // This fit to the parameterization of the nuclear charge form factor 
   // was provided by Rory Miskimen, with attribution to Bernard Frois.
   // It was translated from the original Fortran into C++ by
   // R.T. Jones, July 17, 2021.

   const std::map<double, double>
      A = {{0.1, 0.003845},
           {0.7, 0.009724},
           {1.6, 0.033093},
           {2.1, 0.000120},
           {2.7, 0.083107},
           {3.5, 0.080869},
           {4.2, 0.139957},
           {5.1, 0.260892},
           {6.0, 0.336013},
           {6.6, 0.033637},
           {7.6, 0.018729},
           {8.7, 0.000020}
         };
   const double gamma(1.388);
   const double hbarc(0.197);

   double q = sqrt(Q2_GeV);
   double FF = 0;
   if (fTargetZ == 82) {
      for (auto rA : A) {
          double r = rA.first;
          FF +=  A.at(r) * (sqr(gamma) * cos(q * r/hbarc) +
                         2 * r * hbarc/q * sin(q * r/hbarc))
                        / (sqr(gamma) + 2 * sqr(r)) *
                        exp(-sqr(q) / 4 * sqr(gamma/hbarc));
      }
   }
   else if (fTargetZ > 1) {

      // Take the Fourier transform of the Saxon-Woods charge distribution
      // with the mean charge radius parameterized as a function of Z.
      // This parameterization is taken from the paper by Maximon and Schrack,
      // "The Form Factor of the Fermi Model Spatial Distribution",
      // JOURNAL OF RESEARCH of the Notional Bureau of Standards - 
      // B. Mathematics and Mathematical Physics Vol. 70B, No. 1,
      // January-March 1966. The approximate parameterization for the
      // Saxon-Woods parameters a and c are taken from "Electron Scattering:
      // Form Factors and Nuclear Shapes, Kyle Foster, Nov. 16, 2011.
   
      double adent = 0.55; // fm, corresponds to t = 2.4 fm
      if (fTargetA == 4)   // special case for 4He nucleus
         adent = 0.29; 
      double adent2 = adent * adent;
      double adent3 = adent2 * adent;
      double cdent = 1.07 * pow(fTargetA, 1/3.);  // fm
      double q__fm = q / hbarc; // fm^-1
      double FF0 = 4 * M_PI * cdent / 3 * (sqr(M_PI * adent) + sqr(cdent));
      FF = 4 * sqr(M_PI) * adent3
           / sqr(q__fm * adent * sinh(M_PI * q__fm * adent))
           * (M_PI * q__fm * adent * cosh(M_PI * q__fm * adent)
                               * sin(q__fm * cdent)
                   - q__fm * cdent * sinh(M_PI * q__fm * adent)
                                  * cos(q__fm * cdent));
      int sgn = 1;
      for (int n=1; n < 10; ++n) {
         FF0 += 8. * M_PI * adent3 * sgn *
                 exp(-n * cdent / adent) / pow(n, 3);
         FF += 8. * M_PI * adent3 * sgn *
                  n * exp(-n * cdent / adent) /
                  sqr(n*n + sqr(q__fm * adent));
         sgn *= -1;
      }
      FF /= FF0; // normalize to FF(0) = 1
   }
   else {
      std::cerr << "GlueXBeamConversionProcess::nuclearFormFactor error - "
                << "cannot currently handle target with Z=" << fTargetZ
                << " A=" << fTargetA << ", cannot continue!"
                << std::endl;
      exit(92);
   }
   return FF;
}

double GlueXBeamConversionProcess::nuclearMass_GeV()
{
   if (fTargetA == 1 and fTargetZ == 0)
      return mNeutron;
   else if (fTargetA == 1 and fTargetZ == 1)
      return mProton;
   return G4NistManager::Instance()->GetAtomicMassAmu(fTargetZ) * AMU_GEV;
}

double GlueXBeamConversionProcess::nuclearFermiMomentum_GeV()
{
   // use a simple Fermi gas model of the nucleus to find
   // the RMS momentum of bound nucleons in the nucleus.

   double fermiEnergy(0);
   if (fTargetA == 2)
      fermiEnergy = 9e-3;
   else if (fTargetA == 3)
      fermiEnergy = 20e-3;
   else if (fTargetA > 3)
      fermiEnergy = 33e-3;
   return sqrt(2 * AMU_GEV * fermiEnergy);
}

double GlueXBeamConversionProcess::nuclearBindingEnergy_GeV()
{
   // use the Bethe-Weiztacher formula

   double a=14.0;
   double b=13.0;
   double c=0.585;
   double d=19.3;
   double e=33.0;
   double B_MeV = a - b / pow(fTargetA, 1/3.) 
                    - c * pow(fTargetZ, 2) / pow(fTargetA, 4/3.)
                    - d * pow((fTargetA - 2 * fTargetZ) / fTargetA, 2);
   if ((int(fTargetA + 0.5) % 2) == 0) {
      if ((int(fTargetZ + 0.5) % 2) == 0)
         B_MeV += e / pow(fTargetA, 7/4.);
      else
         B_MeV -= e / pow(fTargetA, 7/4.);
   }
   return B_MeV * 1e-3;
}
