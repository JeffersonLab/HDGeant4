//
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
#include "G4KaonZeroLong.hh"
#include "G4KaonZeroShort.hh"
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

#define CHECK_KINEMATICS 1

// Phi photoproduction total cross section as function of photon energy,
// digitized from Fig. 3 in Wang et al, arXiv:2208.10289 (22 Aug 2022).

double gammaPhiXS_dE(0.1*GeV);
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
double gammaPhiXS_scale_factor(1e3);
double gammaPhiXS_tslope(6.2/(GeV*GeV));
double gammaPhiXS_Emin(1.9*GeV);
double gammaPhiXS_Emax(22.0*GeV);
double neutralKaonPhi_bratio(0.342);
double massTargetNucleon(0.938*GeV);
double massPhiMeson(1.020*GeV);
double massK0Meson(0.4976*GeV);
double NAvagadro(6.023e23);
double cm2_ub(1e-30);

// conversion target material properties
double beryllium_A(9);
double beryllium_density_gcm3(1.85);
double tungsten_A(0.9 * 183.8 + 0.1 * 63.5);
double tungsten_density_gcm3(1 / ((0.9 / 19.38) + (0.1 / 8.96)));

G4Mutex GlueXKlongConversionProcess::fMutex = G4MUTEX_INITIALIZER;
int GlueXKlongConversionProcess::fConfigured = 0;

GlueXKlongConversionProcess::GlueXKlongConversionProcess(
                                       const G4String &name, 
                                       G4ProcessType aType)
 : G4VEmProcess(name, aType),
   isInitialised(false)
{
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
         //fStopBeamAfterConverter = 1;
      }
   }

   fConfigured = 1;

   if (verboseLevel > 0) {
       G4cout << GetProcessName() << " is created " << G4endl;
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
      G4double emin = gammaPhiXS_Emin;
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
  return gammaPhiXS_Emin;
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
   // Use of this process precludes use of the tagger
   GlueXUserEventInformation *event_info;
   const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
   event_info = (GlueXUserEventInformation*)event->GetUserInformation();
   if (event_info) {
      hddm_s::ReactionList rea = event_info->getOutputRecord()->getReactions();
      if (rea.size() > 0) {
         hddm_s::VertexList ver = rea(0).getVertices();
         if (ver.size() > 0) {
            hddm_s::ProductList pro = ver(0).getProducts();
            if (pro.size() == 1)
               rea(0).deleteVertices(1);
         }
      }
      hddm_s::HitViewList vie = event_info->getOutputRecord()->getHitViews();
      if (vie.size() > 0) {
         vie(0).deleteTaggers();
      }
   }

   const G4Step *step = G4ParallelWorldProcess::GetHyperStep();
   double KE = step->GetPostStepPoint()->GetKineticEnergy();
   int iEbin = (KE < gammaPhiXS_Emax)?
                KE / (gammaPhiXS_dE) :
                gammaPhiXS_Emax / gammaPhiXS_dE;
   G4VPhysicalVolume *pvol = step->GetPostStepPoint()->GetPhysicalVolume();
   if (pvol && pvol->GetName() == "KLFT") {
      fPIL = cm / (gammaPhiXS_ub[iEbin] * neutralKaonPhi_bratio *
                   cm2_ub * beryllium_density_gcm3 * NAvagadro);
   }
   else if (pvol && pvol->GetName() == "KLFW") {
      fPIL = cm / (gammaPhiXS_ub[iEbin] * neutralKaonPhi_bratio *
                   cm2_ub * tungsten_density_gcm3 * NAvagadro);
   }
   else {
      fPIL = 1e99;
   }
   *condition = NotForced;
   return -log(G4UniformRand()) * fPIL / gammaPhiXS_scale_factor;
}

G4VParticleChange *GlueXKlongConversionProcess::PostStepDoIt(
                                                  const G4Track &track, 
                                                  const G4Step &step)
{
   pParticleChange->Initialize(track);
   GlueXUserEventInformation *eventinfo;
   const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
   eventinfo = (GlueXUserEventInformation*)event->GetUserInformation();
   double tvtx = step.GetPreStepPoint()->GetGlobalTime();
   G4ThreeVector vtx = step.GetPreStepPoint()->GetPosition();
   G4ThreeVector mom = step.GetPreStepPoint()->GetMomentum();
   G4ThreeVector pol = step.GetPreStepPoint()->GetPolarization();
   eventinfo->AddBeamParticle(1, tvtx, vtx, mom, pol);
   GenerateConversionVertex(track, step);
   pParticleChange->ProposeTrackStatus(fStopAndKill);
   eventinfo->SetKeepEvent(1);
   return pParticleChange;
}

void GlueXKlongConversionProcess::GenerateConversionVertex(const G4Track &track,
                                                             const G4Step &step)
{
   G4ThreeVector mom = step.GetPreStepPoint()->GetMomentum();
   G4ThreeVector pol = step.GetPreStepPoint()->GetPolarization();

   // Generate a new vertex for the klong,kshort pair
   double beamVelocity = GlueXPhotonBeamGenerator::getBeamVelocity();
   double steplength = pParticleChange->GetTrueStepLength();
   G4ThreeVector direction(track.GetMomentumDirection());
   G4ThreeVector x0(track.GetPosition());
   double t0 = track.GetGlobalTime();
   double uvtx = G4UniformRand();
   double lvtx = steplength + fPIL * log(1 - uvtx * (1 - exp(-steplength / fPIL)));
   x0 -= lvtx * direction;
   t0 -= lvtx / beamVelocity;
   G4PrimaryVertex* vertex = new G4PrimaryVertex(x0, t0);

   // Compute the phi kinematics in the gamma,N rest frame
   double KE = track.GetKineticEnergy();
   double mandelS = sqr(KE + massTargetNucleon) - sqr(KE);
   double mandelT = log(G4UniformRand() + 1e-99) / (gammaPhiXS_tslope);
   double EstarPhi = (mandelS + sqr(massPhiMeson) - sqr(massTargetNucleon)) /
                     sqrt(4 * mandelS);
   if (EstarPhi < massPhiMeson) {
      //G4cerr << "Warning in GlueXKlongConversionProcess::GenerateConversionVertex - "
      //       << "phi meson production attempted below threshold, "
      //       << "cannot continue." << G4endl;
      return;
   }
   double qstarPhi = sqrt(sqr(EstarPhi) - sqr(massPhiMeson));
   double qstarInc = sqrt(sqr(mandelS - sqr(massTargetNucleon)) / (4 * mandelS));
   double costhetastarPhi = (mandelT - sqr(qstarInc - EstarPhi) +
                             sqr(qstarInc) + sqr(qstarPhi)) /
                             (2 * qstarInc * qstarPhi);
   double phistarPhi = 2 * M_PI * G4UniformRand();
   if (fabs(costhetastarPhi) > 1) {
      costhetastarPhi /= fabs(costhetastarPhi);
      mandelT = sqr(qstarInc - EstarPhi) - sqr(qstarInc - qstarPhi);
   }
   double sinthetastarPhi = sqrt(1 - sqr(costhetastarPhi));
 
   // Compute the klong kinematics in the phi rest frame
   double EstarKL = massPhiMeson / 2;
   double qstarKL = sqrt(sqr(EstarKL) - sqr(massK0Meson));
   double costhetastarKL = 2 * (G4UniformRand() - 0.5);
   double phistarKL = 2 * M_PI * G4UniformRand();
   if (fabs(costhetastarKL) > 1) {
      costhetastarKL /= fabs(costhetastarKL);
   }
   double sinthetastarKL = sqrt(1 - sqr(costhetastarKL));
   G4LorentzVector klong_p(qstarKL * sinthetastarKL * cos(phistarKL),
                           qstarKL * sinthetastarKL * sin(phistarKL),
                           qstarKL * costhetastarKL, EstarKL);
   G4LorentzVector kshort_p(-qstarKL * sinthetastarKL * cos(phistarKL),
                            -qstarKL * sinthetastarKL * sin(phistarKL),
                            -qstarKL * costhetastarKL, EstarKL);

   // Boost kaons into the gamma,N rest frame
   G4ThreeVector vPhi_reaction(
                 qstarPhi * sinthetastarPhi * cos(phistarPhi) / EstarPhi,
                 qstarPhi * sinthetastarPhi * sin(phistarPhi) / EstarPhi,
                 qstarPhi * costhetastarPhi / EstarPhi);
   klong_p.boost(vPhi_reaction);
   kshort_p.boost(vPhi_reaction);

   // Boost back into the lab frame, then rotate to the align lab axes
   double labbeta = KE / (KE + massTargetNucleon);
   G4ThreeVector vReaction_lab(0, 0, labbeta);
   G4LorentzVector kphoton(0, 0, qstarInc, qstarInc);
   klong_p.boost(vReaction_lab);
   kshort_p.boost(vReaction_lab);
   kphoton.boost(vReaction_lab);
   G4ThreeVector axis = direction.cross(vReaction_lab);
   double alpha = -asin(axis.mag() / labbeta);
   axis /= axis.mag();
   klong_p.rotate(alpha, axis);
   kshort_p.rotate(alpha, axis);
   kphoton.rotate(alpha, axis);
   
#ifdef CHECK_KINEMATICS
   // Check that the results make sense
   G4ThreeVector dphoton = kphoton - mom;
   if (fabs(dphoton.mag()) > 1e-3*GeV) {
      std::cerr << "Warning in GlueXKlongConversionProcess::GenerateConversionVertex - "
                << " transform from phi rest frame to lab failed, gamma momentum=" << mom
                << ", boosted and rotated copy=" << kphoton
                << std::endl;
   }
   G4LorentzVector pPhi = klong_p + kshort_p;
   double mPhi = pPhi.m();
   if (fabs(mPhi - massPhiMeson) > 1e-3*GeV) {
      std::cerr << "Warning in GlueXKlongConversionProcess::GenerateConversionVertex - "
                << " phi mass mismatch, mPhi=" << mPhi 
                << ", massPhiMeson=" << massPhiMeson
                << std::endl;
   }
   G4LorentzVector pbeam(track.GetMomentum(), track.GetKineticEnergy());
   G4LorentzVector pxfer = pPhi - pbeam;
   if (fabs(pxfer.m2() - mandelT) > 1e-3*GeV*GeV) {
      std::cerr << "Warning in GlueXKlongConversionProcess::GenerateConversionVertex - "
                << " generated t mismatch, pxfer.m2()=" << pxfer.m2() 
                <<", mandelT=" << mandelT
                << std::endl;
   }
#endif

   // append secondary vertex to MC record
   G4ParticleDefinition *klong_type = G4KaonZeroLong::Definition();
   G4PrimaryParticle *klong_instance = new G4PrimaryParticle(klong_type,
                      klong_p.px(), klong_p.py(), klong_p.pz(), klong_p.e());
   vertex->SetPrimary(klong_instance);
   G4ParticleDefinition *kshort_type = G4KaonZeroShort::Definition();
   G4PrimaryParticle *kshort_instance = new G4PrimaryParticle(kshort_type,
                      kshort_p.px(), kshort_p.py(), kshort_p.pz(), kshort_p.e());
   vertex->SetPrimary(kshort_instance);

   // create secondaries to continue tracking products
   auto pchange = dynamic_cast<G4ParticleChangeForGamma*>(pParticleChange);
   pchange->InitializeForPostStep(track);
   pchange->AddSecondary(new G4DynamicParticle(klong_type, klong_p));
   pchange->AddSecondary(new G4DynamicParticle(kshort_type, kshort_p));

   GlueXUserEventInformation *event_info;
   const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
   event_info = (GlueXUserEventInformation*)event->GetUserInformation();
   if (event_info) {
      int mech[2];
      char *cmech = (char*)mech;
      snprintf(cmech, 5, "%c%c%c%c", 'K', 'L', 'K', 'S');
      hddm_s::ReactionList rea = event_info->getOutputRecord()->getReactions();
      if (rea.size() > 0) {
         rea(0).setType(372); // peripheral phi production into Kl,Ks
         rea(0).setWeight(1.0);
      }
      event_info->AddPrimaryVertex(*vertex);
      event_info->AddBeamParticle(1, t0, x0, mom, pol);
   }
}
