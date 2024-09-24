//
//
// class implementation for GlueXKlongConversionProcess
//
// author: richard.t.jones at uconn.edu
// version: august 2, 2024
//

#include "GlueXKlongConversionProcess.hh"
#include <G4SystemOfUnits.hh>
#include "G4PhysicalConstants.hh"
#include "G4KaonZeroLong.hh"
#include "G4KaonZeroShort.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4Gamma.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4ParallelWorldProcess.hh"

#include <GlueXUserEventInformation.hh>
#include "GlueXUserOptions.hh"

#include <stdio.h>
#include <iomanip>

#define CHECK_KINEMATICS 1

// Phi photoproduction total cross section as function of photon energy,
// digitized from Fig. 3 in Wang et al, arXiv:2208.10289 (22 Aug 2022).

double GlueXKlongConversionModel::gammaPhiXS_dE(0.1*GeV);
double GlueXKlongConversionModel::gammaPhiXS_ub[220] = {
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
double GlueXKlongConversionModel::gammaPhiXS_scale_factor(1);
double GlueXKlongConversionModel::gammaPhiXS_tslope(6.2/(GeV*GeV));
double GlueXKlongConversionModel::gammaPhiXS_Emin(1.9*GeV);
double GlueXKlongConversionModel::gammaPhiXS_Emax(22.0*GeV);

double neutralKaonPhi_bratio(0.342);
double chargedKaonPhi_bratio(0.492);
double massTargetNucleon(0.938*GeV);
double massPhiMeson(1.020*GeV);
double massK0Meson(0.4976*GeV);
double massKpmMeson(0.493677*GeV);
double NAvagadro(6.023e23);
double cm2_ub(1e-30);

// conversion target material properties
double beryllium_Z(4);
double beryllium_A(9);
double beryllium_density_gcm3(1.85);
double tungsten_Z(74);
double tungsten_A(0.85 * 183.8 + 0.15 * 63.5); // 85% tungsten + 15% copper
double tungsten_density_gcm3(1 / ((0.85 / 19.38) + (0.15 / 8.96)));

GlueXKlongConversionProcess::GlueXKlongConversionProcess(
                                       const G4String &name, 
                                       G4ProcessType aType)
 : G4VEmProcess(name, aType),
   isInitialised(false)
{
   SetStartFromNullFlag(false);
   SetBuildTableFlag(true);
   SetProcessSubType(fGammaConversion);
   SetSplineFlag(false);

   verboseLevel = 0;
   // Verbosity for debugging purposes:
   // 0 = nothing 
   // 1 = calculation of cross sections, file openings...
   // 2 = entering in methods

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
      if (! EmModel()) {
         SetEmModel(new GlueXKlongConversionModel);
      }
      AddEmModel(1, EmModel());
   } 
   SetMinKinEnergy(GlueXKlongConversionModel::gammaPhiXS_Emin);
   SetMaxKinEnergy(GlueXKlongConversionModel::gammaPhiXS_Emax);
}

void GlueXKlongConversionProcess::PrintInfo()
{}         



GlueXKlongConversionModel::GlueXKlongConversionModel()
 : G4VEmModel("KlongTargetConversion"),
   isInitialised(false)
{
   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXBeamConversionProcess constructor - "
             << "GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   std::map<int,double> gphi1020pars;
   if (user_opts->Find("GPHI", gphi1020pars))
   {
      if (gphi1020pars.find(1) != gphi1020pars.end())
         gammaPhiXS_scale_factor = gphi1020pars[1];
      if (gphi1020pars.find(2) != gphi1020pars.end())
         gammaPhiXS_tslope = gphi1020pars[2];
      if (gphi1020pars.find(3) != gphi1020pars.end())
         gammaPhiXS_Emin = gphi1020pars[3];
      if (gphi1020pars.find(4) != gphi1020pars.end())
         gammaPhiXS_Emax = gphi1020pars[4];
   }

   fParticleChange = 0;

   verboseLevel= 0;
   // Verbosity for debugging purposes:
   // 0 = nothing 
   // 1 = calculation of cross sections, file openings...
   // 2 = entering in methods
 
   if (verboseLevel > 0) {
      G4cout << "GlueXKlongConversionModel is constructed " << G4endl;
   }
}

GlueXKlongConversionModel::~GlueXKlongConversionModel()
{}

void GlueXKlongConversionModel::Initialise(const G4ParticleDefinition* particle,
                                           const G4DataVector& cuts)
{
   SetLowEnergyLimit(gammaPhiXS_Emin);
   SetHighEnergyLimit(gammaPhiXS_Emax);

   if (verboseLevel > 1) {
      G4cout << "Calling Initialise() of GlueXKlongConversionModel." << G4endl
             << "Energy range: "
             << LowEnergyLimit() / GeV << " GeV - "
             << HighEnergyLimit() / GeV << " GeV"
             << G4endl;
   }

   if (isInitialised) {
      return;
   }
   fParticleChange = GetParticleChangeForGamma();
   isInitialised = true;
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

   double pil = DBL_MAX;
   double kinE = track.GetKineticEnergy();
   if (kinE > GlueXKlongConversionModel::gammaPhiXS_Emin) {
      pil = G4VEmProcess::PostStepGetPhysicalInteractionLength(track, previousStepSize, condition);
   }
   else {
      *condition = NotForced;
      theNumberOfInteractionLengthLeft = -1.0;
      currentInteractionLength = DBL_MAX;
   }
 
   if (verboseLevel > 0) {
       G4cout << "GlueXKlongConversionProcess::PostStepGetPhysicalInteractionLength"
              << " called with energy " << track.GetKineticEnergy()
              << " returns pil=" << pil << G4endl;
   }
   
   return pil;
}

G4double GlueXKlongConversionModel::ComputeCrossSectionPerAtom(
                              const G4ParticleDefinition* photon,
                                    G4double kinEnergy, 
                                    G4double Z, 
                                    G4double A, 
                                    G4double cut,
                                    G4double emax)
{
   int iEbin = (kinEnergy < gammaPhiXS_Emax)?  kinEnergy / (gammaPhiXS_dE) :
                                               gammaPhiXS_Emax / gammaPhiXS_dE;
   G4double xs_atom = 0;
   if (Z == beryllium_Z) {
      xs_atom = gammaPhiXS_ub[iEbin]*microbarn * beryllium_A;
   }
   else if (Z == tungsten_Z) {
      xs_atom = gammaPhiXS_ub[iEbin]*microbarn * tungsten_A;
   }

   if (verboseLevel > 0) {
       G4cout << "GlueXKlongConversionModel::ComputeCrossSectionPerAtom"
              << " called with energy " << kinEnergy
              << " returns xs_atom=" << xs_atom * gammaPhiXS_scale_factor << G4endl;
   }
   return xs_atom * gammaPhiXS_scale_factor;
}

void GlueXKlongConversionModel::SampleSecondaries(
                         std::vector<G4DynamicParticle*>* secondaries,
                         const G4MaterialCutsCouple* material_couple,
                         const G4DynamicParticle* photon,
                         G4double tmin,
                         G4double maxEnergy)
{
   G4ThreeVector mom = photon->GetMomentum();
   G4ThreeVector pol = photon->GetPolarization();
   G4ThreeVector direction(photon->GetMomentumDirection());

   // Decide whether to decay to the neutral or charged channel
   double udecay = G4UniformRand();
   int neutralkaons(0);
   double massKaons;
   if (udecay < neutralKaonPhi_bratio) {
      neutralkaons = 1;
      massKaons = massK0Meson;
   }
   else if (udecay < neutralKaonPhi_bratio + chargedKaonPhi_bratio) {
      neutralkaons = 0;
      massKaons = massKpmMeson;
   }
   else {
      return;
   }

   // Compute the phi kinematics in the gamma,N rest frame
   double KE = photon->GetKineticEnergy();
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
   double qstarKL = sqrt(sqr(EstarKL) - sqr(massKaons));
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
   G4LorentzVector pbeam(photon->Get4Momentum());
   G4LorentzVector pxfer = pPhi - pbeam;
   if (fabs(pxfer.m2() - mandelT) > 1e-3*GeV*GeV) {
      std::cerr << "Warning in GlueXKlongConversionProcess::GenerateConversionVertex - "
                << " generated t mismatch, pxfer.m2()=" << pxfer.m2() 
                <<", mandelT=" << mandelT
                << std::endl;
   }
#endif

   // create secondaries for the decay products
   if (neutralkaons) {
      G4ParticleDefinition *klong_type = G4KaonZeroLong::Definition();
      G4ParticleDefinition *kshort_type = G4KaonZeroShort::Definition();
      secondaries->push_back(new G4DynamicParticle(klong_type, klong_p));
      secondaries->push_back(new G4DynamicParticle(kshort_type, kshort_p));
   }
   else {
      G4ParticleDefinition *kplus_type = G4KaonPlus::Definition();
      G4ParticleDefinition *kminus_type = G4KaonMinus::Definition();
      secondaries->push_back(new G4DynamicParticle(kplus_type, klong_p));
      secondaries->push_back(new G4DynamicParticle(kminus_type, kshort_p));
   }

   // kill off the converted gamma
   fParticleChange->ProposeTrackStatus(fStopAndKill);
   fParticleChange->SetProposedKineticEnergy(0.);
   fParticleChange->ProposeLocalEnergyDeposit(KE);
}
