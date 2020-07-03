//
// class implementation for GlueXBeamConversionProcess
//
// author: richard.t.jones at uconn.edu
// version: december 24, 2016
//

#include "GlueXBeamConversionProcess.hh"
#include "GlueXPhotonBeamGenerator.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"
#include "GlueXUserOptions.hh"
#include "G4ParallelWorldProcess.hh"

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
#define FORCED_LIH2_BETHE_HEITLER 0

#include <G4SystemOfUnits.hh>
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4RunManager.hh"
#include "G4TrackVector.hh"
#include "G4Gamma.hh"

#include <stdio.h>
#include <iomanip>

#ifdef USING_DIRACXX
#include <TLorentzBoost.h>
#include <TPhoton.h>
#include <TLepton.h>
#include <TCrossSection.h>

void unif01(int n, double *u) { G4Random::getTheEngine()->flatArray(n,u); }

#if USE_ADAPTIVE_SAMPLER
#include "AdaptiveSampler.hh"
AdaptiveSampler sampler(6, &unif01, 1);
int sampler_initialized = 0;
#endif

PairConversionGeneration *GlueXBeamConversionProcess::fPairsGeneration = 0;
#endif

ImportanceSampler GlueXBeamConversionProcess::fPaircohPDF;
ImportanceSampler GlueXBeamConversionProcess::fTripletPDF;
G4double GlueXBeamConversionProcess::fBHpair_mass_min=0;

GlueXBeamConversionProcess::GlueXBeamConversionProcess(const G4String &name, 
                                                       G4ProcessType aType)
 : G4VDiscreteProcess(name, aType)
{
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
       user_opts->Find("BEAM", beampars) &&
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
   }

#if USING_DIRACXX
   if (fPairsGeneration == 0)
      fPairsGeneration = new PairConversionGeneration();
#endif

   fPaircohPDF.Pcut = 60;
   fTripletPDF.Pcut = 15;

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

GlueXBeamConversionProcess::GlueXBeamConversionProcess(
                            GlueXBeamConversionProcess &src)
 : G4VDiscreteProcess(src)
{
   fStopBeamBeforeConverter = src.fStopBeamBeforeConverter;
   fStopBeamAfterConverter = src.fStopBeamAfterConverter;
   fStopBeamAfterTarget = src.fStopBeamAfterTarget;
}

GlueXBeamConversionProcess::~GlueXBeamConversionProcess()
{
#ifdef USING_DIRACXX
#if USE_ADAPTIVE_SAMPLER
   if (sampler_initialized)
      sampler.saveState("BHgen_stats.astate");
#endif
#endif
}

GlueXBeamConversionProcess GlueXBeamConversionProcess::operator=(
                           GlueXBeamConversionProcess &src)
{
   GlueXBeamConversionProcess copy(src);
   return copy;
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
      *condition = Forced;
      return 100*cm;
   }
   else if (track.GetTrackID() == 1 && pvol && pvol->GetName() == "LIH2" &&
       (FORCED_LIH2_BETHE_HEITLER || 
        fStopBeamAfterTarget ))
   {
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
      GenerateBetheHeitlerProcess(step);

      if (verboseLevel > 0) {
         G4cout << "GlueXBeamConversionProcess: beam particle stopped"
                << " in the LiH2 target, Bethe-Heitler conversion forced."
                << G4endl;
      }
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
   // to generate Mpair and u1 generates qR. The algorithm succeeds
   // because the mapping u0->Mpair and u1->qR used here is the
   // same as is used in GenerateBeamPairConversion.

#if USING_DIRACXX

   LDouble_t kin = 9.; // GeV
   TPhoton g0;
   TLepton p1(mElectron);
   TLepton e2(mElectron);
   TLepton e3(mElectron);
   TThreeVectorReal p;
   g0.SetMom(p.SetPolar(kin,0,0));
   g0.SetPol(TThreeVector(0,0,0));
   p1.AllPol();
   e2.AllPol();
   e3.AllPol();

   LDouble_t Epos = kin / 2;
   LDouble_t Mmin = 2 * mElectron;
   LDouble_t Mcut = 5e-3; // 5 MeV cutoff parameter
   LDouble_t um0 = 1 + sqr(Mcut / Mmin);
   LDouble_t qRcut = 1e-3; // 1 MeV/c cutoff parameter

   int Nbins = 50;
   fTripletPDF.Psum = 0;
   fPaircohPDF.Psum = 0;
   for (int i0=0; i0 < Nbins; ++i0) {
      LDouble_t u0 = (i0 + 0.5) / Nbins;
      LDouble_t um = pow(um0, u0);
      LDouble_t Mpair = Mcut / sqrt(um - 1);
      LDouble_t k12star2 = sqr(Mpair / 2) - sqr(mElectron);
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
            fPaircohPDF.density.push_back(0);
            fTripletPDF.density.push_back(0);
            fPaircohPDF.integral.push_back(fPaircohPDF.Psum);
            fTripletPDF.integral.push_back(fTripletPDF.Psum);
            continue;
         }
         LDouble_t k12star = sqrt(k12star2);
         LDouble_t q12mag = sqrt(sqr(E12) - sqr(Mpair));
         LDouble_t costhetastar = (Epos - E12 / 2) * 
                                  Mpair / (k12star * q12mag);
         LDouble_t costhetaR = (sqr(Mpair) / 2 + (kin + mElectron) *
                                (E3 - mElectron)) / (kin * qR);
         if (fabs(costhetastar) > 1 || fabs(costhetaR) > 1) {
            fPaircohPDF.density.push_back(0);
            fTripletPDF.density.push_back(0);
            fPaircohPDF.integral.push_back(fPaircohPDF.Psum);
            fTripletPDF.integral.push_back(fTripletPDF.Psum);
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
         fTripletPDF.Psum += fTripletPDF.Pmax = tripXS * weight;
         fPaircohPDF.Psum += fPaircohPDF.Pmax = pairXS * weight;
         fTripletPDF.density.push_back(fTripletPDF.Pmax);
         fPaircohPDF.density.push_back(fPaircohPDF.Pmax);
         fPaircohPDF.integral.push_back(fPaircohPDF.Psum);
         fTripletPDF.integral.push_back(fTripletPDF.Psum);
      }
   }

   LDouble_t du2 = 1. / sqr(Nbins);
   for (int i0=0, index=0; i0 < Nbins; ++i0) {
      for (int i1=0; i1 < Nbins; ++i1, ++index) {
         LDouble_t randvar = i0 + (i1 + 0.5) / Nbins;
         fTripletPDF.randvar.push_back(randvar);
         fPaircohPDF.randvar.push_back(randvar);
         fTripletPDF.density[index] /= fTripletPDF.Psum * du2;
         fPaircohPDF.density[index] /= fPaircohPDF.Psum * du2;
         fTripletPDF.integral[index] /= fTripletPDF.Psum;
         fPaircohPDF.integral[index] /= fPaircohPDF.Psum;
      }
   }
#endif

   fTripletPDF.Pmax = 0;
   fTripletPDF.Psum = 0;
   fPaircohPDF.Pmax = 0;
   fPaircohPDF.Psum = 0;
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

#if defined DO_TRIPLET_IMPORTANCE_SAMPLE || defined DO_PAIRCOH_IMPORTANCE_SAMPLE
   if (fTripletPDF.density.size() == 0) {
      G4cout << "GlueXBeamConversionProcess::GenerateBeamPairConversion:"
             << G4endl
             << "   Setting up cross section tables, please wait... "
             << std::flush;
      prepareImportanceSamplingPDFs();
      G4cout << "finished." << G4endl;
   }
#endif

   TPhoton gIn;
   TLepton p1(mElectron);
   TLepton e2(mElectron);
   TLepton e3(mElectron);
   p1.AllPol();
   e2.AllPol();
   e3.AllPol();
   const G4Track *track = step.GetTrack();
   LDouble_t kin = track->GetKineticEnergy()/GeV;

   // If we are below pair production threshold, do nothing
   const double mTarget = step.GetPreStepPoint()->GetMaterial()->GetA() * 0.932;
   if (kin < 2 * mElectron * (1 + mElectron / mTarget))
      return;

   G4ThreeVector mom(track->GetMomentum());
   TThreeVectorReal mom0(mom[0]/GeV, mom[1]/GeV, mom[2]/GeV);
   gIn.SetMom(mom0);
   G4ThreeVector pol(track->GetPolarization());
   TThreeVectorReal pol0(pol[0], pol[1], pol[2]);
   gIn.SetPlanePolarization(pol0, pol0.Length());

   // Define an angle and axis that rotates zhat into the direction
   // of the incident gamma, so that the generated kinematics is
   // defined with the incident gamma aligned with zhat, and then
   // rotated at the end into the final spatial direction.
   TThreeVectorReal rockaxis(mom0);
   rockaxis.Cross(TThreeVectorReal(0,0,1));
   double rockangle = asin(rockaxis.Length() / mom0.Length());
   rockaxis /= rockaxis.Length();

   while (true) {
      LDouble_t weight = 1;

      // Generate uniform in E+, phi12, phiR
      LDouble_t Epos = kin * G4UniformRand();
      while (Epos < mElectron) {
         Epos = kin * G4UniformRand();
      }
      weight *= kin - mElectron;
      LDouble_t phi12 = 2*M_PI * G4UniformRand();
      weight *= 2*M_PI;
      LDouble_t phiR = 2*M_PI * G4UniformRand();
      weight *= 2*M_PI;

      LDouble_t u0 = G4UniformRand();
      LDouble_t u1 = G4UniformRand();

#if DO_TRIPLET_IMPORTANCE_SAMPLE

      int i = fTripletPDF.search(u1);
      LDouble_t fi = fTripletPDF.density[i];
      LDouble_t ui = fTripletPDF.integral[i];
      LDouble_t ri = fTripletPDF.randvar[i];
      LDouble_t xi = ri - floor(ri);
      LDouble_t dx = (xi > 0.5)? ri - fTripletPDF.randvar[i-1]:
                              fTripletPDF.randvar[i+1] - ri;
      u1 = xi + dx / 2 - (ui - u1) / (fi * dx);
      u0 = (u0 + floor(ri)) * dx;
      weight /= fi;

#elif DO_PAIRCOH_IMPORTANCE_SAMPLE

      int i = fPaircohPDF.search(u1);
      LDouble_t fi = fPaircohPDF.density[i];
      LDouble_t ui = fPaircohPDF.integral[i];
      LDouble_t ri = fPaircohPDF.randvar[i];
      LDouble_t xi = ri - floor(ri);
      LDouble_t dx = (xi > 0.5)? ri - fPaircohPDF.randvar[i-1]:
                              fPaircohPDF.randvar[i+1] - ri;
      u1 = xi + dx / 2 - (ui - u1) / (fi * dx);
      u0 = (u0 + floor(ri)) * dx;
      weight /= fi;

#endif
   
      // Generate Mpair as 1 / (M [M^2 + Mcut^2])
      LDouble_t Mmin = 2 * mElectron;
      LDouble_t Mcut = 0.005;  // GeV
      LDouble_t um0 = 1 + sqr(Mcut / Mmin);
      LDouble_t um = pow(um0, u0);
      LDouble_t Mpair = Mcut / sqrt(um - 1 + 1e-99);
      weight *= Mpair * (sqr(Mcut) + sqr(Mpair)) * log(um0) / (2 * sqr(Mcut));
   
      // Generate qR^2 with weight 1 / [qR^2 sqrt(qRcut^2 + qR^2)]
      LDouble_t qRmin = sqr(Mpair) /(2 * kin);
      LDouble_t qRcut = 1e-3; // GeV
      LDouble_t uq0 = qRmin / (qRcut + sqrt(sqr(qRcut) + sqr(qRmin)));
      LDouble_t uq = pow(uq0, u1);
      LDouble_t qR = 2 * qRcut * uq / (1 - sqr(uq));
      LDouble_t qR2 = qR * qR;
      weight *= qR2 * sqrt(1 + qR2 / sqr(qRcut)) * (-2 * log(uq0));
   
      // Include overall measure Jacobian factor
      weight *= Mpair / (2 * kin);
   
      // Generate with importance sampling
      LDouble_t Striplet = fTripletPDF.Npassed * 
                        (fTripletPDF.Ntested / (fTripletPDF.Psum + 1e-99));
      LDouble_t Spaircoh = fPaircohPDF.Npassed *
                        (fPaircohPDF.Ntested / (fPaircohPDF.Psum + 1e-99));
      if (Striplet < Spaircoh) {                     // try incoherent generation
         ++fTripletPDF.Ntested;
   
         // Solve for the c.m. momentum of e+ in the pair 1,2 rest frame
         LDouble_t k12star2 = sqr(Mpair / 2) - sqr(mElectron);
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
         if (Epos > E12 - mElectron) {
            // no kinematic solution because Epos > E12 - mElectron, try again
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
         if (q2.Length() < qR) {
            // recoil/pair electrons switched, try again
            continue;
         }

         // Rotate final-state momenta into the lab frame
         p1.SetMom(q1.Rotate(rockaxis, rockangle));
         e2.SetMom(q2.Rotate(rockaxis, rockangle));
         e3.SetMom(q3.Rotate(rockaxis, rockangle));

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
         if (Pfactor > fTripletPDF.Pmax)
            fTripletPDF.Pmax = Pfactor;
         if (Pfactor > fTripletPDF.Pcut) {
            G4cout << "Warning in GenerateBeamPairConversion - Pfactor " 
                   << Pfactor << " exceeds fTripletPDF.Pcut = " 
                   << fTripletPDF.Pcut << G4endl
                   << "  present qR = " << qR << G4endl
                   << "  present Mpair = " << Mpair << G4endl
                   << "  present maximum Pfactor = "
                   << fTripletPDF.Pmax << G4endl
                   << "  current generator efficiency = "
                   << fTripletPDF.Npassed /
                      (fTripletPDF.Ntested + 1e-99)
                   << G4endl;
         }
         fTripletPDF.Psum += Pfactor;
         if (G4UniformRand() * fTripletPDF.Pcut > Pfactor) {
            continue;
         }
         ++fTripletPDF.Npassed;
         break;
      }

      else {                          // try coherent generation
         ++fPaircohPDF.Ntested;
   
         // Solve for the c.m. momentum of e+ in the pair 1,2 rest frame
         // assuming that the atomic target absorbs zero energy
         LDouble_t k12star2 = sqr(Mpair / 2) - sqr(mElectron);
         if (k12star2 < 0) {
            // try again (this should never happen!)
            continue;
         }
         LDouble_t k12star = sqrt(k12star2);
         LDouble_t Eele = kin - Epos;
         if (kin < Mpair) {
            // no kinematic solution because kin < Mpair, try again
            continue;
         }
         else if (Eele < mElectron) {
            // no kinematic solution because Eele < mElectron, try again
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
         p1.SetMom(q1.Rotate(rockaxis, rockangle));
         e2.SetMom(q2.Rotate(rockaxis, rockangle));
         e3.SetMom(TThreeVectorReal(0,0,0));
         LDouble_t diffXS = fPairsGeneration->DiffXS_pair(gIn, p1, e2);
   
         // Use keep/discard algorithm
         LDouble_t Pfactor = diffXS * weight;
         if (Pfactor > fPaircohPDF.Pmax)
            fPaircohPDF.Pmax = Pfactor;
         if (Pfactor > fPaircohPDF.Pcut) {
            G4cout << "Warning in GenerateBeamPairConversion - Pfactor " 
                   << Pfactor << " exceeds fPaircohPDF.Pcut = " 
                   << fPaircohPDF.Pcut << G4endl
                   << "  present qR = " << qR << G4endl
                   << "  present Mpair = " << Mpair << G4endl
                   << "  present maximum Pfactor = "
                   << fPaircohPDF.Pmax << G4endl
                   << "  current generator efficiency = "
                   << fPaircohPDF.Npassed /
                      (fPaircohPDF.Ntested + 1e-99)
                   << G4endl;
         }
         fPaircohPDF.Psum += Pfactor;
         if (G4UniformRand() * fPaircohPDF.Pcut > Pfactor) {
            continue;
         }
         ++fPaircohPDF.Npassed;
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
   x0 -= uvtx * steplength * direction;
   t0 -= uvtx * steplength / beamVelocity;
   G4PrimaryVertex vertex(x0, t0);
   G4ParticleDefinition *positron = G4Positron::Definition();
   G4ParticleDefinition *electron = G4Electron::Definition();
   G4ThreeVector psec1(p1.Mom()[1]*GeV, p1.Mom()[2]*GeV, p1.Mom()[3]*GeV);
   G4ThreeVector psec2(e2.Mom()[1]*GeV, e2.Mom()[2]*GeV, e2.Mom()[3]*GeV);
   G4ThreeVector psec3(e3.Mom()[1]*GeV, e3.Mom()[2]*GeV, e3.Mom()[3]*GeV);
   G4DynamicParticle *sec1 = new G4DynamicParticle(positron, psec1);
   G4DynamicParticle *sec2 = new G4DynamicParticle(electron, psec2);
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
      snprintf(cmech, 5, "%c%c%c%c", 'C', 'O', 'N', 'V');
      event_info->AddSecondaryVertex(secondaries, 1, mech[0]);
   }

#if VERBOSE_PAIRS_SPLITTING
   if ((fTripletPDF.Npassed + fPaircohPDF.Npassed) % 500 == 0) {
      G4cout << std::setprecision(5)
             << "triplet cross section is " << fTripletPDF.Pcut *
             fTripletPDF.Npassed / (fTripletPDF.Ntested + 1e-99) 
             << " +/- " << fTripletPDF.Pcut *
             sqrt(fTripletPDF.Npassed) / (fTripletPDF.Ntested + 1e-99) 
             << " barns, efficiency is " 
             << fTripletPDF.Npassed / (fTripletPDF.Ntested + 1e-99)
             << G4endl
             << "pair cross section is " << fPaircohPDF.Pcut *
             fPaircohPDF.Npassed / (fPaircohPDF.Ntested + 1e-99) 
             << " +/- " << fPaircohPDF.Pcut *
             sqrt(fPaircohPDF.Npassed) / (fPaircohPDF.Ntested + 1e-99) 
             << " barns, efficiency is " 
             << fPaircohPDF.Npassed / (fPaircohPDF.Ntested + 1e-99) 
             << G4endl
             << "counts are "
             << fTripletPDF.Npassed << " / " << fPaircohPDF.Npassed
             << " = "
             << fTripletPDF.Npassed / (fPaircohPDF.Npassed + 1e-99)
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
   // process. The incident gamma ray is stopped and the e+/e- pair vertex
   // is added to the output event record.

#ifndef USING_DIRACXX

   G4cerr << "GlueXBeamConversionProcess::GenerateBetheHeitlerProcess error:"
          << G4endl
          << "  You have enabled Bethe-Heitler conversion in the LIH2 target,"
          << G4endl
          << "  but your have built HDGeant4 without support for the Dirac++"
		  << G4endl
          << "  library. Either rebuild with Dirac++ support or else turn off"
		  << G4endl
          << "  this process. Aborting this run..."
		  << G4endl;
   exit(1);

#else

#if USE_ADAPTIVE_SAMPLER
   if (sampler_initialized == 0) {
      sampler.restoreState("BHgen.astate");
      sampler.setVerbosity(2);
      sampler.reset_stats();
      sampler_initialized = 1;
   }
#endif

   TPhoton gIn;
   TLepton eOut(mElectron);
   TLepton pOut(mElectron);
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

   // If we are below pair production threshold, do nothing
   if (kin < 2 * mElectron * (1 + mElectron / mProton))
      return;

   G4ThreeVector mom(track->GetMomentum());
   TThreeVectorReal mom0(mom[0]/GeV, mom[1]/GeV, mom[2]/GeV);
   gIn.SetMom(mom0);
   G4ThreeVector pol(track->GetPolarization());
   TThreeVectorReal pol0(pol[0], pol[1], pol[2]);
   gIn.SetPlanePolarization(pol0, pol0.Length());

   // Define an angle and axis that rotates zhat into the direction
   // of the incident gamma, so that the generated kinematics is
   // defined with the incident gamma aligned with zhat, and then
   // rotated at the end into the final spatial direction.
   TThreeVectorReal rockaxis(mom0);
   rockaxis.Cross(TThreeVectorReal(0,0,1));
   double rockangle = asin(rockaxis.Length() / mom0.Length());
   rockaxis /= rockaxis.Length();

   // The beam photon energy is already generated by the cobrems
   // generator. Supply a normalized beam energy as the leading
   // elemnt of the u-vector in the call to sample(), and let it
   // return the remaining 5 random numbers needed to fix the
   // kinematics of the pair-production reaction.
   double E0_GeV = GlueXPrimaryGeneratorAction::GetInstance()->
                   GetCobremsGeneration()->getBeamEnergy();
   double u[6];
   u[0] = kin / E0_GeV;
#if USE_ADAPTIVE_SAMPLER
   double weight = sampler.sample(u);
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
   LDouble_t Mthresh = 2 * mElectron;
   LDouble_t Mcut = 0.005;  // GeV
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
   LDouble_t qRcut = 1e-3; // GeV
   LDouble_t uq0 = qRmin / (qRcut + sqrt(sqr(qRcut) + sqr(qRmin)));
   LDouble_t uq = pow(uq0, u[5]);
   LDouble_t qR = 2 * qRcut * uq / (1 - sqr(uq));
   LDouble_t qR2 = qR * qR;
   weight *= qR2 * sqrt(1 + qR2 / sqr(qRcut)) * (-2 * log(uq0));

   // Include overall measure Jacobian factor
   weight *= Mpair / (2 * kin);

   try {
      if (Epos < mElectron) {
         throw std::runtime_error("positron energy less than its rest mass.");
      }

      // Solve for the c.m. momentum of e+ in the pair 1,2 rest frame
      LDouble_t k12star2 = sqr(Mpair / 2) - sqr(mElectron);
      if (k12star2 < 0) {
         throw std::runtime_error("no kinematic solution because k12star2 < 0");
      }
      LDouble_t k12star = sqrt(k12star2);
      LDouble_t Erec = sqrt(qR2 + sqr(mProton));
      LDouble_t E12 = kin + mProton - Erec;
      if (E12 < Mpair) {
         throw std::runtime_error("no kinematic solution because E12 < Mpair");
      }
      LDouble_t q12mag = sqrt(sqr(E12) - sqr(Mpair));
      LDouble_t costhetastar = (Epos - E12 / 2) * Mpair / (k12star * q12mag);
      if (Epos > E12 - mElectron) {
         throw std::runtime_error("no kinematic solution because Epos > E12 - mElectron");
      }
      else if (fabs(costhetastar) > 1) {
         throw std::runtime_error("no kinematic solution because |costhetastar| > 1");
      }

      // Solve for the recoil nucleon kinematics
      LDouble_t costhetaR = (sqr(Mpair) / 2 + (kin + mProton) *
                            (Erec - mProton)) / (kin * qR);
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
      pOut.SetMom(q1.Rotate(rockaxis, rockangle));
      eOut.SetMom(q2.Rotate(rockaxis, rockangle));
      nOut.SetMom(q3.Rotate(rockaxis, rockangle));

      // Check 4-momentum conservation
      TFourVectorReal pIn(gIn.Mom() + nIn.Mom());
      TFourVectorReal pFi(pOut.Mom() + eOut.Mom() + nOut.Mom());
      TFourVectorReal::SetResolution(1e-10);
      if (pIn != pFi) {
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
      }

      // Compute the polarized differential cross section (barnes/GeV^4)
      // returned as d(sigma)/(dE+ dphi+ d^3qR)
      LDouble_t t = sqr(Erec - mProton) - qR2;
      LDouble_t tau = -t / sqr(2 * mProton);
      LDouble_t proton_magnetic_moment = 2.793;
      LDouble_t F1_timelike = 1;
      LDouble_t F2_timelike = 0;
      LDouble_t F1_spacelike = (1 / sqr(1 - t/0.71)) / (1 + tau) *
                               (1 + proton_magnetic_moment * tau);
      LDouble_t F2_spacelike = (1 / sqr(1 - t/0.71)) / (1 + tau) * 
                               (proton_magnetic_moment - 1);
      diffXS = TCrossSection::BetheHeitlerNucleon(gIn, nIn,
                                                  eOut, pOut, nOut,
                                                  F1_spacelike,
                                                  F2_spacelike,
                                                  F1_timelike,
                                                  F2_timelike);
#if USE_ADAPTIVE_SAMPLER
      sampler.feedback(u, weight * diffXS);
#endif
   }
   catch (const std::exception &e) {

      // These events have no cross section, but do not discard
      // them because they are needed to get the right MC integral.
 
      nOut.SetMom(TThreeVectorReal(0,0,1e-12));
      eOut.SetMom(TThreeVectorReal(0,0,1e-12));
      pOut.SetMom(TThreeVectorReal(0,0,1e-12));
#if USE_ADAPTIVE_SAMPLER
      sampler.feedback(u, 0);
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
   G4ParticleDefinition *proton = G4Proton::Definition();
   G4ThreeVector psec1(pOut.Mom()[1]*GeV, pOut.Mom()[2]*GeV, pOut.Mom()[3]*GeV);
   G4ThreeVector psec2(eOut.Mom()[1]*GeV, eOut.Mom()[2]*GeV, eOut.Mom()[3]*GeV);
   G4ThreeVector psec3(nOut.Mom()[1]*GeV, nOut.Mom()[2]*GeV, nOut.Mom()[3]*GeV);
   G4DynamicParticle *sec1 = new G4DynamicParticle(positron, psec1);
   G4DynamicParticle *sec2 = new G4DynamicParticle(electron, psec2);
   G4DynamicParticle *sec3 = new G4DynamicParticle(proton, psec3);
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
   ++passes;
   int tens = 1000;
   while (passes >= tens*10)
      tens *= 10;
   if (passes % tens == 0) {
      std::cout << "sampler reports efficiency " << sampler.getEfficiency()
                << std::endl;
      sampler.saveState("BHgen_stats.astate");
   }
#endif

#endif
}
