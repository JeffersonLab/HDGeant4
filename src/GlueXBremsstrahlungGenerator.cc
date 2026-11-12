//
// GlueXBremsstrahlungGenerator class implementation
//
// author: richard.t.jones at uconn.edu
// version: february 4, 2017
//
// In the context of the Geant4 event-level multithreading model,
// this class is "shared", ie. has no thread-local state. It is
// invoked from the GlueXGeneratorAction class, which is
// responsible for managing the interlocks to ensure thread safety.
//
// Notes:
// This is an auxilliary class that is not used in the main production
// version of HDGeant4. It was written to provide an alternate means to
// calculate statistical properties of the bremsstrahlung beam from the
// underlying QED processes. The primary photon beam generator class is
// GlueXPhotonBeamGenerator, which uses analytic distribution formulas
// to generate distributions in photon energy, angle, and polarization.
// By contrast, this GlueXBremsstrahlungGenerator directly computes the
// QED amplitude for a specific choice for the event kinematics, and
// then employs importance sampling to achieve a true reproduction of
// the parent distribution without needing to know any formulas that
// describe them. This latter method has the advantage that one can
// use it to sample any property of the reaction (eg. the polarization
// of the recoil electron, correlation of the photon polarization with
// emission angle, etc.) simply by histogramming it over a sample.

#include <iostream>
#include <stdlib.h>
#include <math.h>

#define sqr(x) ((x)*(x))

#define IMPORTANCE_SAMPLING_HIST_FILE "bsampling_weights.root"

#if USING_DIRACXX

#include <TROOT.h>
#include <TPhoton.h>
#include <TLepton.h>
#include <TLorentzBoost.h>
#include <TCrossSection.h>
#include <Complex.h>
#include <constants.h>

#include <G4ios.hh>

#include <GlueXBremsstrahlungGenerator.hh>

TTree *GlueXBremsstrahlungGenerator::fTree = 0;
TFile *GlueXBremsstrahlungGenerator::fTreeFile = 0;
G4Mutex GlueXBremsstrahlungGenerator::fMutex = G4MUTEX_INITIALIZER;

GlueXBremsstrahlungGenerator::GlueXBremsstrahlungGenerator(TFile *rootfile)
 : fBeamEnergy(12.), fMinEnergy(3.)
{
   for (int i=0; i < 5; i++) {
      fImportSample[i] = 0;
   }

   G4AutoLock barrier(&fMutex);
   if (fTree == 0)
      return;

#if defined IMPORTANCE_SAMPLING_HIST_FILE
   {
      TFile bsamples(IMPORTANCE_SAMPLING_HIST_FILE);
      fImportSample[0] = (TH1D*)bsamples.Get("u0_weight");
      fImportSample[1] = (TH1D*)bsamples.Get("u1_weight");
      fImportSample[2] = (TH1D*)bsamples.Get("u2_weight");
      fImportSample[3] = (TH1D*)bsamples.Get("u3_weight");
      fImportSample[4] = (TH1D*)bsamples.Get("u4_weight");
      bool found = true;
      for (int i=0; i < 5; ++i) {
         if (fImportSample[i]) {
            fImportSample[i]->SetDirectory(0);
            normalize(fImportSample[i]);
         }
         else {
            found = false;
         }
      }
      if (found) {
         std::cout << "GlueXBremsstrahlungGenerator constructor: "
                   << std::endl
                   << "  importance sampling histograms read from "
                   << IMPORTANCE_SAMPLING_HIST_FILE
                   << std::endl;
      }
   }
#endif

   if (rootfile) {
      if (fTreeFile)
         delete fTreeFile;
      fTreeFile = rootfile;
   }
   else if (fTreeFile == 0) {
      fTreeFile = new TFile("beamtree.root", "recreate");
      if (fTree) {
         delete fTree;
         fTree = 0;
      }
   }
   fTreeFile->cd();

   if (fTree == 0) {
      fTree = new TTree("beam", "polarized bremstrahlung generator");
      fTree->Branch("Ebeam", &Ebeam, "Ebeam/D");
      fTree->Branch("qT", &qT, "qT/D");
      fTree->Branch("qTphi", &qTphi, "qTphi/D");
      fTree->Branch("qL", &qL, "qL/D");
      fTree->Branch("kstar", &kstar, "kstar/D");
      fTree->Branch("M2", &M2, "M2/D");
      fTree->Branch("phi", &phi, "phi/D");
      fTree->Branch("u", &u[0], "u[5]/D");
      fTree->Branch("k", &k[0], "k[4]/D");
      fTree->Branch("pR", &pR[0], "pR[4]/D");
      fTree->Branch("polar_0_90", &polar_0_90, "polar_0_90/D");
      fTree->Branch("polar_45_135", &polar_45_135, "polar_45_135/D");
      fTree->Branch("polar_90_0", &polar_90_0, "polar_90_0/D");
      fTree->Branch("polar_135_45", &polar_135_45, "polar_135_45/D");
      fTree->Branch("diffXS", &diffXS, "diffXS/D");
      fTree->Branch("weight", &weight, "weight/D");
   }
}

GlueXBremsstrahlungGenerator::~GlueXBremsstrahlungGenerator()
{
   G4AutoLock barrier(&fMutex);
   if (fTree && fTreeFile) {
      fTree->Write();
      delete fTree;
      fTree = 0;
   }
   if (fTreeFile) {
      delete fTreeFile;
      fTreeFile = 0;
   }
}

void GlueXBremsstrahlungGenerator::GenerateBeamPhotons(int nevents)
{
   // This method implements the generator as a Monte Carlo loop
   // that runs until nevents bremsstrahlung events have been
   // generated between photon energy fMinEnergy and fBeamEnergy.
   // The results are stored in a ROOT tree called "beam".
   //
   // The key input for this generator is the differential cross section
   // for atomic bremsstrahlung that is computed using the TCrossSection
   // class provided by the Dirac++ package. The Bremsstrahlung method
   // of this class takes as its arguments Dirac++ objects describing
   // the incident electron, the recoil electron and the emitted photon,
   // and returns the differential cross section in microbarns/GeV^4/r.
   // The differential measure is d(sigma)/(dk dphi d^3 q) where k is
   // the energy of the bremsstrahlung photon and phi is the azimuthal angle
   // of the photon.  The polar angle of the photon is fixed by kinematics.
   // It is assumed that eIn.Mom()[0] = eOut.Mom()[0]+gOut.Mom()[0], that
   // the energy carried away by the recoil is zero in the laboratory frame,
   // but it is not checked. This result must be multiplied by the square
   // of the form factor of the atomic target and then integrated over
   // recoil kinematics to generate the photon spectrum.
   //
   // The way this is done below is as follows. First the recoil momentum
   // q is split into two transverse components collectively called qT
   // and longitudinal component qL, then qL is transformed into a Lorentz
   // invariant M2 = (P0 + q)^2 where P0 is the incident electron momentum
   // 4-vector and q is the four-vector (0, qTx, qTy, qL) of the momentum
   // exchange with the target atom: qL = -(M2 -mElectron^2 + qT^2)/(2 E0)
   // for a ultra-relativistic incident electron energy E0. The physical
   // meaning of M2 is the invariant mass-squared of the final electron-
   // photon system. Boosting the photon into this system, its lab energy
   // gets translated into an emission angle theta_cm with respect to the
   // boost vector. In combination with phi_cm = phi, generating a random
   // direction in the M2 rest frame becomes equivalent to generating a
   // random energy for the photon in the lab. The advantage of these new
   // coordinates is that their kinematic bounds are now independent from
   // one another. In terms of these new coordinates, the measure is now
   //
   //  dk dphi d^3 q = [(M2 - mElectron^2)/(4 M2)] d^2 qT dM2 dOmega_cm
   //
   // By this approach, the final photon and recoil electron kinematics
   // are generated first in a reaction rest frame, and then boosted into
   // the lab as a final step.

   int good_event;
   for (int event=0; event < nevents; event += good_event) {
      good_event = 0; 
      weight = 1;
      for (int i=0; i < 5; ++i) {
         if (fImportSample[i]) {
            u[i] = fImportSample[i]->GetRandom();
            int bin = fImportSample[i]->FindBin(u[i]);
            weight /= fImportSample[i]->GetBinContent(bin);
         }
         else {
            u[i] = fRandom.Uniform();
         }
      }
      Ebeam = fBeamEnergy;
      double pbeam = sqrt(sqr(Ebeam) - sqr(mElectron));

      // generate u = q0^2 / (q0^2 + qT^2) ~ Uniform[0,1]
      double q0 = 0.2e-3; // 200 keV/c 
      qT = q0 * sqrt(1/u[0] - 1);
      qTphi = u[1] * 2*M_PI;
      weight *= sqr((sqr(q0) + sqr(qT)) / q0);

      // generate u = m0^2 / (m0^2 + M2 -mElectron^2) ~ Uniform[0,1]
      double m0sqr = sqr(1.e-3); // 1 MeV
      M2 = sqr(mElectron) + m0sqr * (1/u[2] - 1);
      qL = sqrt(sqr(Ebeam) - M2 - sqr(qT)) - pbeam;
      weight *= sqr(m0sqr + M2 - sqr(mElectron)) / m0sqr;

      // generate cm angles theta,phi uniform on the unit sphere
      double costheta = 2 * u[3] - 1;
      double sintheta = sqrt(1 - sqr(costheta));
      phi = u[4] * 2*M_PI;
      weight *= 4*M_PI;

      // compute the weight on the chosen measure
      weight *= (M2 - sqr(mElectron)) / (4 * M2);

      // compute the reaction kinematics based on what was generated above
      TFourVectorReal pin(Ebeam, 0, 0, pbeam);
      TFourVectorReal q(0, qT * cos(qTphi), qT * sin(qTphi), qL);
      TThreeVectorReal pcm(sintheta * cos(phi), sintheta * sin(phi), costheta);
      TFourVectorReal plab(pin + q);
      if (plab.InvariantSqr() < M2 - 1e-9) {
         std::cout << "bad kinematics!!" 
                   << " M2 - plab^2 = " << M2 - plab.InvariantSqr()
                   << std::endl;
      }
      TLorentzBoost toLab(plab);
      toLab.Invert();
      kstar = (M2 - sqr(mElectron)) / (2 * sqrt(M2));
      TFourVectorReal kout(kstar, kstar * pcm);
      kout.Boost(toLab);
      TFourVectorReal pout(sqrt(sqr(mElectron) + sqr(kstar)), -kstar * pcm);
      pout.Boost(toLab);

      if (kout[0] < fMinEnergy)
         continue;

      TLepton ein(pin, mElectron);
      TLepton eout(pout, mElectron);
      TPhoton gout(kout);
      TThreeVectorReal zeroVector(0,0,0);
      TThreeVectorReal posXhat(1,0,0);
      TThreeVectorReal posYhat(0,1,0);
      TThreeVectorReal posUhat(1,1,0);
      TThreeVectorReal posVhat(-1,1,0);
      posUhat /= sqrt(2.);
      posVhat /= sqrt(2.);
      ein.SetPol(zeroVector);
      eout.AllPol();
      gout.AllPol();
      double FF2 = sqr(AtomicFormFactor(q.LengthSqr()));
      diffXS = FF2 * TCrossSection::Bremsstrahlung(ein, eout, gout);
      gout.SetPol(posXhat);
      double polar0 = FF2 * TCrossSection::Bremsstrahlung(ein, eout, gout);
      gout.SetPol(posYhat);
      double polar90 = FF2 * TCrossSection::Bremsstrahlung(ein, eout, gout);
      gout.SetPol(posUhat);
      double polar45 = FF2 * TCrossSection::Bremsstrahlung(ein, eout, gout);
      gout.SetPol(posVhat);
      double polar135 = FF2 * TCrossSection::Bremsstrahlung(ein, eout, gout);

      for (int i=0; i < 4; ++i) {
         k[i] = kout[i];
         pR[i] = pout[i];
      }
      polar_0_90 = polar0 / diffXS;
      polar_45_135 = polar45 / diffXS;
      polar_90_0 = polar90 / diffXS;
      polar_135_45 = polar135 / diffXS;
      {
         G4AutoLock barrier(&fMutex);
         fTree->Fill();
         good_event = 1;
         if (event % 1000 == 0)
            G4cout << fTree->GetEntries() << " events written\r"
                   << std::flush;
      }
   }
}

double GlueXBremsstrahlungGenerator::AtomicFormFactor(double q2)
{
   // This is a generic dipole model of a light atom, atomic number Z. 
   // The argument q2 should be in GeV^2, and the result is a
   // dimensionless coefficient that multiplies the charge of
   // the electron. Both electron cloud and nuclear contributions
   // to the atomic electric form factor are included.

   double Z = 13; // aluminum
   double beta = 111 * pow(Z, -1/3.) / mElectron;	// ff cutoff in /GeV
   double Fff = 1 / (1 + sqr(beta) * q2);
   return Z * (1 - Fff);
}

void GlueXBremsstrahlungGenerator::normalize(TH1D *hist)
{
   double sum = hist->Integral();
   double dx = hist->GetBinWidth(1);
   double normfactor = 1 / (sum * dx);
   int nbins = hist->GetNbinsX();
   for (int i=1; i <= nbins; ++i) {
      hist->SetBinContent(i, hist->GetBinContent(i) * normfactor);
   }
}

#endif
