//
// class implementation for GlueXPhotonBeamGenerator
//
// author: richard.t.jones at uconn.edu
// version: december 24, 2016

#include "GlueXPhotonBeamGenerator.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserOptions.hh"
#include "HddmOutput.hh"
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>

#include <JANA/JApplication.h>
#include <JANA/JCalibration.h>

// This flag allows HDGeant4 to be used as a Monte Carlo event
// generator of a coherent bremsstrahlung beam. When running in
// "generate not simulate" mode, the output hddm file is filled
// with MC information describing photons incident on the 
// primary collimator, ie. before collimation. This mode is
// enabled when the flag GENBEAM 'precol' is present together
// with the BEAM card in the control.in file. See comments at
// the head of GlueXBeamConversionProcess.cc for other options
// related to the GENBEAM control.in flag.

int GlueXPhotonBeamGenerator::fGenerateNotSimulate = 0;
int GlueXPhotonBeamGenerator::fBeamBackgroundTagOnly = 0;

double GlueXPhotonBeamGenerator::fBeamBucketPeriod = 4. * ns;
double GlueXPhotonBeamGenerator::fBeamStartZ = -24 * m;
double GlueXPhotonBeamGenerator::fBeamDiameter = 0.5 * cm;
double GlueXPhotonBeamGenerator::fBeamVelocity = 2.99792e8 * m/s;
double GlueXPhotonBeamGenerator::fBeamOffset[2] = {0,0};

int GlueXPhotonBeamGenerator::fForceFixedPolarization = false;
double GlueXPhotonBeamGenerator::fFixedPolarization = 0;
double GlueXPhotonBeamGenerator::fFixedPolarization_phi = 0;

// This utility function is useful in the debugger,
// but do not use it for actual simulation operatons.

int g4random_seeds(CLHEP::HepRandomEngine *eng=0) {
   const long int *seed = G4Random::getTheSeeds();
   G4cout << "current seeds are " << seed[0] << ", " << seed[1] << G4endl;
   printf("seeds array pointer is %lx\n", (unsigned long)seed);
   CLHEP::HepRandomEngine *engine = CLHEP::HepRandom::getTheEngine();
   printf("randoms generator engine pointer is %lx\n", (unsigned long)engine);
   if (eng != 0)
      CLHEP::HepRandom::setTheEngine(eng);
   return 0;
}

GlueXPhotonBeamGenerator::GlueXPhotonBeamGenerator(CobremsGeneration *gen)
 : fCobrems(gen),
   fTagger(0)
{
   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXPhotonBeamGenerator constructor - "
             << "GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   std::map<int,std::string> infile;
   std::map<int,double> beampars;
   std::map<int,std::string> genbeampars;
   if (!user_opts->Find("INFILE", infile) &&
       user_opts->Find("BEAM", beampars) &&
       user_opts->Find("GENBEAM", genbeampars))
   {
      if (genbeampars.find(1) != genbeampars.end() && 
         (genbeampars[1] == "precol" ||
          genbeampars[1] == "PRECOL" ||
          genbeampars[1] == "Precol" ||
          genbeampars[1] == "PreCol" ))
      {
         fGenerateNotSimulate = 1;
      }
      else {
         fGenerateNotSimulate = -1;
      }
      GlueXUserEventInformation::setWriteNoHitEvents(1);
   }
   std::map<int, int> bgratepars;
   std::map<int, int> bggatepars;
   std::map<int, int> bgtagonlypars;
   if (user_opts->Find("BEAM", beampars) &&
       user_opts->Find("BGRATE", bgratepars) &&
       user_opts->Find("BGGATE", bggatepars))
   {
      if (user_opts->Find("BGTAGONLY", bgtagonlypars)) {
         fBeamBackgroundTagOnly = bgtagonlypars[1];
      }
      else {
         fBeamBackgroundTagOnly = 0;
      }
   }

   prepareImportanceSamplingPDFs();

   // Create interface for interactive commands

   fMessenger = new G4GenericMessenger(this, "/PhotonBeam/",
                                       "Photon beam generator control");
   fMessenger->DeclareMethod("enableFixedPolarization", 
                             &GlueXPhotonBeamGenerator::enableFixedPolarization,
       "Tell the photon beam generator to force a fixed polarization\n"
       " on all beam photons created by the simulation.\n"
       " Arguments are:\n"
       "      polarization: value in the range [0,1]\n"
       "      phi: azimuthal angle (degrees) [0, 180]");
   fMessenger->DeclareMethod("disableFixedPolarization",
                             &GlueXPhotonBeamGenerator::disableFixedPolarization,
       "Tell the photon beam generator not to force a fixed polarization\n"
       " on beam photons created by the simulation.");
}

GlueXPhotonBeamGenerator::~GlueXPhotonBeamGenerator()
{
   delete fMessenger; 
}

void GlueXPhotonBeamGenerator::prepareImportanceSamplingPDFs()
{
   // Construct lookup tables representing the PDFs used for
   // importance-sampling the coherent bremsstrahlung kinematics.

   const int Ndim = 500;
   double Emin = fCobrems->getPhotonEnergyMin() * GeV;
   double Emax = fCobrems->getBeamEnergy() * GeV;

   // Compute approximate PDF for dNc/dx
   double xmin = Emin / Emax;
   double dx = (1 - xmin) / Ndim;
   double xarr[Ndim], yarr[Ndim];
   double ymax = 0;
   for (int i=0; i < Ndim; ++i) {
      xarr[i] = xmin + (i + 0.5) * dx;
      yarr[i] = fCobrems->Rate_dNcdx(xarr[i]);
      yarr[i] = (yarr[i] > 0)? yarr[i] : 0;
      ymax = (yarr[i] > ymax)? yarr[i] : ymax;
   }
   fCobrems->applyBeamCrystalConvolution(Ndim, xarr, yarr);
   double pcut = 5;
   double psum = 0;
   for (int i=0; i < Ndim; ++i) {
      double yplusbl = yarr[i] + ymax / pcut;
      psum += yplusbl;
      fCoherentPDFx.randvar.push_back(xarr[i]);
      fCoherentPDFx.density.push_back(yplusbl);
      fCoherentPDFx.integral.push_back(psum);
   }
   for (int i=0; i < Ndim; ++i) {
      fCoherentPDFx.density[i] /= psum * dx;
      fCoherentPDFx.integral[i] /= psum;
   }
   fCoherentPDFx.Pcut = 3 * psum / Ndim;

   // Compute approximate PDF for dNi/dlogx
   double logxmin = log(xmin);
   double dlogx = -logxmin / Ndim;
   double dNidlogx;
   double qsum = 0;
   for (int i=0; i < Ndim; ++i) {
      double logx = logxmin + (i + 0.5) * dlogx;
      double x = exp(logx);
      dNidlogx = fCobrems->Rate_dNidxdt2(x, 0) * x;
      dNidlogx = (dNidlogx > 0)? dNidlogx : 0;
      qsum += dNidlogx;
      fIncoherentPDFlogx.randvar.push_back(logx);
      fIncoherentPDFlogx.density.push_back(dNidlogx);
      fIncoherentPDFlogx.integral.push_back(qsum);
   }
   for (int i=0; i < Ndim; ++i) {
      fIncoherentPDFlogx.density[i] /= qsum * dlogx;
      fIncoherentPDFlogx.integral[i] /= qsum;
   }
   fIncoherentPDFlogx.Pcut = qsum / Ndim;
 
   // Compute approximate PDF for dNi/dy
   fIncoherentPDFtheta02 = 1.8;
   double ymin = 1e-3;
   double dy = (1 - ymin) / Ndim;
   double dNidxdy;
   double tsum = 0;
   for (int i=0; i < Ndim; ++i) {
      double y = ymin + (i + 0.5) * dy;
      double theta2 = fIncoherentPDFtheta02 * (1 / y - 1);
      dNidxdy = fCobrems->Rate_dNidxdt2(0.5, theta2) *
                fIncoherentPDFtheta02 / (y*y);
      dNidxdy = (dNidxdy > 0)? dNidxdy : 0;
      tsum += dNidxdy;
      fIncoherentPDFy.randvar.push_back(y);
      fIncoherentPDFy.density.push_back(dNidxdy);
      fIncoherentPDFy.integral.push_back(tsum);
   }
   for (int i=0; i < Ndim; ++i) {
      fIncoherentPDFy.density[i] /= tsum * dy;
      fIncoherentPDFy.integral[i] /= tsum;
   }
   fIncoherentPDFy.Pcut = 1.5 * tsum / Ndim;
   fCoherentPDFx.Pmax = 0;
   fCoherentPDFx.Psum = 0;
   fIncoherentPDFlogx.Pmax = 0;
   fIncoherentPDFlogx.Psum = 0;
   fIncoherentPDFy.Pmax = 0;
   fIncoherentPDFy.Psum = 0;
}

void GlueXPhotonBeamGenerator::GeneratePrimaryVertex(G4Event* anEvent)
{
   GenerateBeamPhoton(anEvent, 0);
}

void GlueXPhotonBeamGenerator::GenerateBeamPhoton(G4Event* anEvent, double t0)
{
   // Generates a single beam photon according to the coherent bremsstrahlung
   // model defined by class CobremsGeneration.  The photon begins its lifetime
   // just upstream of the primary collimator (WARNING: position is hard-wired
   // in the code below) and is tracked by the simulation from there forward.
   // Its time t0 should identify its beam bucket, ie. the time the photon
   // would reach the midplane of the target. To enable beam motion spreading,
   // define the beam box size below.

   // The algorithm below generates coherent bremsstrahlung photons using a
   // importance-sampling technique. This algorithm requires that we prepare
   // an approximate probability density function for the generated photons.
   // The function is not in general equal to the true physical PDF, which
   // varies from event to event depending on the direction of the incident
   // electron, and also the exact angle of the crystal planes at the point
   // of scattering which moves because of the mosaic spread of the crystal.
   // The important thing is that the approximate PDF be reasonably close to
   // the average over all beam particles and the crystal mosaic, and that
   // deviations from event to event are sufficiently small that rejection
   // sampling can be used to take them into account with high efficiency.
   //
   // The kinematics of bremsstrahlung are described by three independent
   // variables (x, theta, phi) where x is the photon energy in units of the
   // incident electron energy, and theta,phi are the polar,azimuthal angles
   // of the photon in a lab frame tilted so that the incident electron comes
   // in along the z axis. Polar angle theta is represented by dimensionless
   // variable y = theta0^2 / (theta^2 + theta0^2) where contant theta0 is
   // chosen to optimize the uniformity of the PDF in y. On each event,
   // a new random tuple (x, phi, y) is generated on the interval x:[0,1],
   // phi:[0,2pi], y:[0,1] using a split-and-recombine strategy. One side 
   // of the split covers the coherent process and the other side covers the
   // incoherent one.
   //
   //  1) coherent process - the PDF here is continuous in x,phi according
   //     the dNc/(dx dphi), and the dependence on y is a sequence of delta 
   //     functions corresponding to the different planes that contribute to
   //     the scattering at the given value of x. Here we take advantage of
   //     the fact that the marginal distribution dNc/dx is proportional to
   //     dNc/(dx dphi) at phi=pi/4. This allows us to decompose the generation
   //     into two stages, first generating x from dNc/dx and then generating
   //     phi from dNc/(dx dphi) at fixed x. The x generation step is performed
   //     using importance sampling based on the average PDF stored in table
   //     fCoherentPDF, followed by rejection sampling based on the value of
   //     dNc/(dx dphi) computed for the particular kinematics of each event.
   //     The y value is obtained by sampling the weighted list of q2 values
   //     that contributed the to q-sum in the calculation of dNc/(dx dphi).
   //
   //  2) incoherent process - the PDF here is continuous in x,phi,y but it
   //     is uniform in phi, so it is effectively a 2D distribution. Here we
   //     take advantage of the fact that x and y are independent variables
   //     to a good approximation, which allows us to generate x using
   //     importance sampling from an approximation to dNi/(dx dtheta^2) at
   //     theta=0 and y ~ uniform [0,1], then employ rejection sampling based
   //     on the exact PDF dNi/(dx dtheta2) to get a true sample.
   //
   // Recombination after the split is very simple. First we integrate over
   // phi in both cases to obtain values dNc/dx and dNi/(dx dy). It turns
   // out that in spite of the fact that the y-dependence is discrete in the
   // coherent case and continuous in the incoherent case, the sum over the
   // probabilities for all values of y in dNc/dx is always normalized to 1
   // independently for all values of x. Hence we can treat y as a psuedo
   // coordinate y' ~ Unif[0,1] and form a 2D PDF dNc/(dx dy') which is
   // numerically equal to dNc/dx, do the rejection sampling in combination
   // with that applied to dNi/(dx dy) and then replace the fake variable y'
   // with the true y that was sampled as described above.

   assert (fCoherentPDFx.integral.size() != 0);

   // The GlueXUserEventInformation constructor can set the random number
   // seeds for this event, so this must happen at here at the top.
   GlueXUserEventInformation *event_info;
   if (t0 == 0) {
      event_info = new GlueXUserEventInformation();
      anEvent->SetUserInformation(event_info);
   }
   else {
      event_info = (GlueXUserEventInformation*)anEvent->GetUserInformation();
   }
   assert (event_info != 0);

   double phiMosaic = 2*M_PI * G4UniformRand();
   double rhoMosaic = sqrt(-2 * log(G4UniformRand()));
   rhoMosaic *= fCobrems->getTargetCrystalMosaicSpread() * radian;
   double thxMosaic = rhoMosaic * cos(phiMosaic);
   double thyMosaic = rhoMosaic * sin(phiMosaic);

   double xemittance = fCobrems->getBeamEmittance() * m*radian;
   double yemittance = xemittance / 2.5; // nominal, should be checked
   double xspotsize = fCobrems->getCollimatorSpotrms() * m;
   double yspotsize = xspotsize; // nominal, should be checked
   double phiBeam = 2*M_PI * G4UniformRand();
   double rhoBeam = sqrt(-2 * log(G4UniformRand()));
   double thxBeam = (xemittance / xspotsize) * rhoBeam * cos(phiBeam);
   double thyBeam = (yemittance / yspotsize) * rhoBeam * sin(phiBeam);

   double raddz = fCobrems->getTargetThickness() * m;
   double varMS = fCobrems->Sigma2MS(raddz/m * G4UniformRand());
   double rhoMS = sqrt(-2 * varMS * log(G4UniformRand()));
   double phiMS = 2*M_PI * G4UniformRand();
   double thxMS = rhoMS * cos(phiMS);
   double thyMS = rhoMS * sin(phiMS);

   double targetThetax = fCobrems->getTargetThetax() * radian;
   double targetThetay = fCobrems->getTargetThetay() * radian;
   double targetThetaz = fCobrems->getTargetThetaz() * radian;
   double thetax = thxBeam + thxMS - targetThetax - thxMosaic;
   double thetay = thyBeam + thyMS - targetThetay - thyMosaic;
   double thetaz = -targetThetaz;
   fCobrems->setTargetOrientation(thetax/radian, thetay/radian, thetaz/radian);

   // Generate with importance sampling
   double x = 0;
   double phi = 0;
   double theta2 = 0;
   double polarization = 0;
   double polarization_phi = 0;
   double coherent_fraction = fCoherentPDFx.Pcut;
   coherent_fraction /= fCoherentPDFx.Pcut + fIncoherentPDFy.Pcut;
   while (true) {
      double splitrand = G4UniformRand();
      if (splitrand < coherent_fraction) { // try coherent generation
         ++fCoherentPDFx.Ntested;

         double u = G4UniformRand();
         int i = fCoherentPDFx.search(u);
         double ui = fCoherentPDFx.integral[i];
         double u_i = (i > 0)? fCoherentPDFx.integral[i-1] : 0;
         double xi = fCoherentPDFx.randvar[i];
         double dx = (i > 0)? xi - fCoherentPDFx.randvar[i-1]:
                              fCoherentPDFx.randvar[i+1] - xi;
         x = xi + dx * (0.5 - (ui - u) / (ui - u_i));
         assert (x > 0);
         double dNcdxPDF = (ui - u_i) / dx;
         double dNcdx = fCobrems->Rate_dNcdx(x);
         double Pfactor = dNcdx / dNcdxPDF;
         if (Pfactor > fCoherentPDFx.Pmax)
            fCoherentPDFx.Pmax = Pfactor;
         if (Pfactor > fCoherentPDFx.Pcut) {
            G4cout << "Warning in GenerateBeamPhoton - Pfactor " << Pfactor
                   << " exceeds fCoherentPDFx.Pcut = " << fCoherentPDFx.Pcut
                   << G4endl
                   << "  present x = " << x << G4endl
                   << "  present maximum Pfactor = "
                   << fCoherentPDFx.Pmax << G4endl
                   << "  current generator efficiency = "
                   << fCoherentPDFx.Npassed /
                      (fCoherentPDFx.Ntested + 1e-99)
                   << G4endl;
         }
         ++fCoherentPDFx.Psum += Pfactor;
         if (G4UniformRand() * fCoherentPDFx.Pcut > Pfactor) {
            continue;
         }
         ++fCoherentPDFx.Npassed;

         double freq;
         double fmax = dNcdx / M_PI;
         while (true) {
            phi = 2*M_PI * G4UniformRand();
            freq = fCobrems->Rate_dNcdxdp(x, phi);
            assert (freq < fmax);
            if (G4UniformRand() * fmax < freq)
               break;
         }
         double uq = freq * G4UniformRand();
         int j = ImportanceSampler::search(uq, fCobrems->fQ2weight);
         theta2 = fCobrems->fQ2theta2[j];
         polarization = fCobrems->Polarization(x, theta2, phi);
         polarization_phi = M_PI / 2;
         break;
      }
      else {                                     // try incoherent generation
         ++fIncoherentPDFlogx.Ntested;

         double ux = G4UniformRand();
         int i = fIncoherentPDFlogx.search(ux);
         double ui = fIncoherentPDFlogx.integral[i];
         double u_i = (i > 0)? fIncoherentPDFlogx.integral[i-1] : 0;
         double logxi = fIncoherentPDFlogx.randvar[i];
         double dlogx = (i > 0)? logxi - fIncoherentPDFlogx.randvar[i-1]:
                                 fIncoherentPDFlogx.randvar[i+1] - logxi;
         double logx = logxi + dlogx * (0.5 - (ui - ux) / (ui - u_i));
         assert (logx < 0);
         x = exp(logx);
         double dNidxdyPDF = (ui - u_i) / x / dlogx;
         double uy = G4UniformRand();
         int j = fIncoherentPDFy.search(uy);
         double uj = fIncoherentPDFy.integral[j];
         double u_j = (j > 0)? fIncoherentPDFy.integral[j-1] : 0;
         double yj = fIncoherentPDFy.randvar[j];
         double dy = (j > 0)? yj - fIncoherentPDFy.randvar[j-1]:
                             fIncoherentPDFy.randvar[j+1] - yj;
         double y = yj + dy * (0.5 - (uj - uy) / (uj - u_j));
         assert (y > 0);
         dNidxdyPDF *= (uj - u_j) / dy;
         theta2 = fIncoherentPDFtheta02 * (1 / (y + 1e-99) - 1);
         double dNidxdy = fCobrems->Rate_dNidxdt2(x, theta2) *
                          fIncoherentPDFtheta02 / (y*y + 1e-99);
         double Pfactor = dNidxdy / dNidxdyPDF;
         if (Pfactor > fIncoherentPDFlogx.Pmax)
            fIncoherentPDFlogx.Pmax = Pfactor;
         if (Pfactor > fIncoherentPDFy.Pcut) {
            G4cout << "Warning in GenerateBeamPhoton - Pfactor " << Pfactor
                   << " exceeds fIncoherentPDFy.Pcut = " 
                   << fIncoherentPDFy.Pcut
                   << G4endl
                   << "  present x = " << x << G4endl
                   << "  present y = " << y << G4endl
                   << "  present maximum Pfactor = "
                   << fIncoherentPDFlogx.Pmax << G4endl
                   << "  current generator efficiency = "
                   << fIncoherentPDFlogx.Npassed /
                      (fIncoherentPDFlogx.Ntested + 1e-99)
                   << G4endl;
         }
         fIncoherentPDFlogx.Psum += Pfactor;
         if (G4UniformRand() * fIncoherentPDFy.Pcut > Pfactor) {
            continue;
         }
         ++fIncoherentPDFlogx.Npassed;

         phi = 2*M_PI * G4UniformRand();
         polarization = fCobrems->AbremsPolarization(x, theta2, phi);
         polarization_phi = phi - M_PI / 2;
         break;
      }
   }

#if VERBOSE_COBREMS_SPLITTING
   if (fIncoherentPDFlogx.Npassed / 100 * 100 == fIncoherentPDFlogx.Npassed) {
      G4cout << "coherent rate is "
             << fCoherentPDFx.Psum / (fCoherentPDFx.Ntested + 1e-99)
             << ", efficiency is "
             << fCoherentPDFx.Npassed / (fCoherentPDFx.Ntested + 1e-99)
             << G4endl
             << "incoherent rate is "
             << fIncoherentPDFlogx.Psum / (fIncoherentPDFlogx.Ntested + 1e-99)
             << ", efficiency is "
             << fIncoherentPDFlogx.Npassed / (fIncoherentPDFlogx.Ntested + 1e-99)
             << G4endl
             << "counts are "
             << fCoherentPDFx.Npassed << " / " << fIncoherentPDFlogx.Npassed
             << " = "
             << fCoherentPDFx.Npassed / (fIncoherentPDFlogx.Npassed + 1e-99)
             << G4endl;
   }
#endif

   // Put the radiator back the way you found it
   fCobrems->setTargetOrientation(targetThetax/radian,
                                  targetThetay/radian,
                                  targetThetaz/radian);

   // Define the particle kinematics and polarization in lab coordinates
   G4ParticleDefinition *part = GlueXPrimaryGeneratorAction::GetParticle("gamma");
   double Emax = fCobrems->getBeamEnergy() * GeV;
   double Erms = fCobrems->getBeamErms() * GeV;
   double Ebeam = Emax + Erms * G4RandGauss::shoot();
   double theta = sqrt(theta2) * electron_mass_c2 / Emax;
   double alphax = thxBeam + thxMS + theta * cos(phi);
   double alphay = thyBeam + thyMS + theta * sin(phi);
   double pabs = Ebeam * x;
   double px = pabs * alphax;
   double py = pabs * alphay;
   double pz = sqrt(pabs*pabs - px*px - py*py);
   double colphi = 2*M_PI * G4UniformRand();
   double vspotrms = fCobrems->getCollimatorSpotrms() * m;
   double colrho = vspotrms * sqrt(-2 * log(G4UniformRand()));
   double colDist = fCobrems->getCollimatorDistance() * m;
   double radx = colrho * cos(colphi) - colDist * thxBeam;
   double rady = colrho * sin(colphi) - colDist * thyBeam;
   double colx = radx + colDist * alphax;
   double coly = rady + colDist * alphay;
#if defined BEAM_BOX_SIZE
   colx += BEAM_BOX_SIZE * (G4UniformRand() - 0.5);
   coly += BEAM_BOX_SIZE * (G4UniformRand() - 0.5);
#endif
   colx += fBeamOffset[0];
   coly += fBeamOffset[1];
   G4ThreeVector vtx(colx, coly, fBeamStartZ);
   if (fForceFixedPolarization) {
      polarization = fFixedPolarization;
      polarization_phi = fFixedPolarization_phi;
   }
   G4ThreeVector pol(polarization * cos(polarization_phi),
                     polarization * sin(polarization_phi),
                     -(px * polarization * cos(polarization_phi) +
                       py * polarization * sin(polarization_phi)) / pz);
   // use upper half-space to define the polarization plane
   if (pol[2] < 0)
      pol = -pol;
   G4ThreeVector mom(px, py, pz);

   // If beam photon is primary particle, use it to initialize event info
   double targetCenterZ = GlueXPrimaryGeneratorAction::getTargetCenterZ();
   int bg = 1;
   double tvtx;
   if (t0 == 0) {
      tvtx = (vtx[2] - targetCenterZ) / fBeamVelocity;
      tvtx -= GenerateTriggerTime(anEvent);
      if (fGenerateNotSimulate == 0) {
         event_info->AddBeamParticle(1, tvtx, vtx, mom, pol);
      }
      bg = 0;
   }
   else {
      tvtx = fBeamBucketPeriod * floor(t0 / fBeamBucketPeriod + 0.5);
      tvtx += (vtx[2] - targetCenterZ) / fBeamVelocity;
      if (fBeamBackgroundTagOnly) {
         double ttag = tvtx + (targetCenterZ - vtx[2]) / fBeamVelocity;
         GenerateTaggerHit(anEvent, pabs, ttag, bg);
         return;
      }
   }

   // Generate new primary for the beam photon
   G4PrimaryVertex* vertex = new G4PrimaryVertex(vtx, tvtx);
   G4PrimaryParticle* photon = new G4PrimaryParticle(part, px, py, pz);
   photon->SetPolarization(pol);
   vertex->SetPrimary(photon);
   if (fGenerateNotSimulate < 1) {
      anEvent->AddPrimaryVertex(vertex);
   }

   // Include information about the radiating electron, but do not track it
   G4ParticleDefinition *ve = GlueXPrimaryGeneratorAction::GetParticle("e-");
   double vpx = Ebeam * thxBeam;
   double vpy = Ebeam * thyBeam;
   double vpz = sqrt(Ebeam*Ebeam - vpx*vpx - vpy*vpy);
   G4PrimaryParticle *velectron = new G4PrimaryParticle(ve, vpx, vpy, vpz);
   double colvx = radx + colDist * thxBeam;
   double colvy = rady + colDist * thyBeam;
   colvx += fBeamOffset[0];
   colvy += fBeamOffset[1];
   G4ThreeVector vvtx(colvx, colvy, fBeamStartZ);
   G4PrimaryVertex *vvertex = new G4PrimaryVertex(vvtx, tvtx);
   vvertex->SetPrimary(velectron);
   event_info->AddPrimaryVertex(*vvertex);
   delete vvertex;

   // If running in event generation only mode, default is not to
   // save the event to the output file. This will be set back to
   // true if the beam particle makes it to the reference plane.
   if (fGenerateNotSimulate == -1) {
      event_info->SetKeepEvent(0);
   }

   // If bg beam particle, append to MC record
   if (bg || fGenerateNotSimulate == 1) {
      event_info->AddPrimaryVertex(*vertex);
   }

   if (fGenerateNotSimulate == 0) {

      // Register a tagger hit for each beam photon

      double ttag = tvtx + (targetCenterZ - vtx[2]) / fBeamVelocity;
      GenerateTaggerHit(anEvent, pabs, ttag, bg);
   }
}

double GlueXPhotonBeamGenerator::GenerateTriggerTime(const G4Event *event)
{
   // The primary interaction vertex time is referenced to a clock
   // whose t=0 is synchronized to the crossing of a beam bunch
   // through the clock reference plane. This beam bunch may not
   // contain the beam particle whose interaction generated the
   // vertex, but it represents best-guess based on the arrival time
   // of the L1 trigger signal. The spread in the L1 relative to the
   // interacting bunch time is parameterized as a Gaussian.

   static int last_run_number = 0;
   int run_number = HddmOutput::getRunNo();
   if (run_number != last_run_number) {
      fBeamBucketPeriod = getBeamBucketPeriod(run_number);
      double refZ = getRFreferencePlaneZ(run_number);
      GlueXPrimaryGeneratorAction::setRFreferencePlaneZ(refZ);
      last_run_number = run_number;
   }
   double L1sigmat = GlueXPrimaryGeneratorAction::getL1triggerTimeSigma();
   double t0 = L1sigmat * G4RandGauss::shoot();
   return fBeamBucketPeriod * floor(t0 / fBeamBucketPeriod + 0.5);
}

int GlueXPhotonBeamGenerator::GenerateTaggerHit(const G4Event *event,
                                                double energy, 
                                                double time,
                                                int bg)
{
   // Register a tagger hit

   int runNo = HddmOutput::getRunNo();
   if (fTagger == 0) {
      fTagger = new GlueXPseudoDetectorTAG(runNo);
   }
   else if (fTagger->getRunNo() != runNo) {
      delete fTagger;
      fTagger = new GlueXPseudoDetectorTAG(runNo);
   }
   return fTagger->addTaggerPhoton(event, energy, time, bg);
}

void GlueXPhotonBeamGenerator::GenerateRFsync(const G4Event *event)
{
   // Append a RF sync element to the output record, assuming that
   // GenerateTriggerTime has already been called for this event.

   int runNo = HddmOutput::getRunNo();
   if (fTagger == 0) {
      fTagger = new GlueXPseudoDetectorTAG(runNo);
   }
   else if (fTagger->getRunNo() != runNo) {
      delete fTagger;
      fTagger = new GlueXPseudoDetectorTAG(runNo);
   }
   double tsync = 512*ns * G4UniformRand();
   tsync = fBeamBucketPeriod * floor(tsync / fBeamBucketPeriod + 0.5);
   fTagger->addRFsync(event, tsync);
}

double GlueXPhotonBeamGenerator::getBeamBucketPeriod(int runno)
{
   // Look up the beam bucket period for this run in ccdb
   // unless the user has already set the value by hand.

   if (runno > 0) {
      jana::JCalibration *jcalib = japp->GetJCalibration(runno);
      G4cout << "JCalibration context: " << jcalib->GetContext()
             << G4endl;
      std::map<std::string, double> result;
      std::string map_key("/PHOTON_BEAM/RF/beam_period");
      if (jcalib->Get(map_key, result)) {
         G4cerr << "Error in GlueXPhotonBeamGenerator::getBeamBucketPeriod"
                << " - error fetching " << map_key << " from ccdb, "
                << "cannot continue." << G4endl;
         exit(-1);
      }
      else if (result.find("beam_period") != result.end()) {
         fBeamBucketPeriod = result["beam_period"] * ns;
      }
      else {
         G4cerr << "Error in GlueXPhotonBeamGenerator::getBeamBucketPeriod"
                << " - error finding value for " << map_key
                << " in ccdb, cannot continue." << G4endl;
         exit(-1);
      }
   }
   return fBeamBucketPeriod;
}

double GlueXPhotonBeamGenerator::getRFreferencePlaneZ(int runno)
{
   // Look up the reference plane Z for this run in ccdb
   // unless the user has already set the value by hand.

   double refZ = 65 * cm;

   if (runno > 0) {
      jana::JCalibration *jcalib = japp->GetJCalibration(runno);
      G4cout << "JCalibration context: " << jcalib->GetContext()
             << G4endl;
      std::map<std::string, double> result;
      std::string map_key("/PHOTON_BEAM/RF/reference_plane_z");
      if (jcalib->Get(map_key, result)) {
         G4cerr << "Error in GlueXPhotonBeamGenerator::getRFreferencePlaneZ"
                << " - error fetching " << map_key << " from ccdb, "
                << "keeping default value " << refZ / cm << " cm." << G4endl;
      }
      else if (result.find("z_position") != result.end()) {
         refZ = result["z_position"] * cm;
         G4cout << "Info: RF reference plane set to " << refZ / cm << " cm."
                << G4endl;
      }
      else {
         G4cerr << "Error in GlueXPhotonBeamGenerator::getRFreferencePlaneZ"
                << " - error finding value for " << map_key
                << " in ccdb, cannot continue." << G4endl;
         exit(-1);
      }
   }
   return refZ;
}
