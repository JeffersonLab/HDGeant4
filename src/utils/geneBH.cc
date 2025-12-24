//
// geneBH.cc : electron Bethe Heitler event generator
//
// author: richard.t.jones at uconn.edu
// version: july 5, 2021
//
// usage: genBH -n <#> [options] <output_file.hddm>
//   where options may include any of the following
//      -t <#> : number of threads to run, default 1
//      -E <val> : energy of incident electron (GeV), default 11.0 GeV
//      -m <val> : generate e+e- pairs with invariant mass >= val (GeV/c^2)
//      -r <val> : set initial random number seed to val
//      -R <run> : set simulation run number, default 9000
//
// notes:
//  1) If option -m is not given then the full invariant mass spectrum
//     of the pairs is generated, from 2*mElectron to the rootS-mTarget.
//
//  2) The event count option -n <#> must be present, but it does not
//     have to be first among options, as shown in the usage pattern.
//
// physics:
//  1) The generator is based on the tree-level e+e- pair production
//     process from elastic scattering from an atomic target.
//  2) A complete treatment of the tree-level cross section is used,
//     taking into account polarization of the incident electron and
//     the outgoing leptons. Default behavior is to average over the
//     iniital spin and sum over final spins, but there is no need to
//     do this if specific polarization in the beam or final state
//     particles is of interest.
//  3) The target atom is treated as a zero spin, with a parameterized
//     nuclear charge form factor. Screening by atomic electrons is
//     described by an atomic form factor.
//  4) The weight factors written into the output hddm <reaction> tag
//     must be carried with the event. They have units of microbarns,
//     defined such that the total Bethe Heitler cross section can be
//     estimated by Monte Carlo integration as the average value of
//     these weights. To get the average right, do not drop events
//     from the output that have zero weight, as these contribute to
//     the denominator of the average, although not the numerator.
//  5) Support has been added for pair production from electrons in
//     the target, in addition to coherent scattering from the atom
//     as a whole which is dominated by the nucleus. This process is
//     sometimes called triplet production rather than pair because
//     there are 3 leptons that emerge from the target along with
//     the incident electron.

#include <HDDM/hddm_mc_s.hpp>

#include <TLepton.h>
#include <TLorentzBoost.h>
#include <TCrossSection.h>
#include <constants.h>
#include <sqr.h>

#include <TRandom2.h>

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include <thread>
#include <stdexcept>

#define TRIPLET_FRACTION 0.1

double Ebeam(3.4);
LDouble_t M12_minimum(0);
int seedVal(0);
uint runno(9000);

LDouble_t mLepton(mElectron);
Particle_t pLepton(Positron);
Particle_t nLepton(Electron);
int nLeptonPDG(11);
int pLeptonPDG(-11);

// Define the properties of the target
Particle_t Target(Ta181);
int TargetPDG(1000721810);
LDouble_t mTarget(180.94788 * 0.93149410242);
LDouble_t TargetZ(73);

LDouble_t FFatomic(LDouble_t qrec) {
   // A crude parameterization of the electron cloud form factor
   // for an atom with atomic number TargetZ, taken from ref.
   // G4BetheHeitler5DModel.cc (Bernard).
   const LDouble_t beta = 2.17e5 / pow(TargetZ, 1/3.);
   return 1 / (1 + sqr(beta * qrec));
}
LDouble_t FFnuclear(LDouble_t qrec) {
   // This needs to be modified to describe lepton pairs
   // with masses above 100 MeV/c^2.
   return 1;
}

std::atomic<int> eventNo(0);

typename hddm_mc_s::ostream *esink;

void usage()
{
   std::cout <<
      "Usage: geneBH -n <#> [options] <output_file.hddm>\n"
      "  where options may include any of the following\n"
      "     -t <#> : number of threads to run, default 1\n"
      "     -E <val> : energy of incident electron (GeV), default 11.0\n"
      "     -m <val> : generate e+e- pairs with invariant mass >= val (GeV/c^2)\n"
      "     -r <val> : set initial random number seed to val\n"
      "     -R <run> : set simulation run number, default 9000\n"
      << std::endl;
   exit(1);
}

int generate(int nevents)
{
   LDouble_t pbeam = sqrt(sqr(Ebeam) - sqr(mElectron));

   thread_local TRandom2 *randoms;
   randoms = new TRandom2(0);
   randoms->SetSeed(seedVal);

   double sum(0);
   double sum2(0);
   for (int n = 0; n < nevents; ++n) {
      LDouble_t weight = 1;

      // generate Mpair with weight (1/M) / (Mcut^2 + M^2)
      LDouble_t Mthresh = 2 * mLepton;
      LDouble_t Mcut = 10 * mLepton;      // 5 MeV cutoff parameter for e+e-
      LDouble_t um0 = 1 + sqr(Mcut / Mthresh);
      LDouble_t umax = 1;
      if (M12_minimum > Mthresh)
         umax = log(1 + sqr(Mcut / M12_minimum)) / log(um0);
      LDouble_t um = pow(um0, randoms->Uniform(umax));
      LDouble_t Mpair = Mcut / sqrt(um - 1);
      weight *= Mpair * (sqr(Mcut) + sqr(Mpair))
                      * log(um0) / (2 * sqr(Mcut));
      weight *= umax;

      // generate nu with uniform weight on [Mpair, Ebeam - mLepton]
      LDouble_t nu = Mpair + randoms->Uniform(Ebeam - Mpair - mLepton);
      weight *= Ebeam - Mpair - mLepton;

      // generate Q^2 with weight (1/Q^2) / sqrt(Qcut^2 + Q^2)
      LDouble_t qmin = sqrt(sqr(Ebeam) - sqr(mLepton)) -
                       sqrt(sqr(Ebeam - nu) - sqr(mLepton));
      LDouble_t Qmin = sqrt(sqr(qmin) - sqr(nu));
      LDouble_t Qcut = 2 * mLepton;      // 1 MeV/c cutoff parameter for e+e-
      LDouble_t uQ0 = Qmin / (Qcut + sqrt(sqr(Qcut) + sqr(Qmin)));
      LDouble_t uQ = pow(uQ0, randoms->Uniform(1));
      LDouble_t Q2 = sqr(2 * Qcut * uQ/(1 - sqr(uQ)));
      weight *= Q2 * sqrt(1 + Q2 / sqr(Qcut))
                      * (-2 * log(uQ0));
 
      // generate qR^2 with weight (1/qR^2) / sqrt(qRcut^2 + qR^2)
      LDouble_t qRmin = sqrt(sqr(nu) + Q2) - sqrt(sqr(nu) - sqr(Mpair));
      LDouble_t qRcut = 2 * mLepton;      // 1 MeV/c cutoff parameter for e+e-
      LDouble_t uq0 = qRmin / (qRcut + sqrt(sqr(qRcut) + sqr(qRmin)));
      LDouble_t uq = pow(uq0, randoms->Uniform(1));
      LDouble_t qR2 = sqr(2 * qRcut * uq/(1 - sqr(uq)));
      weight *= qR2 * sqrt(1 + qR2 / sqr(qRcut))
                      * (-2 * log(uq0));

      // generate Epos uniform on [0,nu]
      LDouble_t Epos = randoms->Uniform(nu);
      weight *= nu;
   
      // generate phi12 uniform on [0,2pi]
      LDouble_t phi12 = randoms->Uniform(2*PI_);
      weight *= 2*PI_;

      // generate phiR uniform on [0,2pi]
      LDouble_t phiR = randoms->Uniform(2*PI_);
      weight *= 2*PI_;
 
      // generate phiq uniform on [0,2pi]
      LDouble_t phiq = randoms->Uniform(2*PI_);
      weight *= 2*PI_;

      // fork between coherent pair and triplet process
      LDouble_t target_mass;
      if (randoms->Uniform(1) < TRIPLET_FRACTION) {
         weight /= TRIPLET_FRACTION;
         target_mass = mElectron;
      }
      else {
         weight /= 1 - TRIPLET_FRACTION + 1e-99;
         target_mass = mTarget;
      }
   
      // overall measure Jacobian factors
      weight *= Mpair / (2 * sqrt(sqr(nu) + Q2));
      weight *= (Ebeam - nu) / (2 * pbeam);

      // compute the differential cross section
      TLepton eIn(mElectron), eOut(mElectron);
      TLepton lnOut(mLepton), lpOut(mLepton);
      TLepton tIn(mElectron), tOut(mElectron);
      TFourVectorReal qRecoil;
      LDouble_t diffXS(0);

      try {
         TFourVectorReal pIn(Ebeam,0,0,pbeam);

         // Solve for the rest of the kinematics:
         //   q is referenced to the beam direction (0,0,1)
         //   qR is referenced to the direction of q
         //   p+ is referenced to the direction of q12=q-qR

         LDouble_t qin = sqrt(sqr(nu) + Q2);
         LDouble_t costhetaq = (Q2 + 2 * Ebeam * nu) / (2 * pbeam * qin);
         if (costhetaq > 1) {
            throw std::runtime_error("no kinematic solution because costhetaq > 1");
         }
         LDouble_t sinthetaq = sqrt(1 - sqr(costhetaq));
         TFourVectorReal q(nu, qin * sinthetaq * cos(phiq),
                               qin * sinthetaq * sin(phiq),
                               qin * costhetaq);

         LDouble_t qR = sqrt(qR2);
         LDouble_t Erec = sqrt(sqr(target_mass) + qR2);
         LDouble_t costhetaqR = (2 * (target_mass + nu) * (Erec - target_mass) +
                                sqr(Mpair) + Q2) / (2 * qin * qR);
         if (costhetaqR > 1) {
            throw std::runtime_error("no kinematic solution because costhetaR > 1");
         }
         LDouble_t sinthetaqR = sqrt(1 - sqr(costhetaqR));
         qRecoil = TFourVectorReal(Erec - target_mass,
                                   qR * sinthetaqR * cos(phiR),
                                   qR * sinthetaqR * sin(phiR),
                                   qR * costhetaqR);
         qRecoil.Rotate(0, -asin(sinthetaq), -phiq);

         LDouble_t pStar2 = sqr(Mpair / 2) - sqr(mLepton);
         if (pStar2 < 0) {
            throw std::runtime_error("no kinematic solution because pStar2 < 0");
         }
         LDouble_t pStar = sqrt(pStar2);
         LDouble_t E12 = nu + target_mass - Erec;
         LDouble_t p12mag = sqrt(sqr(E12) - sqr(Mpair));
         LDouble_t costhetastar = (Epos - E12 / 2) * Mpair / (pStar * p12mag);
         if (fabs(costhetastar) > 1) {
            throw std::runtime_error("no kinematic solution because costhetastar < 1");
         }
         LDouble_t sinthetastar = sqrt(1 - sqr(costhetastar));
         TThreeVectorReal k12(pStar * sinthetastar * cos(phi12),
                              pStar * sinthetastar * sin(phi12),
                              pStar * costhetastar);
         TFourVectorReal p1(Mpair / 2, k12);
         TFourVectorReal ppair(q - qRecoil);
         if (fabs(ppair[0] - E12) > E12 * 1e-12) {
            throw std::runtime_error("no kinematic solution because ppair[0] != E12");
         }
         else if (fabs(ppair.Length() - p12mag) > p12mag * 1e-12) {
            throw std::runtime_error("no kinematic solution because |ppair| != |p12|");
         }
         TLorentzBoost toLab(-ppair[1] / ppair[0],
                             -ppair[2] / ppair[0],
                             -ppair[3] / ppair[0]);
         p1.Boost(toLab);

         eIn.SetMom(pIn);
         eOut.SetMom(pIn - q);
         lpOut.SetMom(p1);
         lnOut.SetMom(ppair - p1);

         // Set the initial,final polarizations
         eIn.SetPol(TThreeVectorReal(0,0,0));
         eOut.AllPol();
         lpOut.AllPol();
         lnOut.AllPol();

         if (target_mass == mElectron) {
            tIn.SetMom(TThreeVectorReal(0,0,0));
            tIn.SetPol(TThreeVectorReal(0,0,0));
            tOut.SetMom((TThreeVectorReal)qRecoil);
            tOut.AllPol();
            diffXS = TCrossSection::eTripletProduction(eIn, eOut, lpOut, lnOut, tIn, tOut);
            diffXS *= TargetZ;

            // Include the atomic form factor
            diffXS *= 1 - sqr(FFatomic(qR));
 
            // Suppress over-counting due to identical electrons in final state
            if (eOut.Mom()[0] < lnOut.Mom()[0] || lnOut.Mom()[0] < tOut.Mom()[0])
               weight = 0;
         }
         else {
            diffXS = TCrossSection::ePairProduction(eIn, eOut, lpOut, lnOut);
            diffXS *= sqr(TargetZ);

            // Include the atomic and nuclear form factors
            diffXS *= sqr(1 - FFatomic(qR));
            diffXS *= sqr(FFnuclear(qR));
 
            // Suppress double counting of the forward electron, assumes e+e-
            if (eOut.Mom()[0] < lnOut.Mom()[0])
               weight = 0;
         }
      }
      catch (const std::exception &e) {

         // These events have no cross section, but do not discard
         // them because they are needed to get the right MC integral.
 
         eOut.SetMom(TThreeVectorReal(0,0,1e-12));
         lpOut.SetMom(TThreeVectorReal(0,0,1e-12));
         lnOut.SetMom(TThreeVectorReal(0,0,1e-12));
         diffXS = 0;
      }

      // Keep statistics on the total cross section
      sum += weight * diffXS;
      sum2 += sqr(weight * diffXS);
      if (n > 0 && n / 10000 * 10000 == n) {
         std::cout << "est. total cross section after " 
                   << n << " events : " << sum / n 
                   << " +/- " << sqrt(sum2 - sqr(sum) / n) / n
                   << " ub" << std::endl;
      }

      // Save kinematics to output event record
      hddm_mc_s::HDDM record;
      hddm_mc_s::PhysicsEventList events = record.addPhysicsEvents();
      events(0).setEventNo(++eventNo);
      events(0).setRunNo(runno);
      hddm_mc_s::ReactionList reactions = events(0).addReactions();
      reactions(0).setType(235586);
      reactions(0).setWeight(weight * diffXS);
      // Write 4 random seeds to the hddm file.
      // seed1 = saved by generator, saved for future reference, not to be reused
      // seed2 = saved by simulator (hdgeant or hdgeant4), saved for future reference, not to be reused
      // seed3 = saved by mcsmear, not to be reused
      // seed4 = saved by analyzer, not to be reused
      hddm_mc_s::RandomList ranl = reactions(0).addRandoms();
      ranl().setSeed1(randoms->GetSeed());
      ranl().setSeed2(gRandom->Integer(std::numeric_limits<int32_t>::max()));
      ranl().setSeed3(gRandom->Integer(std::numeric_limits<int32_t>::max()));
      ranl().setSeed4(gRandom->Integer(std::numeric_limits<int32_t>::max()));
      hddm_mc_s::BeamList beams = reactions(0).addBeams();
      hddm_mc_s::MomentumList bmoms = beams(0).addMomenta();
      hddm_mc_s::PropertiesList bprops = beams(0).addPropertiesList();
      beams(0).setType(Electron);
      bmoms(0).setE(Ebeam);
      bmoms(0).setPx(0);
      bmoms(0).setPy(0);
      bmoms(0).setPz(pbeam);
      bprops(0).setCharge(-1);
      bprops(0).setMass(mElectron);
      hddm_mc_s::TargetList targs = reactions(0).addTargets();
      hddm_mc_s::MomentumList tmoms = targs(0).addMomenta();
      hddm_mc_s::PropertiesList tprops = targs(0).addPropertiesList();
      targs(0).setType(Target);
      tmoms(0).setE(mTarget);
      tmoms(0).setPx(0);
      tmoms(0).setPy(0);
      tmoms(0).setPz(0);
      tprops(0).setCharge(0);
      tprops(0).setMass(mTarget);
      hddm_mc_s::VertexList verts = reactions(0).addVertices();
      hddm_mc_s::ProductList prods = verts(0).addProducts(4);
      prods(0).setDecayVertex(0);
      prods(0).setId(1);
      prods(0).setMech(0);
      prods(0).setParentid(0);
      prods(0).setPdgtype(pLeptonPDG);
      prods(0).setType(pLepton);
      prods(0).setDecayVertex(0);
      prods(1).setId(2);
      prods(1).setMech(0);
      prods(1).setParentid(0);
      prods(1).setPdgtype(nLeptonPDG);
      prods(1).setType(nLepton);
      prods(1).setDecayVertex(0);
      prods(2).setId(3);
      prods(2).setMech(0);
      prods(2).setParentid(0);
      prods(2).setPdgtype(11);
      prods(2).setType(Electron);
      prods(3).setId(4);
      prods(3).setMech(0);
      prods(3).setParentid(0);
      if (target_mass == mElectron) {
         prods(3).setPdgtype(11);
         prods(3).setType(Electron);
      }
      else {
         prods(3).setPdgtype(TargetPDG);
         prods(3).setType(Target);
      }
      for (int i=0; i < 4; ++i) {
         prods(i).addMomenta();
         prods(0).getMomentum().addMomentum_doubles();
         prods(i).addPropertiesList();
      }
      prods(0).getMomentum().setE(lpOut.Mom()[0]);
      prods(0).getMomentum().setPx(lpOut.Mom()[1]);
      prods(0).getMomentum().setPy(lpOut.Mom()[2]);
      prods(0).getMomentum().setPz(lpOut.Mom()[3]);
      prods(0).getMomentum().getMomentum_double().setE(lpOut.Mom()[0]);
      prods(0).getMomentum().getMomentum_double().setPx(lpOut.Mom()[1]);
      prods(0).getMomentum().getMomentum_double().setPy(lpOut.Mom()[2]);
      prods(0).getMomentum().getMomentum_double().setPz(lpOut.Mom()[3]);
      prods(0).getProperties().setMass(mLepton);
      prods(0).getProperties().setCharge(1);
      prods(1).getMomentum().setE(lnOut.Mom()[0]);
      prods(1).getMomentum().setPx(lnOut.Mom()[1]);
      prods(1).getMomentum().setPy(lnOut.Mom()[2]);
      prods(1).getMomentum().setPz(lnOut.Mom()[3]);
      prods(1).getMomentum().getMomentum_double().setE(lnOut.Mom()[0]);
      prods(1).getMomentum().getMomemtum_double().setPx(lnOut.Mom()[1]);
      prods(1).getMomentum().getMomentum_double().setPy(lnOut.Mom()[2]);
      prods(1).getMomentum().getMomentum_double().setPz(lnOut.Mom()[3]);
      prods(1).getProperties().setMass(mLepton);
      prods(1).getProperties().setCharge(-1);
      prods(2).getMomentum().setE(eOut.Mom()[0]);
      prods(2).getMomentum().setPx(eOut.Mom()[1]);
      prods(2).getMomentum().setPy(eOut.Mom()[2]);
      prods(2).getMomentum().setPz(eOut.Mom()[3]);
      prods(2).getMomentum().getMomentum_double().setE(eOut.Mom()[0]);
      prods(2).getMomentum().getMomentum_double().setPx(eOut.Mom()[1]);
      prods(2).getMomentum().getMomentum_double().setPy(eOut.Mom()[2]);
      prods(2).getMomentum().getMomentum_double().setPz(eOut.Mom()[3]);
      prods(2).getProperties().setMass(mElectron);
      prods(2).getProperties().setCharge(-1);
      prods(3).getMomentum().setE(qRecoil[0]);
      prods(3).getMomentum().setPx(qRecoil[1]);
      prods(3).getMomentum().setPy(qRecoil[2]);
      prods(3).getMomentum().setPz(qRecoil[3]);
      prods(3).getMomentum().getMomentum_double().setE(qRecoil[0]);
      prods(3).getMomentum().getMomentum_double().setPx(qRecoil[1]);
      prods(3).getMomentum().getMomentum_double().setPy(qRecoil[2]);
      prods(3).getMomentum().getMomentum_double().setPz(qRecoil[3]);
      prods(3).getProperties().setMass(target_mass);
      if (target_mass == mElectron) {
         prods(3).getProperties().setCharge(-1);
      }
      else {
         prods(3).getProperties().setCharge(TargetZ);
      }
      verts(0).addOrigins();
      verts(0).getOrigin().setT(Q2);
      verts(0).getOrigin().setVx(0);
      verts(0).getOrigin().setVy(0);
      verts(0).getOrigin().setVz(0);
      *esink << record;
   }

   std::cout << "est. total cross section after " 
             << nevents << " events : " << sum / nevents 
             << " +/- " << sqrt(sum2 - sqr(sum) / nevents) / nevents
             << " ub" << std::endl;

   delete randoms;
   return 0;
}

int main(int argc, char *argv[])
{
   int nevents(0);
   int nthreads(1);

   std::string outfname;
   for (int iarg=1; iarg < argc; ++iarg) {
      const char *arg(argv[iarg]);
      if (arg[0] != '-') {
         outfname.assign(arg);
         break;
      }
      else if (arg[1] == 'n') {
         if (strlen(arg) > 2)
            nevents = std::atoi(arg+2);
         else
            nevents = atoi(argv[++iarg]);
      }
      else if (arg[1] == 't') {
         if (strlen(arg) > 2)
            nthreads = std::atoi(arg+2);
         else
            nthreads = atoi(argv[++iarg]);
      }
      else if (arg[1] == 'E') {
         if (strlen(arg) > 2)
            Ebeam = std::atof(arg+2);
         else
            Ebeam = atof(argv[++iarg]);
      }
      else if (arg[1] == 'm') {
         if (strlen(arg) > 2)
            M12_minimum = std::atof(arg+2);
         else
            M12_minimum = atof(argv[++iarg]);
      }
      else if (arg[1] == 'r') {
         if (strlen(arg) > 2)
            seedVal = std::atoi(arg+2);
         else
            seedVal = atoi(argv[++iarg]);
      }
      else if (arg[1] == 'R') {
         if (strlen(arg) > 2)
            runno = std::atoi(arg+2);
         else
            runno = atoi(argv[++iarg]);
      }
      else {
         usage();
      }
   }
   if (nevents == 0 || outfname.size() == 0)
      usage();

   std::ofstream fout(outfname);
   if (!fout.is_open()) {
      std::cerr << "genBH error - cannot open " 
                << outfname << " for output, giving up."
                << std::endl;
      exit(2);
   }
   esink = new hddm_mc_s::ostream(fout);

   int nevents_per_thread = nevents / nthreads;
   if (nevents_per_thread * nthreads < nevents)
      nevents_per_thread += 1;
   std::vector<std::thread> threads;
   for (int n=0; n < nevents;) {
      int ni(nevents_per_thread);
      ni = (ni > nevents - n)? nevents - n : ni;
      threads.push_back(std::thread(generate,ni));
      n += ni;
   }

   for (int i=0; i < (int)threads.size(); ++i) {
      threads[i].join();
   }
}
