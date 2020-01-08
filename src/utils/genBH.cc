//
// genBH.cc : Bethe Heitler event generator
//
// author: richard.t.jones at uconn.edu
// version: march 15, 2018
//
// usage: genBH -n <#> [options] <output_file.hddm>
//   where options may include any of the following
//      -t <#> : number of threads to run, default 1
//      -E <val> : energy of incident photon (GeV), default 9.0
//      -e <val> : use bremsstrahlung spectrum with endpoint e (GeV)
//
// notes:
//  1) If option -e is specified, the value of the E option is taken
//     as the minimum photon energy in the generated bremsstrahlung 
//     spectrum.
//
//  2) The event count option -n <#> must be present, but it does not
//     have to be first among options, as shown in the usage pattern.
//
// physics:
//  1) The generator is based on the tree-level e+e- pair production
//     process from elastic scattering from a free nucleon.
//  2) A complete treatment of the tree-level cross section is used,
//     taking into account polarization of the incident photon and
//     nucleon. It also supports polarization selection of the outgoing
//     nucleon and leptons. Default behavior is to sum over final spins
//     but there is no need to to this if specific polarization in the
//     final state particles is of interest.
//  3) The nucleon target structure is represented by the Dirac and
//     Pauli form factors. The case of scattering from free electrons
//     or coherent scattering from atoms is covered by other tools in
//     the Dirac++ package (see Triplets.C, or or Pairs.C).
//  4) The weight factors written into the output hddm <reaction> tag
//     must be carried with the event. They have units of microbarns,
//     defined such that the total Bethe Heitler cross section can be
//     estimated by Monte Carlo integration as the average value of
//     these weights.
//

#include <HDDM/hddm_mc_s.hpp>

#include <CobremsGeneration.hh>
#include <TPhoton.h>
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

double Ephoton(9.0);
double Endpoint(0.);
std::atomic<int> eventNo(0);

typename hddm_mc_s::ostream *esink;

void usage()
{
   std::cout <<
      "Usage: genBH -n <#> [options] <output_file.hddm>\n"
      "  where options may include any of the following\n"
      "     -t <#> : number of threads to run, default 1\n"
      "     -E <val> : energy of incident photon (GeV), default 9.0\n"
      "     -e <val> : use bremsstrahlung spectrum with endpoint e (GeV)\n"
      << std::endl;
   exit(1);
}

int generate(int nevents)
{
   thread_local TRandom2 *randoms;
   randoms = new TRandom2(0);

   double sum(0);
   double sum2(0);
   for (int n = 0; n < nevents; ++n) {
      LDouble_t weight = 1;
      LDouble_t kin = Ephoton;

      // generate Epos uniform on [0,kin]
      LDouble_t Epos = randoms->Uniform(kin);
      weight *= kin;
   
      // generate phi12 uniform on [0,2pi]
      LDouble_t phi12 = randoms->Uniform(2*PI_);
      weight *= 2*PI_;

      // generate phiR uniform on [0,2pi]
      LDouble_t phiR = randoms->Uniform(2*PI_);
      weight *= 2*PI_;
   
      // generate Mpair with weight (1/M) / (Mcut^2 + M^2)
      LDouble_t Mmin = 2 * mElectron;
      LDouble_t Mcut = 5e-3;                  // 5 MeV cutoff parameter
      LDouble_t um0 = 1 + sqr(Mcut / Mmin);
      LDouble_t um = pow(um0, randoms->Uniform(1));
      LDouble_t Mpair = Mcut / sqrt(um - 1);
      weight *= Mpair * (sqr(Mcut) + sqr(Mpair))
                      * log(um0) / (2 * sqr(Mcut));

      // generate qR^2 with weight (1/qR^2) / sqrt(qRcut^2 + qR^2)
      LDouble_t qRmin = sqr(Mpair )/ (2 * kin);
      LDouble_t qRcut = 1e-3;                 // 1 MeV/c cutoff parameter
      LDouble_t uq0 = qRmin / (qRcut + sqrt(sqr(qRcut) + sqr(qRmin)));
      LDouble_t uq = pow(uq0, randoms->Uniform(1));
      LDouble_t qR2 = sqr(2 * qRcut * uq/(1 - sqr(uq)));
      weight *= qR2 * sqrt(1 + qR2 / sqr(qRcut))
                      * (-2 * log(uq0));

      // overall measure Jacobian factor
      weight *= Mpair / (2 * kin);

      // compute the differential cross section
      TPhoton gIn;
      TLepton eOut(mElectron), pOut(mElectron);
      TLepton nIn(mProton), nOut(mProton);
      LDouble_t diffXS(0);

      try {
         nIn.SetMom(TThreeVectorReal(0,0,0));
         gIn.SetMom(TThreeVectorReal(0,0,kin));

         // Solve for the rest of the kinematics
         LDouble_t qR = sqrt(qR2);
         LDouble_t Erec = sqrt(sqr(mProton) + qR2);
         //LDouble_t Eele = kin + mProton - Erec - Epos;
         LDouble_t costhetaR = (2 * (mProton + kin) * (Erec - mProton) +
                                sqr(Mpair)) / (2 * kin * qR);
         if (costhetaR > 1) {
            throw std::runtime_error("no kinematic solution because costhetaR > 1");
         }
         LDouble_t sinthetaR = sqrt(1 - sqr(costhetaR));
         TThreeVectorReal qRecoil(qR * sinthetaR * cos(phiR),
                                  qR * sinthetaR * sin(phiR),
                                  qR * costhetaR);
         nOut.SetMom(qRecoil);

         LDouble_t pStar2 = sqr(Mpair / 2) - sqr(mElectron);
         if (pStar2 < 0) {
            throw std::runtime_error("no kinematic solution because pStar2 < 0");
         }
         LDouble_t pStar = sqrt(pStar2);
         LDouble_t E12 = kin + mProton - Erec;
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
         TLorentzBoost toLab(qRecoil[1] / E12,
                             qRecoil[2] / E12,
                             (qRecoil[3] - kin) / E12);
         p1.Boost(toLab);
         pOut.SetMom(p1);
         TThreeVectorReal p2(gIn.Mom() - qRecoil - p1);
         eOut.SetMom(p2);

         // Set the initial,final polarizations
         gIn.SetPol(TThreeVectorReal(1,0,0));
         nIn.SetPol(TThreeVectorReal(0,0,0));
         eOut.AllPol();
         pOut.AllPol();
         nOut.AllPol();

         // Compute the polarized differential cross section
         LDouble_t F1_spacelike = 1;
         LDouble_t F2_spacelike = 0;
         LDouble_t F1_timelike = 1;
         LDouble_t F2_timelike = 0;
         diffXS = TCrossSection::BetheHeitlerNucleon(gIn, nIn,
                                                     eOut, pOut, nOut,
                                                     F1_spacelike,
                                                     F2_spacelike,
                                                     F1_timelike,
                                                     F2_timelike);
      }
      catch (const std::exception &e) {

         // These events have no cross section, but do not discard
         // them because they are needed to get the right MC integral.
 
         nOut.SetMom(TThreeVectorReal(0,0,0));
         eOut.SetMom(TThreeVectorReal(0,0,0));
         pOut.SetMom(TThreeVectorReal(0,0,0));
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
      hddm_mc_s::ReactionList reactions = events(0).addReactions();
      reactions(0).setType(235586);
      reactions(0).setWeight(weight * diffXS);
      hddm_mc_s::BeamList beams = reactions(0).addBeams();
      hddm_mc_s::MomentumList bmoms = beams(0).addMomenta();
      hddm_mc_s::PropertiesList bprops = beams(0).addPropertiesList();
      beams(0).setType(Gamma);
      bmoms(0).setE(kin);
      bmoms(0).setPx(0);
      bmoms(0).setPy(0);
      bmoms(0).setPz(kin);
      bprops(0).setCharge(0);
      bprops(0).setMass(0);
      hddm_mc_s::TargetList targs = reactions(0).addTargets();
      hddm_mc_s::MomentumList tmoms = targs(0).addMomenta();
      hddm_mc_s::PropertiesList tprops = targs(0).addPropertiesList();
      targs(0).setType(Proton);
      tmoms(0).setE(mProton);
      tmoms(0).setPx(0);
      tmoms(0).setPy(0);
      tmoms(0).setPz(0);
      tprops(0).setCharge(1);
      tprops(0).setMass(mProton);
      hddm_mc_s::VertexList verts = reactions(0).addVertices();
      hddm_mc_s::ProductList prods = verts(0).addProducts(3);
      prods(0).setDecayVertex(0);
      prods(0).setId(1);
      prods(0).setMech(0);
      prods(0).setParentid(0);
      prods(0).setPdgtype(-11);
      prods(0).setType(Positron);
      prods(0).setDecayVertex(0);
      prods(1).setId(2);
      prods(1).setMech(0);
      prods(1).setParentid(0);
      prods(1).setPdgtype(11);
      prods(1).setType(Electron);
      prods(1).setDecayVertex(0);
      prods(2).setId(3);
      prods(2).setMech(0);
      prods(2).setParentid(0);
      prods(2).setPdgtype(2212);
      prods(2).setType(Proton);
      for (int i=0; i < 3; ++i) {
         prods(i).addMomenta();
         prods(i).addPropertiesList();
      }
      prods(0).getMomentum().setE(pOut.Mom()[0]);
      prods(0).getMomentum().setPx(pOut.Mom()[1]);
      prods(0).getMomentum().setPy(pOut.Mom()[2]);
      prods(0).getMomentum().setPz(pOut.Mom()[3]);
      prods(0).getProperties().setMass(mElectron);
      prods(0).getProperties().setCharge(1);
      prods(1).getMomentum().setE(eOut.Mom()[0]);
      prods(1).getMomentum().setPx(eOut.Mom()[1]);
      prods(1).getMomentum().setPy(eOut.Mom()[2]);
      prods(1).getMomentum().setPz(eOut.Mom()[3]);
      prods(1).getProperties().setMass(mElectron);
      prods(1).getProperties().setCharge(-1);
      prods(2).getMomentum().setE(nOut.Mom()[0]);
      prods(2).getMomentum().setPx(nOut.Mom()[1]);
      prods(2).getMomentum().setPy(nOut.Mom()[2]);
      prods(2).getMomentum().setPz(nOut.Mom()[3]);
      prods(2).getProperties().setMass(mProton);
      prods(2).getProperties().setCharge(+1);
      verts(0).addOrigins();
      verts(0).getOrigin().setT(0);
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
            Ephoton = std::atof(arg+2);
         else
            Ephoton = atof(argv[++iarg]);
      }
      else if (arg[1] == 'e') {
         if (strlen(arg) > 2)
            Endpoint = std::atof(arg+2);
         else
            Endpoint = atof(argv[++iarg]);
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
