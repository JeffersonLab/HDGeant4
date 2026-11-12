//
// PairConversionGeneration class implementation
//
// author: richard.t.jones at uconn.edu
// version: december 17, 2016
//
// notes:
//
// This class computes differential pair production rates with full
// polarization degrees of freedom taken into account for both the
// photon and the leptons. Both incoherent pair production with a
// recoil electron and nuclear/coherent atomic scattering without a
// free recoil electron are included. Calculations are performed by
// computing directly the Feynman amplitudes for all 8 leading-order
// graphs and summing them together. No approximations are used.
//
// dependencies:
//
// 1. This class depends on the Dirac++ package by the same author.
// 2. The Dirac++ package uses extended precision floating point arithmetic
//    called "long double" in c++. On Intel cpu architectures, this is 
//    implemented using 80-bit floating-point registers known as "x87".
//    Extended precision is needed to keep rounding errors in the sum over
//    cancelling Feynman amplitudes from producing excessive rounding
//    errors in the computation of the differential cross section.
//
// units:
// Any length is in m; energy,momentum,mass in GeV (c=1); angles in
// radians; time in seconds; cross section in barns.

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/python.hpp>

#ifdef USING_DIRACXX
#define BOOST_PYTHON_WRAPPING 1

#include <PairConversionGeneration.hh>

#include "Complex.h"
#include "TCrossSection.h"
#include "TLorentzBoost.h"

#include "G4ios.hh"

#include <iostream>

inline unsigned int sqr(unsigned int x) { return x*x; }
inline Int_t sqr(Int_t x) { return x*x; }
inline Double_t sqr(Double_t x) { return x*x; }
inline LDouble_t sqr(LDouble_t x) { return x*x; }
inline Complex_t sqr(Complex_t x) { return x*x; }

const TThreeVectorReal zeroVector(0,0,0);
const TThreeVectorReal posXhat(1,0,0);
const TThreeVectorReal negXhat(-1,0,0);
const TThreeVectorReal posYhat(0,1,0);
const TThreeVectorReal negYhat(0,-1,0);
const TThreeVectorReal posZhat(0,0,1);
const TThreeVectorReal negZhat(0,0,-1);

PairConversionGeneration::PairConversionGeneration(std::vector<double> Z,
                                                   std::vector<double> /* A */,
                                                   std::vector<double> /* w */)
{
   if (Z.size() == 1) {
       fConverterZ = Z[0];
   }
   else {
      std::cerr << "PairConversionGeneration constructor error: "
                << "converter material is a compound or mixture, " << std::endl
                << "only pure elements are supported by the "
                << "present algorithm, cannot continue." << std::endl;
      exit(13);
   }
}

PairConversionGeneration::~PairConversionGeneration()
{}

LDouble_t PairConversionGeneration::FFatomic(LDouble_t qRecoil)
{
   // return the atomic form factor of the pair converter
   // normalized to unity at zero momentum transfer qRecoil (GeV/c).
   // Lengths are in Angstroms in this function.

   LDouble_t ff(0);
   if (fConverterZ == 1) {
      LDouble_t a0Bohr = 0.529177 / 1.97327e-6;
      ff = 1 / pow(1 + pow(a0Bohr * qRecoil / 2, 2), 2);
   }
   else if (fConverterZ == 4) {
      // parameterization for 4Be given by online database at
      // http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction
      //                    /atomicformfactors/formfactors.php
      LDouble_t acoeff[] = {1.5919, 1.1278, 0.5391, 0.7029};
      LDouble_t bcoeff[] = {43.6427, 1.8623, 103.483, 0.5420};
      LDouble_t ccoeff[] = {0.0385};
      LDouble_t q_invA = qRecoil / 1.97327e-6;
      ff = ccoeff[0];
      for (int i=0; i < 4; ++i) {
         ff += acoeff[i] * exp(-bcoeff[i] * pow(q_invA / (4 * M_PI), 2));
      }
      ff /= fConverterZ;
   }
   else if (fConverterZ > 1 && fConverterZ < 93) {
      // parameterization implemented by Bernard et al
      // in Geant4 class G4BetheHeitler5DModel.
      const LDouble_t beta = 2.17e5 / pow(fConverterZ, 1/3.);
	  ff = 1 / (1 + sqr(beta * qRecoil));
   }
   else {
      std::cerr << "PairConversionGeneration::FFatomic error: "
                << "no model currently implemented for element "
                << "Z=" << fConverterZ << std::endl;
      exit(13);
   }
   return ff;
}

LDouble_t PairConversionGeneration::DiffXS_pair(const TPhoton &gIn, 
                                                const TLepton &pOut,
                                                const TLepton &eOut)
{
   // Calculates the lepton pair production cross section for a
   // gamma ray off an atom at a particular recoil momentum vector q.
   // The cross section is returned as d(sigma)/(dE dphi d^3q) where E is
   // the energy of the final-state lepton and phi is its azimuthal angle.
   // The polar angles of the pair are fixed by momentum conservation.
   // It is assumed that gIn.Mom()[0] = eOut.Mom()[0]+pOut.Mom()[0], that
   // the energy carried away by the recoil is zero in the laboratory frame,
   // but it is not checked.  The calculation is performed in the lab frame.
   // This cross section is only a partial result, because it does not
   // include the integral d^3 q over the form factor of the target.  This
   // depends on the crystal structure of the target atom, and so is left to
   // be carried out by more specialized code.  Units are barns/GeV^4.

   TPhoton g0(gIn);
   TLepton p1(pOut);
   TLepton e2(eOut);

   // Set the initial,final polarizations
   p1.AllPol();
   e2.AllPol();

   // Multiply the basic cross section by the converter atomic form factor
   LDouble_t result = TCrossSection::PairProduction(g0, e2, p1);
   TFourVectorReal qR(gIn.Mom() - eOut.Mom() - pOut.Mom());
   result *= sqr(fConverterZ * (1 - FFatomic(qR.Length())));
   return result * 1e-6;

   // The unpolarized Bethe-Heitler cross section is given here for comparison
   LDouble_t kin = gIn.Mom()[0];
   LDouble_t Epos = pOut.Mom()[0];
   LDouble_t Eneg = eOut.Mom()[0];
   LDouble_t mLepton = eOut.Mass();
   LDouble_t delta = 136 * mLepton / pow(fConverterZ, 0.33333) *
                     kin / (Eneg * Epos);
   LDouble_t aCoul = sqr(alphaQED * fConverterZ);
   LDouble_t fCoul = aCoul * (1 / (1 + aCoul) + 0.20206 - 0.0369 * aCoul +
                              0.0083 * pow(aCoul, 2) - 0.002 * pow(aCoul, 3));
   LDouble_t xsi = log(1440 / pow(fConverterZ, 0.66667)) / 
                   (log(183 / pow(fConverterZ, 0.33333) - fCoul));
   LDouble_t FofZ = (8./3.) * log(fConverterZ) + ((kin < 0.05)? 0 : 8 * fCoul);
   LDouble_t Phi1 = 20.867 - 3.242 * delta + 0.625 * sqr(delta);
   LDouble_t Phi2 = 20.209 - 1.930 * delta - 0.086 * sqr(delta);
   LDouble_t Phi0 = 21.12 - 4.184 * log(delta + 0.952);
   if (delta > 1) {
      Phi1 = Phi2 = Phi0;
   }
   result = hbarcSqr / sqr(mLepton) * pow(alphaQED, 3) / kin
            * fConverterZ * (fConverterZ + xsi)
            * (
                 (sqr(Eneg) + sqr(Epos)) / sqr(kin) * (Phi1 - FofZ/2) +
                 (2./3.) * (Eneg * Epos) / sqr(kin) * (Phi2 - FofZ/2)
              );
   return result * 1e-6;
}

LDouble_t PairConversionGeneration::DiffXS_triplet(const TPhoton &gIn,
                                                   const TLepton &pOut,
                                                   const TLepton &eOut2,
                                                   const TLepton &eOut3)
{
   // Calculates the e+e- pair production rate on a free electron target,
   // including incident photon polarization effects, for a given set of
   // kinematics.  The kinematics are specified by the initial photon
   // energy kin, the mass of the e+e- pair M, the recoil momentum vector
   // qR, the azimuthal angle of the plane containing the e+e- pair phi+,
   // and the energy of the pair positron E+. The returned value is the
   // differential cross section measured in barns/GeV^4 per atom, with
   // fConverterZ electrons per atom, differential in
   //    (d^3 qR dphi+ dE+) = (M / 2 kin) (dM dqR^2 dphiR dphi+ dE+).

   TPhoton g0(gIn);
   TLepton e0(zeroVector, mElectron);
   TLepton p1(pOut);
   TLepton e2(eOut2);
   TLepton e3(eOut3);

   // Avoid double-counting due to identical fs electrons by requiring
   // that e3 (recoil e-) have a lower momentum magnitude than e2.
   if (e2.Mass() == e3.Mass() && e3.Mom().Length() > e2.Mom().Length())
      return 0;

   // Set the initial,final polarizations
   e0.SetPol(zeroVector);
   p1.AllPol();
   e2.AllPol();
   e3.AllPol();

   // Correct the basic cross section by the converter atomic form factor
   LDouble_t result = TCrossSection::TripletProduction(g0, e0, p1, e2, e3);
   result *= fConverterZ * (1 - sqr(FFatomic(e3.Mom().Length())));
   return result * 1e-6;
}

#endif
