//
// PairConversionGeneration class header
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
// free recoil electron are included.
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

#ifdef USING_DIRACXX
#ifndef PairConversionGeneration_h
#define PairConversionGeneration_h 1

#include "TPhoton.h"
#include "TLepton.h"

#include <vector>

#if BOOST_PYTHON_WRAPPING
#include <boost/python.hpp>
#endif

class PairConversionGeneration {
 public:
   PairConversionGeneration(std::vector<double> Z,
                            std::vector<double> A,
                            std::vector<double> w);
   ~PairConversionGeneration();

   LDouble_t FFatomic(LDouble_t qRecoil);
   LDouble_t DiffXS_pair(const TPhoton &gIn, 
                         const TLepton &pOut, const TLepton &eOut);
   LDouble_t DiffXS_triplet(const TPhoton &gIn, const TLepton &pOut,
                            const TLepton &eOut2, const TLepton &eOut3);

   const TThreeVectorReal &GetPolarization();
   unsigned int GetConverterZ();
   void SetConverterZ(unsigned int Z);

 protected:
   unsigned int fConverterZ;  // atomic number of converter material

 private:
   PairConversionGeneration(const PairConversionGeneration &src);
   PairConversionGeneration &operator=(const PairConversionGeneration &src);
};

inline unsigned int PairConversionGeneration::GetConverterZ()
{
   return fConverterZ;
}

inline void PairConversionGeneration::SetConverterZ(unsigned int Z)
{
   fConverterZ = Z;
}

#endif
#endif
