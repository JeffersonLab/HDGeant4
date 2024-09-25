//
// G4TRandom
//
// author: richard.t.jones at uconn.edu
// version: july 18, 2020
//
// This is thin a wrapper around the Geant4 random number generator
// that provides random numbers needed by ROOT objects within the
// Geant4 context. This allows importance sampling to be applied
// using TH1::GetRandom() without losing the essential feature of
// the simulation that a given seed will result in the same output
// event regardless of which thread works on it or the order that
// events are simulated.

#include <TRandom.h>
#include <Randomize.hh>

class G4TRandom : public TRandom {
 public:
   G4TRandom() { gRandom = this; }
   virtual ~G4TRandom() {}
   virtual  Int_t    Binomial(Int_t ntot, Double_t prob);
   virtual  Double_t BreitWigner(Double_t mean=0, Double_t gamma=1);
   virtual  void     Circle(Double_t &x, Double_t &y, Double_t r);
   virtual  Double_t Exp(Double_t tau); 
   virtual  Double_t Gaus(Double_t mean=0, Double_t sigma=1);
   virtual  UInt_t   GetSeed() const;
   virtual  UInt_t   Integer(UInt_t imax);
   virtual  Double_t Landau(Double_t mean=0, Double_t sigma=1);
   virtual  ULong64_t Poisson(Double_t mean);
   virtual  Double_t PoissonD(Double_t mean);
   virtual  void     Rannor(Float_t &a, Float_t &b);
   virtual  void     Rannor(Double_t &a, Double_t &b);
   virtual  void     ReadRandom(const char *filename);
   virtual  void     SetSeed(ULong_t seed=0);
   virtual  Double_t Rndm();
   // keep for backward compatibility
   virtual  Double_t Rndm(Int_t ) { return Rndm(); }
   virtual  void     RndmArray(Int_t n, Float_t *array);
   virtual  void     RndmArray(Int_t n, Double_t *array);
   virtual  void     Sphere(Double_t &x, Double_t &y, Double_t &z, Double_t r);
   virtual  Double_t Uniform(Double_t x1=1);
   virtual  Double_t Uniform(Double_t x1, Double_t x2);
   virtual  void     WriteRandom(const char *filename) const;
};

inline Int_t G4TRandom::Binomial(Int_t ntot, Double_t prob) {
   return CLHEP::RandBinomial::shoot(ntot, prob);
}

inline Double_t G4TRandom::BreitWigner(Double_t mean, Double_t gamma) {
   return CLHEP::RandBreitWigner::shoot(mean, gamma);
}

inline void G4TRandom::Circle(Double_t &x, Double_t &y, Double_t r) {
   double phi = CLHEP::RandFlat::shoot() * 2* M_PI;
   x = r * cos(phi);
   y = r * sin(phi);
}

inline Double_t G4TRandom::Exp(Double_t tau) {
   return CLHEP::RandExponential::shoot(tau);
}

inline Double_t G4TRandom::Gaus(Double_t mean, Double_t sigma) {
   return CLHEP::RandGauss::shoot(mean, sigma);
}

inline UInt_t G4TRandom::GetSeed() const {
   // not implemented
   return -1;
}

inline UInt_t G4TRandom::Integer(UInt_t imax) {
   return CLHEP::RandFlat::shootInt(imax);
}

inline Double_t G4TRandom::Landau(Double_t mean, Double_t sigma) {
   return CLHEP::RandLandau::shoot() * sigma + mean;
}

inline ULong64_t G4TRandom::Poisson(Double_t mean) {
   return CLHEP::RandPoisson::shoot(mean);
}

inline Double_t G4TRandom::PoissonD(Double_t mean) {
   // if you are looking for a smooth continuous interpolation
   // of the discrete Poisson distribution, this is not it!
   return CLHEP::RandPoisson::shoot(mean);
}

inline void G4TRandom::Rannor(Float_t &a, Float_t &b) {
   double aa=0, bb=0;
   CLHEP::RandGauss::shoot(aa, bb);
   a = aa;
   b = bb;
}

inline void G4TRandom::Rannor(Double_t &a, Double_t &b) {
   CLHEP::RandGauss::shoot(a, b);
}

inline void G4TRandom::ReadRandom(const char *filename) {
   // not implemented
}

inline void G4TRandom::SetSeed(ULong_t seed) {
   // not implemented
}

inline Double_t G4TRandom::Rndm() {
   return CLHEP::RandFlat::shoot();
}

inline void G4TRandom::RndmArray(Int_t n, Float_t *array) {
   double *buf = new double[n];
   CLHEP::RandFlat::shootArray(n, buf);
   for (int i=0; i < n; ++i)
      array[i] = buf[i];
   delete [] buf;
}

inline void G4TRandom::RndmArray(Int_t n, Double_t *array) {
   CLHEP::RandFlat::shootArray(n, array);
}

inline void G4TRandom::Sphere(Double_t &x, Double_t &y, Double_t &z, Double_t r) {
   double costheta = CLHEP::RandFlat::shoot() * 2 - 1;
   double phi = CLHEP::RandFlat::shoot() * 2 * M_PI;
   double sintheta = sqrt(1 - costheta*costheta);
   x = r * sintheta * cos(phi);
   y = r * sintheta * sin(phi);
   z = r * costheta;
}

inline Double_t G4TRandom::Uniform(Double_t x1) {
   return CLHEP::RandFlat::shoot(x1);
}

inline Double_t G4TRandom::Uniform(Double_t x1, Double_t x2) {
   return CLHEP::RandFlat::shoot(x1, x2);
}

inline void G4TRandom::WriteRandom(const char *filename) const {
   // not implemented
}
