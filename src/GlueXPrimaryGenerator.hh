//
// GlueXPrimaryGenerator class header
//
// author: richard.t.jones at uconn.edu
// version: december 24, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "shared", ie. has no thread-local state. It is
// invoked from the GlueXPrimaryGeneratorAction class, which is
// responsible for managing the interlocks to ensure thread safety.

#ifndef GlueXPrimaryGenerator_H
#define GlueXPrimaryGenerator_H

#include <G4VPrimaryGenerator.hh>
#include <G4Event.hh>
#include <HDDM/hddm_s.hpp>

class GlueXPrimaryGenerator: public G4VPrimaryGenerator
{
 public:
   GlueXPrimaryGenerator(hddm_s::istream *hddm_source);
   virtual ~GlueXPrimaryGenerator();

   virtual void GeneratePrimaryVertex(G4Event *event);

 protected:
   hddm_s::istream *fHDDMistream;

 private:
   GlueXPrimaryGenerator(const GlueXPrimaryGenerator &src) = delete;
   GlueXPrimaryGenerator &operator=(const GlueXPrimaryGenerator &src) = delete;
};

#endif
