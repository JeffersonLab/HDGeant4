//
// GlueXHitPSCpaddle - class header
//
// author: richard.t.jones at uconn.edu
// version: october 26, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitPSCpaddle_h
#define GlueXHitPSCpaddle_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitPSCpaddle : public G4VHit
{
 public:
   GlueXHitPSCpaddle() {}
   GlueXHitPSCpaddle(G4int arm, G4int module);
   GlueXHitPSCpaddle(const GlueXHitPSCpaddle &src);
   int operator==(const GlueXHitPSCpaddle &right) const;
   GlueXHitPSCpaddle &operator+=(const GlueXHitPSCpaddle &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int arm_;             // =0 for right/south, =1 for left/north
   G4int module_;          // counter number, from 1 increasing with x
   G4int overflow_;

   struct hitinfo_t {
      G4double dE_GeV;     // energy deposition (GeV)
      G4double t_ns;       // pulse leading-edge time (ns)
      G4double itrack_;    // track index of first particle making this hit
      G4double ptype_G3;   // G3 type of first particle making this hit
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(arm_, module_); }
   static G4int GetKey(G4int arm, G4int module) {
      return (arm << 10) + module;
   }
};

typedef G4THitsMap<GlueXHitPSCpaddle> GlueXHitsMapPSCpaddle;

extern G4ThreadLocal G4Allocator<GlueXHitPSCpaddle>* GlueXHitPSCpaddleAllocator;

inline void* GlueXHitPSCpaddle::operator new(size_t)
{
   if (!GlueXHitPSCpaddleAllocator)
      GlueXHitPSCpaddleAllocator = new G4Allocator<GlueXHitPSCpaddle>;
   return (void *) GlueXHitPSCpaddleAllocator->MallocSingle();
}

inline void GlueXHitPSCpaddle::operator delete(void *aHit)
{
   GlueXHitPSCpaddleAllocator->FreeSingle((GlueXHitPSCpaddle*) aHit);
}

#endif
