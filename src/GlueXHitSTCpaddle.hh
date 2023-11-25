//
// GlueXHitSTCpaddle - class header
//
// author: richard.t.jones at uconn.edu
// version: october 21, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitSTCpaddle_h
#define GlueXHitSTCpaddle_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitSTCpaddle : public G4VHit
{
 public:
   GlueXHitSTCpaddle() {}
   GlueXHitSTCpaddle(G4int sector);
   GlueXHitSTCpaddle(const GlueXHitSTCpaddle &src);
   int operator==(const GlueXHitSTCpaddle &right) const;
   GlueXHitSTCpaddle &operator+=(const GlueXHitSTCpaddle &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int sector_;          // counter number, from 1 at/past phi=0
   G4int overflow_;

   struct hitinfo_t {
      G4double dE_MeV;     // energy deposition (MeV)
      G4double t_ns;       // pulse leading-edge time (ns)
      G4double itrack_;    // track index of first particle making this hit
      G4double ptype_G3;   // G3 type of first particle making this hit
      G4double t0_ns;      // time of passage of the track making this hit
      G4double z_cm;       // z coordinate of the hit in global refsys
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(sector_); }
   static G4int GetKey(G4int sector) {
      return sector;
   }
};

typedef G4THitsMap<GlueXHitSTCpaddle> GlueXHitsMapSTCpaddle;

extern G4ThreadLocal G4Allocator<GlueXHitSTCpaddle>* GlueXHitSTCpaddleAllocator;

inline void* GlueXHitSTCpaddle::operator new(size_t)
{
   if (!GlueXHitSTCpaddleAllocator)
      GlueXHitSTCpaddleAllocator = new G4Allocator<GlueXHitSTCpaddle>;
   return (void *) GlueXHitSTCpaddleAllocator->MallocSingle();
}

inline void GlueXHitSTCpaddle::operator delete(void *aHit)
{
   GlueXHitSTCpaddleAllocator->FreeSingle((GlueXHitSTCpaddle*) aHit);
}

#endif
