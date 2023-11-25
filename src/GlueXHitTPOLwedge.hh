//
// GlueXHitTPOLwedge - class header
//
// author: richard.t.jones at uconn.edu
// version: october 21, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitTPOLwedge_h
#define GlueXHitTPOLwedge_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitTPOLwedge : public G4VHit
{
 public:
   GlueXHitTPOLwedge() {}
   GlueXHitTPOLwedge(G4int sector, G4int ring=0);
   GlueXHitTPOLwedge(const GlueXHitTPOLwedge &src);
   int operator==(const GlueXHitTPOLwedge &right) const;
   GlueXHitTPOLwedge &operator+=(const GlueXHitTPOLwedge &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int ring_;            // ring number, from 1, inner to outer
   G4int sector_;          // sector number, from 1 at/past phi=0
   G4int overflow_;

   struct hitinfo_t {
      G4double dE_MeV;     // energy deposition (MeV)
      G4double t_ns;       // pulse leading-edge time (ns)
      G4double itrack_;    // track index of first particle making this hit
      G4double ptype_G3;   // G3 type of first particle making this hit
      G4double t0_ns;      // time of passage of the track making this hit
      G4double r_cm;       // r coordinate of the hit in global refsys
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(sector_, ring_); }
   static G4int GetKey(G4int sector, G4int ring=0) {
      return (ring << 10) + sector;
   }
};

typedef G4THitsMap<GlueXHitTPOLwedge> GlueXHitsMapTPOLwedge;

extern G4ThreadLocal G4Allocator<GlueXHitTPOLwedge>* GlueXHitTPOLwedgeAllocator;

inline void* GlueXHitTPOLwedge::operator new(size_t)
{
   if (!GlueXHitTPOLwedgeAllocator)
      GlueXHitTPOLwedgeAllocator = new G4Allocator<GlueXHitTPOLwedge>;
   return (void *) GlueXHitTPOLwedgeAllocator->MallocSingle();
}

inline void GlueXHitTPOLwedge::operator delete(void *aHit)
{
   GlueXHitTPOLwedgeAllocator->FreeSingle((GlueXHitTPOLwedge*) aHit);
}

#endif
