//
// GlueXHitCEREtube - class header
//
// author: richard.t.jones at uconn.edu
// version: october 21, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitCEREtube_h
#define GlueXHitCEREtube_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitCEREtube : public G4VHit
{
 public:
   GlueXHitCEREtube() {}
   GlueXHitCEREtube(G4int sector);
   GlueXHitCEREtube(const GlueXHitCEREtube &src);
   int operator==(const GlueXHitCEREtube &right) const;
   GlueXHitCEREtube &operator+=(const GlueXHitCEREtube &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int sector_;          // counter number, from 1 at/past phi=0
   G4int overflow_;

   struct hitinfo_t {
      G4double pe_;        // number of photoelectrons in hit
      G4double t_ns;       // pulse leading-edge time (ns)
      G4int itrack_;       // track number creating hit
      G4double z0_cm;      // z coordinate where photon was created
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(sector_); }
   static G4int GetKey(G4int sector) {
      return sector;
   }
};

typedef G4THitsMap<GlueXHitCEREtube> GlueXHitsMapCEREtube;

extern G4ThreadLocal G4Allocator<GlueXHitCEREtube>* GlueXHitCEREtubeAllocator;

inline void* GlueXHitCEREtube::operator new(size_t)
{
   if (!GlueXHitCEREtubeAllocator)
      GlueXHitCEREtubeAllocator = new G4Allocator<GlueXHitCEREtube>;
   return (void *) GlueXHitCEREtubeAllocator->MallocSingle();
}

inline void GlueXHitCEREtube::operator delete(void *aHit)
{
   GlueXHitCEREtubeAllocator->FreeSingle((GlueXHitCEREtube*) aHit);
}

#endif
