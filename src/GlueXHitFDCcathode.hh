//
// GlueXHitFDCcathode - class header
//
// author: richard.t.jones at uconn.edu
// version: october 11, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitFDCcathode_h
#define GlueXHitFDCcathode_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitFDCcathode : public G4VHit
{
 public:
   GlueXHitFDCcathode() {}
   GlueXHitFDCcathode(G4int chamber, G4int plane, G4int strip);
   GlueXHitFDCcathode(const GlueXHitFDCcathode &src);
   int operator==(const GlueXHitFDCcathode &right) const;
   GlueXHitFDCcathode &operator+=(const GlueXHitFDCcathode &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int chamber_ = 0; // chamber number, from upstream starting at 1
   G4int plane_ = 0;   // cathode plane number, from upstream starting at 1
   G4int strip_ = 0;   // cathode strip number, ordered by u starting at 1

   struct hitinfo_t {
      G4double q_fC = 0;     // pulse integral (fC)
      G4double t_ns = 0;     // pulse leading-edge time (ns)
      G4double itrack_ = 0;  // track index of first particle making this hit
      G4double ptype_G3 = 0; // G3 type of first particle making this hit
      G4double u_cm = 0;     // location of hit wire in adjacent wire plane
      G4double v_cm = 0;     // location of avalanche along wire in adjacent plane
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(chamber_, plane_, strip_); }
   static G4int GetKey(G4int chamber, G4int plane, G4int strip) {
      return (chamber << 20) + (plane << 10) + strip;
   }
};

typedef G4THitsMap<GlueXHitFDCcathode> GlueXHitsMapFDCcathode;

extern G4ThreadLocal G4Allocator<GlueXHitFDCcathode>* GlueXHitFDCcathodeAllocator;

inline void* GlueXHitFDCcathode::operator new(size_t)
{
   if (!GlueXHitFDCcathodeAllocator)
      GlueXHitFDCcathodeAllocator = new G4Allocator<GlueXHitFDCcathode>;
   return (void *) GlueXHitFDCcathodeAllocator->MallocSingle();
}

inline void GlueXHitFDCcathode::operator delete(void *aHit)
{
   GlueXHitFDCcathodeAllocator->FreeSingle((GlueXHitFDCcathode*) aHit);
}

#endif
