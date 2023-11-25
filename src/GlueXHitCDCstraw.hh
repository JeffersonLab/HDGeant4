//
// GlueXHitCDCstraw - class header
//
// author: richard.t.jones at uconn.edu
// version: august 28, 2015
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitCDCstraw_h
#define GlueXHitCDCstraw_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class GlueXHitCDCstraw : public G4VHit
{
 public:
   GlueXHitCDCstraw() {}
   GlueXHitCDCstraw(G4int ring, G4int sector);
   GlueXHitCDCstraw(const GlueXHitCDCstraw &src);
   int operator==(const GlueXHitCDCstraw &right) const;
   GlueXHitCDCstraw &operator+=(const GlueXHitCDCstraw &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int ring_;            // straw ring, numbered from inside starting at 1
   G4int sector_;          // straw number within ring, from 1 at/past phi=0
   G4int overflow_;

   struct hitinfo_t {
      G4double track_;     // G4 track number making this hit
      G4double q_fC;       // pulse integral (fC)
      G4double t_ns;       // pulse leading-edge time (ns)
      G4double d_cm;       // distance (cm) of closest cluster to the wire
      G4double itrack_;    // track index of first particle making this hit
      G4double ptype_G3;   // G3 type of first particle making this hit
      G4double t0_ns;      // time at start of track segment
      G4double t1_ns;      // time at end of track segment
      G4ThreeVector x0_g;  // global coordinates of start of track segment
      G4ThreeVector x1_g;  // global coordinates of end of track segment
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(ring_, sector_); }
   static G4int GetKey(G4int ring, G4int sector) {
      return (ring << 20) + sector;
   }
};

typedef G4THitsMap<GlueXHitCDCstraw> GlueXHitsMapCDCstraw;

extern G4ThreadLocal G4Allocator<GlueXHitCDCstraw>* GlueXHitCDCstrawAllocator;

inline void* GlueXHitCDCstraw::operator new(size_t)
{
   if (!GlueXHitCDCstrawAllocator)
      GlueXHitCDCstrawAllocator = new G4Allocator<GlueXHitCDCstraw>;
   return (void *) GlueXHitCDCstrawAllocator->MallocSingle();
}

inline void GlueXHitCDCstraw::operator delete(void *aHit)
{
   GlueXHitCDCstrawAllocator->FreeSingle((GlueXHitCDCstraw*) aHit);
}

#endif
