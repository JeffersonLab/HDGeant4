//
// GlueXHitUPVbar - class header
//
// author: richard.t.jones at uconn.edu
// version: october 29, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitUPVbar_h
#define GlueXHitUPVbar_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitUPVbar : public G4VHit
{
 public:
   GlueXHitUPVbar() {}
   GlueXHitUPVbar(G4int layer, G4int row);
   GlueXHitUPVbar(const GlueXHitUPVbar &src);
   int operator==(const GlueXHitUPVbar &right) const;
   GlueXHitUPVbar &operator+=(const GlueXHitUPVbar &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int layer_;           // =0 (downstr, vertical) or =1 (upstr, horizontal)
   G4int row_;             // row number, bottom-top, south-north (see hdds)
   G4int overflow_;

   struct hitinfo_t {
      G4int end_;          // end=0: top, north/left; end=1: bottom, south/right
      G4double E_GeV;      // energy deposition by track(GeV)
      G4double t_ns;       // pulse leading-edge time (ns)
      G4double xlocal_cm;  // x coordinate of the hit in local refsys (cm)
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(layer_, row_); }
   static G4int GetKey(G4int layer, G4int row) {
      return (layer << 16) + row;
   }
};

typedef G4THitsMap<GlueXHitUPVbar> GlueXHitsMapUPVbar;

extern G4ThreadLocal G4Allocator<GlueXHitUPVbar>* GlueXHitUPVbarAllocator;

inline void* GlueXHitUPVbar::operator new(size_t)
{
   if (!GlueXHitUPVbarAllocator)
      GlueXHitUPVbarAllocator = new G4Allocator<GlueXHitUPVbar>;
   return (void *) GlueXHitUPVbarAllocator->MallocSingle();
}

inline void GlueXHitUPVbar::operator delete(void *aHit)
{
   GlueXHitUPVbarAllocator->FreeSingle((GlueXHitUPVbar*) aHit);
}

#endif
