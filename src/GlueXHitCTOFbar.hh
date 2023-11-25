//
// GlueXHitCTOFbar - class header
//
// author: staylor at jlab.org
// version: october 25, 2021
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitCTOFbar_h
#define GlueXHitCTOFbar_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitCTOFbar : public G4VHit
{
 public:
   GlueXHitCTOFbar() {}
   GlueXHitCTOFbar(G4int bar);
   GlueXHitCTOFbar(const GlueXHitCTOFbar &src);
   int operator==(const GlueXHitCTOFbar &right) const;
   GlueXHitCTOFbar &operator+=(const GlueXHitCTOFbar &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int bar_;             // bar number (see hdds)
   G4int overflow_;

   struct hitinfo_t {
      G4int end_;          // end=0: top, end=1: bottom
      G4double dE_GeV;     // energy deposition by track(GeV)
      G4double t_ns;       // pulse leading-edge time (ns)
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(bar_); }
   static G4int GetKey(G4int bar) {
      return bar;
   }
};

typedef G4THitsMap<GlueXHitCTOFbar> GlueXHitsMapCTOFbar;

extern G4ThreadLocal G4Allocator<GlueXHitCTOFbar>* GlueXHitCTOFbarAllocator;

inline void* GlueXHitCTOFbar::operator new(size_t)
{
   if (!GlueXHitCTOFbarAllocator)
      GlueXHitCTOFbarAllocator = new G4Allocator<GlueXHitCTOFbar>;
   return (void *) GlueXHitCTOFbarAllocator->MallocSingle();
}

inline void GlueXHitCTOFbar::operator delete(void *aHit)
{
   GlueXHitCTOFbarAllocator->FreeSingle((GlueXHitCTOFbar*) aHit);
}

#endif
