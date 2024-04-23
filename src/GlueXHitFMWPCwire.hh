//
// GlueXHitFMWPCwire - class header
//
// author: richard.t.jones at uconn.edu
// version: november 28, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitFMWPCwire_h
#define GlueXHitFMWPCwire_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitFMWPCwire : public G4VHit
{
 public:
   GlueXHitFMWPCwire() {}
   GlueXHitFMWPCwire(G4int layer, G4int wire);
   GlueXHitFMWPCwire(const GlueXHitFMWPCwire &src);
   int operator==(const GlueXHitFMWPCwire &right) const;
   GlueXHitFMWPCwire &operator+=(const GlueXHitFMWPCwire &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int layer_;           // layer number, from 1 upstream to down
   G4int wire_;            // wire number, from 1 low to high u

   struct hitinfo_t {
      G4double dE_keV;     // energy loss (keV)
      G4double t_ns;       // pulse leading-edge time (ns)
      G4double d_cm;      // doca to wire 
      G4int itrack_;       // number of track creating the hit
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(layer_, wire_); }
   static G4int GetKey(G4int layer, G4int wire) {
      return (layer << 10) + wire + 1;
   }
};

typedef G4THitsMap<GlueXHitFMWPCwire> GlueXHitsMapFMWPCwire;

extern G4ThreadLocal G4Allocator<GlueXHitFMWPCwire>* GlueXHitFMWPCwireAllocator;

inline void* GlueXHitFMWPCwire::operator new(size_t)
{
   if (!GlueXHitFMWPCwireAllocator)
      GlueXHitFMWPCwireAllocator = new G4Allocator<GlueXHitFMWPCwire>;
   return (void *) GlueXHitFMWPCwireAllocator->MallocSingle();
}

inline void GlueXHitFMWPCwire::operator delete(void *aHit)
{
   GlueXHitFMWPCwireAllocator->FreeSingle((GlueXHitFMWPCwire*) aHit);
}

#endif
