//
// GlueXHitGCALblock - class header
//
// author: richard.t.jones at uconn.edu
// version: october 28, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitGCALblock_h
#define GlueXHitGCALblock_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitGCALblock : public G4VHit
{
 public:
   GlueXHitGCALblock() {}
   GlueXHitGCALblock(G4int module);
   GlueXHitGCALblock(const GlueXHitGCALblock &src);
   int operator==(const GlueXHitGCALblock &right) const;
   GlueXHitGCALblock &operator+=(const GlueXHitGCALblock &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int module_;          // GCal module, from 1 starting at/after phi=0
   G4int overflow_;

   struct hitinfo_t {
      G4double E_GeV;      // energy deposition (GeV)
      G4double t_ns;       // pulse leading-edge time (ns)
      G4double zlocal_cm;  // z coordinate of hit (cm)
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(module_); }
   static G4int GetKey(G4int module) {
      return (module + 1) << 16;
   }
};

typedef G4THitsMap<GlueXHitGCALblock> GlueXHitsMapGCALblock;

extern G4ThreadLocal G4Allocator<GlueXHitGCALblock>* GlueXHitGCALblockAllocator;

inline void* GlueXHitGCALblock::operator new(size_t)
{
   if (!GlueXHitGCALblockAllocator)
      GlueXHitGCALblockAllocator = new G4Allocator<GlueXHitGCALblock>;
   return (void *) GlueXHitGCALblockAllocator->MallocSingle();
}

inline void GlueXHitGCALblock::operator delete(void *aHit)
{
   GlueXHitGCALblockAllocator->FreeSingle((GlueXHitGCALblock*) aHit);
}

#endif
