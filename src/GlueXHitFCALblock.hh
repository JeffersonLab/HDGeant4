//
// GlueXHitFCALblock - class header
//
// author: richard.t.jones at uconn.edu
// version: october 25, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitFCALblock_h
#define GlueXHitFCALblock_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitFCALblock : public G4VHit
{
 public:
   GlueXHitFCALblock() {}
   GlueXHitFCALblock(G4int column, G4int row);
   GlueXHitFCALblock(const GlueXHitFCALblock &src);
   int operator==(const GlueXHitFCALblock &right) const;
   GlueXHitFCALblock &operator+=(const GlueXHitFCALblock &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int column_;          // FCal block column, from 1 increasing x
   G4int row_;             // FCAL block row, from 1 increasing y
   G4int overflow_;

   struct hitinfo_t {
      G4double E_GeV;      // energy deposition (GeV)
      G4double t_ns;       // pulse leading-edge time (ns)
      G4double dE_lightguide_GeV;  // light guide energy deposition (GeV)
      G4double t_lightguide_ns;    // light guide pulse leading-edge time (ns)
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(column_, row_); }
   static G4int GetKey(G4int column, G4int row) {
      return ((row + 1) << 16) + column + 1;
   }
};

typedef G4THitsMap<GlueXHitFCALblock> GlueXHitsMapFCALblock;

extern G4ThreadLocal G4Allocator<GlueXHitFCALblock>* GlueXHitFCALblockAllocator;

inline void* GlueXHitFCALblock::operator new(size_t)
{
   if (!GlueXHitFCALblockAllocator)
      GlueXHitFCALblockAllocator = new G4Allocator<GlueXHitFCALblock>;
   return (void *) GlueXHitFCALblockAllocator->MallocSingle();
}

inline void GlueXHitFCALblock::operator delete(void *aHit)
{
   GlueXHitFCALblockAllocator->FreeSingle((GlueXHitFCALblock*) aHit);
}

#endif
