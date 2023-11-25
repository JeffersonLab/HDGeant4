//
// GlueXHitCCALblock - class header
//
// author: richard.t.jones at uconn.edu
// version: october 28, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitCCALblock_h
#define GlueXHitCCALblock_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitCCALblock : public G4VHit
{
 public:
   GlueXHitCCALblock() {}
   GlueXHitCCALblock(G4int column, G4int row);
   GlueXHitCCALblock(const GlueXHitCCALblock &src);
   int operator==(const GlueXHitCCALblock &right) const;
   GlueXHitCCALblock &operator+=(const GlueXHitCCALblock &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int column_;          // CCal block column, from 1 increasing x
   G4int row_;             // CCAL block row, from 1 increasing y
   G4int overflow_;

   struct hitinfo_t {
      G4double E_GeV;      // energy deposition (GeV)
      G4double t_ns;       // pulse leading-edge time (ns)
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(column_, row_); }
   static G4int GetKey(G4int column, G4int row) {
      return ((row + 1) << 16) + column + 1;
   }
};

typedef G4THitsMap<GlueXHitCCALblock> GlueXHitsMapCCALblock;

extern G4ThreadLocal G4Allocator<GlueXHitCCALblock>* GlueXHitCCALblockAllocator;

inline void* GlueXHitCCALblock::operator new(size_t)
{
   if (!GlueXHitCCALblockAllocator)
      GlueXHitCCALblockAllocator = new G4Allocator<GlueXHitCCALblock>;
   return (void *) GlueXHitCCALblockAllocator->MallocSingle();
}

inline void GlueXHitCCALblock::operator delete(void *aHit)
{
   GlueXHitCCALblockAllocator->FreeSingle((GlueXHitCCALblock*) aHit);
}

#endif
