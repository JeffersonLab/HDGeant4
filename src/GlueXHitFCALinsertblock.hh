//
// GlueXHitFCALinsertblock - class header
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitFCALinsertblock_h
#define GlueXHitFCALinsertblock_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitFCALinsertblock : public G4VHit
{
 public:
   GlueXHitFCALinsertblock() {}
   GlueXHitFCALinsertblock(G4int column, G4int row);
   GlueXHitFCALinsertblock(const GlueXHitFCALinsertblock &src);
   int operator==(const GlueXHitFCALinsertblock &right) const;
   GlueXHitFCALinsertblock &operator+=(const GlueXHitFCALinsertblock &right);

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

typedef G4THitsMap<GlueXHitFCALinsertblock> GlueXHitsMapFCALinsertblock;

extern G4ThreadLocal G4Allocator<GlueXHitFCALinsertblock>* GlueXHitFCALinsertblockAllocator;

inline void* GlueXHitFCALinsertblock::operator new(size_t)
{
   if (!GlueXHitFCALinsertblockAllocator)
      GlueXHitFCALinsertblockAllocator = new G4Allocator<GlueXHitFCALinsertblock>;
   return (void *) GlueXHitFCALinsertblockAllocator->MallocSingle();
}

inline void GlueXHitFCALinsertblock::operator delete(void *aHit)
{
   GlueXHitFCALinsertblockAllocator->FreeSingle((GlueXHitFCALinsertblock*) aHit);
}

#endif
