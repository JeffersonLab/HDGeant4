//
// GlueXHitPStile - class header
//
// author: richard.t.jones at uconn.edu
// version: october 26, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitPStile_h
#define GlueXHitPStile_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitPStile : public G4VHit
{
 public:
   GlueXHitPStile() {}
   GlueXHitPStile(G4int arm, G4int column);
   GlueXHitPStile(const GlueXHitPStile &src);
   int operator==(const GlueXHitPStile &right) const;
   GlueXHitPStile &operator+=(const GlueXHitPStile &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int arm_;             // =0 for south (left), =1 for north (right)
   G4int column_;          // column number, from 1 increasing with x
   G4int overflow_;

   struct hitinfo_t {
      G4double dE_GeV;     // energy deposition (GeV)
      G4double t_ns;       // pulse leading-edge time (ns)
      G4double itrack_;    // track index of first particle making this hit
      G4double ptype_G3;   // G3 type of first particle making this hit
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(arm_, column_); }
   static G4int GetKey(G4int arm, G4int column) {
      return (arm << 10) + column;
   }
};

typedef G4THitsMap<GlueXHitPStile> GlueXHitsMapPStile;

extern G4ThreadLocal G4Allocator<GlueXHitPStile>* GlueXHitPStileAllocator;

inline void* GlueXHitPStile::operator new(size_t)
{
   if (!GlueXHitPStileAllocator)
      GlueXHitPStileAllocator = new G4Allocator<GlueXHitPStile>;
   return (void *) GlueXHitPStileAllocator->MallocSingle();
}

inline void GlueXHitPStile::operator delete(void *aHit)
{
   GlueXHitPStileAllocator->FreeSingle((GlueXHitPStile*) aHit);
}

#endif
