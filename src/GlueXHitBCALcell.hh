//
// GlueXHitBCALcell - class header
//
// author: richard.t.jones at uconn.edu
// version: october 24, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitBCALcell_h
#define GlueXHitBCALcell_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitBCALcell : public G4VHit
{
 public:
   GlueXHitBCALcell() {}
   GlueXHitBCALcell(G4int module, G4int layer, G4int sector);
   GlueXHitBCALcell(const GlueXHitBCALcell &src);
   int operator==(const GlueXHitBCALcell &right) const;
   GlueXHitBCALcell &operator+=(const GlueXHitBCALcell &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int module_;          // BCal module number, from 1 at/past phi=0
   G4int layer_;           // cell layer number, from 1 inner - outer
   G4int sector_;          // cell sector number, from 1 advancing phi

   struct hitinfo_t {
      G4double E_GeV;       // energy deposition (GeV)
      G4double t_ns;        // pulse leading-edge time (ns)
      G4double zlocal_cm;   // z coordinate of the hit in local refsys
      G4double incidentId_; // id of particle that generated this shower
      G4double Eup_GeV;     // upstream end energy deposition (GeV)
      G4double Edown_GeV;   // downstream end energy deposition (GeV)
      G4double tup_ns;      // upstream end pulse leading-edge time (ns)
      G4double tdown_ns;    // downstream end pulse leading-edge time (ns)
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(module_, layer_, sector_); }
   static G4int GetKey(G4int module, G4int layer, G4int sector) {
      return (module << 16) + (layer << 9) + sector;
   }
};

typedef G4THitsMap<GlueXHitBCALcell> GlueXHitsMapBCALcell;

extern G4ThreadLocal G4Allocator<GlueXHitBCALcell>* GlueXHitBCALcellAllocator;

inline void* GlueXHitBCALcell::operator new(size_t)
{
   if (!GlueXHitBCALcellAllocator)
      GlueXHitBCALcellAllocator = new G4Allocator<GlueXHitBCALcell>;
   return (void *) GlueXHitBCALcellAllocator->MallocSingle();
}

inline void GlueXHitBCALcell::operator delete(void *aHit)
{
   GlueXHitBCALcellAllocator->FreeSingle((GlueXHitBCALcell*) aHit);
}

#endif
