//
// GlueXHitGEMTRDhit - class header
//
// author: staylor at jlab.org
// version: february 7, 2022
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitGEMTRDhit_h
#define GlueXHitGEMTRDhit_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitGEMTRDhit : public G4VHit
{
 public:
  GlueXHitGEMTRDhit() {}
  GlueXHitGEMTRDhit(G4int layer);
  GlueXHitGEMTRDhit(const GlueXHitGEMTRDhit &src);
  int operator==(const GlueXHitGEMTRDhit &right) const;
  GlueXHitGEMTRDhit &operator+=(const GlueXHitGEMTRDhit &right);
  
  void *operator new(size_t);
  void operator delete(void *aHit);
  
  void Draw() const;
  void Print() const;
  
  // no reason to hide hit data
  
  G4int layer_;           // layer number: 0 above beam line, 1 below
  
  struct hitinfo_t {
    G4double t_ns;       // time of track passing into active region (ns)
    G4double q_fC;       // cluster charge (fC)
    G4double x_cm;       // cluster position: x(cm)
    G4double y_cm;       //                   y(cm)
    G4double d_cm;       // Drift distance (cm)
    G4int itrack_;       // number of track creating the hit
  };
  std::vector<hitinfo_t> hits;

  G4int GetKey() const { return GetKey(layer_); }
  static G4int GetKey(G4int layer) {
    return layer + 1;
  }
};

typedef G4THitsMap<GlueXHitGEMTRDhit> GlueXHitsMapGEMTRDhit;

extern G4ThreadLocal G4Allocator<GlueXHitGEMTRDhit>* GlueXHitGEMTRDhitAllocator;

inline void* GlueXHitGEMTRDhit::operator new(size_t)
{
   if (!GlueXHitGEMTRDhitAllocator)
      GlueXHitGEMTRDhitAllocator = new G4Allocator<GlueXHitGEMTRDhit>;
   return (void *) GlueXHitGEMTRDhitAllocator->MallocSingle();
}

inline void GlueXHitGEMTRDhit::operator delete(void *aHit)
{
   GlueXHitGEMTRDhitAllocator->FreeSingle((GlueXHitGEMTRDhit*) aHit);
}

#endif
