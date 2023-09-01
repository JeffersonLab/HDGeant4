//
// GlueXHitITOFhit - class header
//
// author: staylor at jlab.org
// version: february 7, 2022
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitITOFhit_h
#define GlueXHitITOFhit_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitITOFhit : public G4VHit
{
 public:
  GlueXHitITOFhit() {}
  GlueXHitITOFhit(G4int id);
  GlueXHitITOFhit(const GlueXHitITOFhit &src);
  int operator==(const GlueXHitITOFhit &right) const;
  GlueXHitITOFhit &operator+=(const GlueXHitITOFhit &right);
  
  void *operator new(size_t);
  void operator delete(void *aHit);
  
  void Draw() const;
  void Print() const;
  
  // no reason to hide hit data
  
  G4int id_;
  struct hitinfo_t {
    G4double t_ns;       // time of track passing into active region (ns)
    G4double dE_GeV;     // energy deposition by track(GeV)
    G4double x_cm;       // hit position: x(cm)
    G4double y_cm;       //               y(cm)
  };
  std::vector<hitinfo_t> hits;

  G4int GetKey() const { return GetKey(id_); } 
  static G4int GetKey(G4int id) {
      return id;
   }
  
};

typedef G4THitsMap<GlueXHitITOFhit> GlueXHitsMapITOFhit;

extern G4ThreadLocal G4Allocator<GlueXHitITOFhit>* GlueXHitITOFhitAllocator;

inline void* GlueXHitITOFhit::operator new(size_t)
{
   if (!GlueXHitITOFhitAllocator)
      GlueXHitITOFhitAllocator = new G4Allocator<GlueXHitITOFhit>;
   return (void *) GlueXHitITOFhitAllocator->MallocSingle();
}

inline void GlueXHitITOFhit::operator delete(void *aHit)
{
   GlueXHitITOFhitAllocator->FreeSingle((GlueXHitITOFhit*) aHit);
}

#endif
