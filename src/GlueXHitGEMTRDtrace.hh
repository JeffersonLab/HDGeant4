//
// GlueXHitGEMTRDtrace - class header
//
// author: staylor at jlab.org
// version: february 7, 2022
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitGEMTRDtrace_h
#define GlueXHitGEMTRDtrace_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitGEMTRDtrace : public G4VHit
{
 public:
  GlueXHitGEMTRDtrace() {}
  GlueXHitGEMTRDtrace(G4int layer);
  GlueXHitGEMTRDtrace(const GlueXHitGEMTRDtrace &src);
  int operator==(const GlueXHitGEMTRDtrace &right) const;
  GlueXHitGEMTRDtrace &operator+=(const GlueXHitGEMTRDtrace &right);
  
  void *operator new(size_t);
  void operator delete(void *aHit);
  
  void Draw() const;
  void Print() const;
  
  // no reason to hide hit data
  
  G4int layer_;           // layer number: 0 above beam line, 1 below
  
  struct traceinfo_t {
    G4double dE_keV;     // energy loss (keV)
    G4double t_ns;       // time of track passing into active region (ns)
    G4double x_cm;       // Entrance position: x(cm)
    G4double y_cm;       //                    y(cm)
    G4double dxdz;       // Track direction
    G4double dydz;
    G4int itrack_;       // number of track creating the trace
  };
  std::vector<traceinfo_t> traces;

  G4int GetKey() const { return GetKey(layer_); }
  static G4int GetKey(G4int layer) {
    return layer + 1;
  }
};

typedef G4THitsMap<GlueXHitGEMTRDtrace> GlueXHitsMapGEMTRDtrace;

extern G4ThreadLocal G4Allocator<GlueXHitGEMTRDtrace>* GlueXHitGEMTRDtraceAllocator;

inline void* GlueXHitGEMTRDtrace::operator new(size_t)
{
   if (!GlueXHitGEMTRDtraceAllocator)
      GlueXHitGEMTRDtraceAllocator = new G4Allocator<GlueXHitGEMTRDtrace>;
   return (void *) GlueXHitGEMTRDtraceAllocator->MallocSingle();
}

inline void GlueXHitGEMTRDtrace::operator delete(void *aHit)
{
   GlueXHitGEMTRDtraceAllocator->FreeSingle((GlueXHitGEMTRDtrace*) aHit);
}

#endif
