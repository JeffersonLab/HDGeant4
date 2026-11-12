//
// GlueXHitFDCpoint - class header
//
// author: richard.t.jones at uconn.edu
// version: october 11, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitFDCpoint_h
#define GlueXHitFDCpoint_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitFDCpoint : public G4VHit
{
 public:
   GlueXHitFDCpoint() {}
   GlueXHitFDCpoint(int chamber);
   GlueXHitFDCpoint(const GlueXHitFDCpoint &src);
   int operator==(const GlueXHitFDCpoint &right) const;
   GlueXHitFDCpoint &operator+=(const GlueXHitFDCpoint &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   G4int chamber_ = 0; // FDC chamber, count from upstream starting at 1

   G4double E_GeV = 0;       // energy of primary track at point
   G4double dEdx_GeV_cm = 0; // dE/dx (GeV/cm) of this track inside straw
   G4double dradius_cm = 0;  // track distance (cm) of closest approach to wire
   G4bool primary_ = false;  // true if track belongs to from a primary particle
   G4int ptype_G3 = 0;       // G3 type of particle making this track
   G4double px_GeV = 0;      // momentum (GeV/c) of track at point, x component
   G4double py_GeV = 0;      // momentum (GeV/c) of track at point, y component
   G4double pz_GeV = 0;      // momentum (GeV/c) of track at point, z component
   G4double x_cm = 0;        // global x coordinate of track at point (cm)
   G4double y_cm = 0;        // global y coordinate of track at point (cm)
   G4double z_cm = 0;        // global z coordinate of track at point (cm)
   G4double t_ns = 0;        // time of track crossing at point (ns)
   G4int track_ = 0;         // Geant4 track ID of particle making this track
   G4int trackID_ = 0;       // GlueX-assigned track ID of particle making this track

   G4int GetKey() const { return (track_ << 20) + int(t_ns * 100); }
};

typedef G4THitsMap<GlueXHitFDCpoint> GlueXHitsMapFDCpoint;

extern G4ThreadLocal G4Allocator<GlueXHitFDCpoint>* GlueXHitFDCpointAllocator;

inline void* GlueXHitFDCpoint::operator new(size_t)
{
   if (!GlueXHitFDCpointAllocator)
      GlueXHitFDCpointAllocator = new G4Allocator<GlueXHitFDCpoint>;
   return (void *) GlueXHitFDCpointAllocator->MallocSingle();
}

inline void GlueXHitFDCpoint::operator delete(void *aHit)
{
   GlueXHitFDCpointAllocator->FreeSingle((GlueXHitFDCpoint*) aHit);
}

#endif
