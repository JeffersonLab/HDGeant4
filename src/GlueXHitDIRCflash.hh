//
// GlueXHitDIRCflash - class header
//
// author: richard.t.jones at uconn.edu
// version: november 21, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.
//
///////////////////////////////////////////////////////////////////
// NOTICE - This is not a real hits class for the DIRC detector.
// It was created during the DIRC conceptual design phase to measure
// the impact of the DIRC on GlueX physics channels. This needs to
// be retired as soon as the technical design is complete, and in
// its place a new class should be introduced with a name like
// GlueXHitDIRCtube that represents a single detected DIRC photon.
//
// DO NOT EDIT THIS FILE TO IMPLEMENT A REALISTIC DIRC HITS SCHEME,
// clone it instead, and give the new class a more sensible name!
// You also need to add a tag in event.xml to contain the new hit.
///////////////////////////////////////////////////////////////////

#ifndef GlueXHitDIRCflash_h
#define GlueXHitDIRCflash_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitDIRCflash : public G4VHit
{
 public:
   GlueXHitDIRCflash(G4int bar=0);
   GlueXHitDIRCflash(const GlueXHitDIRCflash &src);
   int operator==(const GlueXHitDIRCflash &right) const;
   GlueXHitDIRCflash &operator+=(const GlueXHitDIRCflash &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int bar_;             // place holder, may not be used

   struct hitinfo_t {
      G4double E_GeV;      // track particle total energy (GeV)
      G4double t_ns;       // pulse leading-edge time (ns)
      G4double x_cm;       // x coordinate where flash was created
      G4double y_cm;       // y coordinate where flash was created
      G4double z_cm;       // z coordinate where flash was created
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(bar_); }
   static G4int GetKey(G4int bar) {
      return bar;
   }
};

typedef G4THitsMap<GlueXHitDIRCflash> GlueXHitsMapDIRCflash;

extern G4ThreadLocal G4Allocator<GlueXHitDIRCflash>* GlueXHitDIRCflashAllocator;

inline void* GlueXHitDIRCflash::operator new(size_t)
{
   if (!GlueXHitDIRCflashAllocator)
      GlueXHitDIRCflashAllocator = new G4Allocator<GlueXHitDIRCflash>;
   return (void *) GlueXHitDIRCflashAllocator->MallocSingle();
}

inline void GlueXHitDIRCflash::operator delete(void *aHit)
{
   GlueXHitDIRCflashAllocator->FreeSingle((GlueXHitDIRCflash*) aHit);
}

#endif
