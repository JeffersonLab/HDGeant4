//
// GlueXExceptionHandler class header
//
// author: richard.t.jones at uconn.edu
// version: september 21, 2026
//
// Custom G4VExceptionHandler that rate-limits the well-known, harmless
// G4Navigator::GetLocalExitNormal() "GeomNav0003 / NOT at a Boundary"
// warning. This fires when the field-propagation boundary locator queries
// the exit normal at a trial point where the navigator does not consider
// itself to be sitting exactly on a boundary -- which happens routinely
// when a track's path is nearly parallel (grazing/tangential) to a
// geometry boundary. It is purely diagnostic (always JustWarning severity)
// and the caller already falls back gracefully, so nothing is lost by
// throttling how much of it reaches the console. Every other exception or
// warning is passed through unchanged to Geant4's normal handling.
//
// In the context of the Geant4 event-level multithreading model, one
// instance of this class must be installed on each thread (master and
// workers), since G4StateManager/exception handler registration is
// thread-local. The occurrence counter is shared (process-wide) across
// all threads via a single std::atomic, so GetGeomNav0003Count() reports
// the true total across the whole run regardless of which thread asks.

#ifndef GlueXExceptionHandler_h
#define GlueXExceptionHandler_h 1

#include <atomic>
#include "G4VExceptionHandler.hh"
#include "G4ExceptionHandler.hh"
#include "G4ExceptionSeverity.hh"
#include "globals.hh"

class GlueXExceptionHandler : public G4VExceptionHandler
{
 public:
   GlueXExceptionHandler();
   ~GlueXExceptionHandler() {}

   virtual G4bool Notify(const char *originOfException,
                          const char *exceptionCode,
                          G4ExceptionSeverity severity,
                          const char *description);

   static G4long GetGeomNav0003Count() { return fGeomNav0003count.load(); }

 private:
   G4ExceptionHandler fDefaultHandler;
   static std::atomic<G4long> fGeomNav0003count;
};

#endif
