//
// GlueXExceptionHandler class implementation
//
// author: richard.t.jones at uconn.edu
// version: september 21, 2026

#include "GlueXExceptionHandler.hh"
#include "G4StateManager.hh"

std::atomic<G4long> GlueXExceptionHandler::fGeomNav0003count(0);

GlueXExceptionHandler::GlueXExceptionHandler()
{
   // Constructing fDefaultHandler (above) registers *it* as the active
   // G4VExceptionHandler, because its own base-class constructor calls
   // G4StateManager::SetExceptionHandler(&fDefaultHandler). Re-register
   // this object here so that Notify() below -- which filters one
   // specific warning and delegates everything else to fDefaultHandler --
   // is the one Geant4 actually calls from now on.
   G4StateManager::GetStateManager()->SetExceptionHandler(this);
}

G4bool GlueXExceptionHandler::Notify(const char *originOfException,
                                      const char *exceptionCode,
                                      G4ExceptionSeverity severity,
                                      const char *description)
{
   if (severity == JustWarning &&
       G4String(exceptionCode) == "GeomNav0003" &&
       G4String(originOfException) == "G4Navigator::GetLocalExitNormal()")
   {
      G4long n = ++fGeomNav0003count;

      // Print occurrences 1-10, then every 10th up to 100, every 100th up
      // to 1000, every 1000th up to 10000, and so on: the counter itself
      // still counts every occurrence exactly, but console output stays
      // bounded even on a run with millions of them.
      G4long step = 1, m = n - 1;
      while (m >= 10) {
         m /= 10;
         step *= 10;
      }
      if (n % step == 0) {
         G4cout << "GlueXExceptionHandler note: G4Navigator asked for the "
                   "exit normal to a geometry boundary while not sitting "
                   "exactly on one (GeomNav0003, from "
                   "G4Navigator::GetLocalExitNormal()). This is expected "
                   "and harmless: it happens when a track's path is nearly "
                   "parallel to a boundary surface, and tracking continues "
                   "normally." << G4endl
                << "  Occurrence #" << n << " of this warning so far this "
                   "run. To prevent console flooding, only occurrences "
                   "1-10, then every 10th up to 100, every 100th up to "
                   "1000, every 1000th up to 10000 (etc) are printed; the "
                   "count above is exact regardless." << G4endl;
      }
      return false;   // never abort the run/event for this one
   }

   return fDefaultHandler.Notify(originOfException, exceptionCode,
                                  severity, description);
}
