//
// GlueXDetectorMessenger class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
//
// In the context of the Geant4 event-level multithreading model,
// this class is "shared", ie. has no thread-local state.

#ifndef GlueXDetectorMessenger_h
#define GlueXDetectorMessenger_h 1

#include "G4UImessenger.hh"

class GlueXDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3Vector;

class GlueXDetectorMessenger: public G4UImessenger
{
  public:
    GlueXDetectorMessenger(GlueXDetectorConstruction*);
    GlueXDetectorMessenger(const GlueXDetectorMessenger& src);
    GlueXDetectorMessenger &operator=(const GlueXDetectorMessenger& src);
   ~GlueXDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    G4String GetCurrentValue(G4UIcommand*);
    
  private:
    GlueXDetectorConstruction* myDetector;
    
    G4UIdirectory*             hdgeant4Dir;
    G4UIdirectory*             detDir;
    G4UIcmdWithADoubleAndUnit* FieldCmd;
    G4UIcmdWithADoubleAndUnit* StepMaxCmd;    
    G4UIcmdWithoutParameter*   OpenGeomCmd;
    G4UIcmdWithoutParameter*   CloseGeomCmd;
    G4UIcmdWith3Vector*        RadiatorAnglesCmd;
    G4UIcmdWith3Vector*        TargetNuclPolCmd;
    G4UIcmdWith3Vector*        TargetElecPolCmd;
};

#endif
