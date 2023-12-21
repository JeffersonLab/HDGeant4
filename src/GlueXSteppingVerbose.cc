//
// GlueXSteppingVerbose class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXSteppingVerbose.hh"

#include "G4TransportationManager.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4Navigator.hh"
#include "G4UnitsTable.hh"

#define G4setw std::setw

G4Mutex GlueXSteppingVerbose::fMutex = G4MUTEX_INITIALIZER;

void GlueXSteppingVerbose::StepInfo()
{
   G4AutoLock barrier(&fMutex);
   CopyState();
   G4int prec = G4cout.precision(5);
 
   if (verboseLevel >= 1) {
     if (verboseLevel >= 4 ) VerboseTrack();
     if (verboseLevel >= 3) {
       G4cout << G4endl
              << G4setw( 5) << "#Step#"     << " "
              << G4setw( 8) << "X"          << "    "
              << G4setw( 8) << "Y"          << "    "
              << G4setw( 8) << "Z"          << "    "
              << G4setw(11) << "KineE"      << " "
              << G4setw(11) << "dEStep"     << " "
              << G4setw(12) << "StepLeng"
              << G4setw(12) << "TrakLeng"
              << G4setw(12) << "Volume"    << "  "
              << G4setw(12) << "Process"   << "  "
              << G4setw(12) << "Status"
              << G4endl;
     }
 
     const G4Step *step = G4ParallelWorldProcess::GetHyperStep();
     G4VPhysicalVolume *pvol = step->GetPostStepPoint()->GetPhysicalVolume();
     G4String volname = (pvol)? pvol->GetName() : "NULL";
     int copyno = (pvol)? pvol->GetCopyNo() : 0;
     G4cout << G4setw(5) << fTrack->GetCurrentStepNumber() << " "
            << G4setw(8) << G4BestUnit(fTrack->GetPosition().x(),"Length")
            << G4setw(8) << G4BestUnit(fTrack->GetPosition().y(),"Length")
            << G4setw(8) << G4BestUnit(fTrack->GetPosition().z(),"Length")
            << G4setw(8) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
            << G4setw(8) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
            << G4setw(8) << G4BestUnit(fStep->GetStepLength(),"Length")
            << G4setw(8) << G4BestUnit(fTrack->GetTrackLength(),"Length")
            << G4setw(12) << volname << ":" << copyno;
 
     const G4VProcess* process
                       = fStep->GetPostStepPoint()->GetProcessDefinedStep();
     G4String procName = " UserLimit";
     if (process) procName = process->GetProcessName();
     if (fStepStatus == fWorldBoundary) procName = "OutOfWorld";
     G4cout << "   " << G4setw(10) << procName;
     G4String stepstat;
     G4StepStatus stepStatus = fStep->GetPostStepPoint()->GetStepStatus();
     if (stepStatus == fWorldBoundary)
        stepstat = "WorldBoundary";
     else if (stepStatus == fGeomBoundary)
        stepstat = "GeomBoundary";
     else if (stepStatus == fAtRestDoItProc)
        stepstat = "AtRestDoItProc";
     else if (stepStatus == fAlongStepDoItProc)
        stepstat = "AlongStepDoItProc";
     else if (stepStatus == fPostStepDoItProc)
        stepstat = "PostStepDoItProc";
     else if (stepStatus == fUserDefinedLimit)
        stepstat = "UserDefinedLimit";
     else if (stepStatus == fExclusivelyForcedProc)
        stepstat = "ExclusivelyForcedProc";
     else
        stepstat == "Undefined";
     G4cout << "   " << G4setw(10) << stepstat;
     G4cout << G4endl;
 
     if (verboseLevel == 2) {
       G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
                        fN2ndariesAlongStepDoIt +
                        fN2ndariesPostStepDoIt;
       if(tN2ndariesTot>0){
    G4cout << "\n    :----- List of secondaries ----------------"
           << G4endl;
         G4cout.precision(4);
    for(size_t lp1=(*fSecondary).size()-tN2ndariesTot;
                         lp1<(*fSecondary).size(); lp1++){
      G4cout << "   "
             << G4setw(13)
             << (*fSecondary)[lp1]->GetDefinition()->GetParticleName()
             << ":  energy ="
             << G4setw(6)
             << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
             << "  time ="
             << G4setw(6)
             << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(),"Time")
             << G4endl;
    }
 
    G4cout << "    :------------------------------------------\n"
           << G4endl;
       }
     }
 
   }

   G4cout.precision(prec);
}

void GlueXSteppingVerbose::TrackingStarted()
{
   G4AutoLock barrier(&fMutex);
   CopyState();
   G4int prec = G4cout.precision(5);

   if (verboseLevel > 0) {
     G4cout << G4setw( 5) << "Step#"      << " "
            << G4setw( 8) << "X"          << "    "
            << G4setw( 8) << "Y"          << "    "
            << G4setw( 8) << "Z"          << "    "
            << G4setw(11) << "KineE"      << " "
            << G4setw(11) << "dEStep"     << " "
            << G4setw(12) << "StepLeng"
            << G4setw(12) << "TrakLeng"
            << G4setw(12) << "Volume"     << "  "
            << G4setw(12) << "Process"    << "  "
            << G4setw(12) << "Status"
            << G4endl;
 
     // Navigation is not yet initialized for this track, so drill down
     // into the list of navigators and find the layer where this point lies.
 
     G4TouchableHistory *touch = 0;
     G4TransportationManager *tmanager =
                         G4TransportationManager::GetTransportationManager();
     std::vector<G4VPhysicalVolume*>::iterator iter = 
                                               tmanager->GetWorldsIterator();
     for (int world = tmanager->GetNoWorlds() - 1; world >= 0; --world) {
       G4Navigator *navigator = tmanager->GetNavigator(iter[world]);
       G4VPhysicalVolume *pvol = navigator->
                  LocateGlobalPointAndSetup(fTrack->GetPosition(),0,false);
       G4LogicalVolume *lvol = (pvol)? pvol->GetLogicalVolume() : 0;
       if (lvol && lvol->GetMaterial()) {
         touch = navigator->CreateTouchableHistory();
         break;
       }
     }
     G4VPhysicalVolume *pvol = (touch)? touch->GetVolume() : 0;
     G4String volname = (pvol)? pvol->GetName() : "NULL";
     int copyno = (touch)? touch->GetCopyNumber() : 0;
     G4cout << G4setw(5) << fTrack->GetCurrentStepNumber() << " "
            << G4setw(8) << G4BestUnit(fTrack->GetPosition().x(),"Length")
            << G4setw(8) << G4BestUnit(fTrack->GetPosition().y(),"Length")
            << G4setw(8) << G4BestUnit(fTrack->GetPosition().z(),"Length")
            << G4setw(8) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
            << G4setw(8) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
            << G4setw(8) << G4BestUnit(fStep->GetStepLength(),"Length")
            << G4setw(8) << G4BestUnit(fTrack->GetTrackLength(),"Length")
            << G4setw(12) << volname << ":" << copyno
            << "   initStep"
            << G4endl;
   }

   G4cout.precision(prec);
}
