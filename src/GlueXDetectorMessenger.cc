//
// GlueXDetectorMessenger class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXDetectorMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "GlueXDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4GeometryManager.hh"
#include "G4RunManager.hh"
#include "GlueXPhysicsList.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXBeamConversionProcess.hh"

GlueXDetectorMessenger::GlueXDetectorMessenger(GlueXDetectorConstruction* myDet)
:myDetector(myDet)
{ 
  hdgeant4Dir = new G4UIdirectory("/hdgeant4/");
  hdgeant4Dir->SetGuidance("UI commands specific to hdgeant4");
  
  detDir = new G4UIdirectory("/hdgeant4/det/");
  detDir->SetGuidance("detector control");
  
  FieldCmd = new G4UIcmdWithADoubleAndUnit("/hdgeant4/det/setField",this);
  FieldCmd->SetGuidance("Define magnetic field.");
  FieldCmd->SetGuidance("Magnetic field will be in z direction.");
  FieldCmd->SetParameterName("Bz",false);
  FieldCmd->SetUnitCategory("Magnetic flux density");
  FieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
        
  StepMaxCmd = new G4UIcmdWithADoubleAndUnit("/hdgeant4/det/stepMax",this);
  StepMaxCmd->SetGuidance("Define a step max");
  StepMaxCmd->SetParameterName("stepMax",true,true);
  StepMaxCmd->SetUnitCategory("Length");
  StepMaxCmd->AvailableForStates(G4State_Idle);

  OpenGeomCmd = new G4UIcmdWithoutParameter("/hdgeant4/openGeometry",this);
  OpenGeomCmd->SetGuidance("Re-open the geometry.");
  OpenGeomCmd->SetGuidance("This command is only needed "
                           "if the geometry has already been closed.");
  OpenGeomCmd->AvailableForStates(G4State_Idle);

  CloseGeomCmd = new G4UIcmdWithoutParameter("/hdgeant4/closeGeometry",this);
  CloseGeomCmd->SetGuidance("Explicitly close the geometry.");
  CloseGeomCmd->AvailableForStates(G4State_Idle);

  RadiatorAnglesCmd = new G4UIcmdWith3Vector("/hdgeant4/diamondAngles",this);
  RadiatorAnglesCmd->SetGuidance("Control the diamond radiator orientation.");
  RadiatorAnglesCmd->SetGuidance("Arguments are theta_x, theta_y, theta_z");
  RadiatorAnglesCmd->SetGuidance(" in radians, applied in the order theta_x,");
  RadiatorAnglesCmd->SetGuidance(" followed by theta_y, followed by theta_z.");
  RadiatorAnglesCmd->SetParameterName("theta_x","theta_y","theta_z",true,true);
  RadiatorAnglesCmd->AvailableForStates(G4State_Idle);

  TargetNuclPolCmd = new G4UIcmdWith3Vector("/hdgeant4/targetNuclearPolarization",this);
  TargetNuclPolCmd->SetGuidance("Control the primary target nuclear polarization.");
  TargetNuclPolCmd->SetGuidance("Arguments are pol_x, pol_y, pol_z");
  TargetNuclPolCmd->SetGuidance(" in fractions, normalized to a unit vector");
  TargetNuclPolCmd->SetGuidance(" for 100% polarization.");
  TargetNuclPolCmd->SetParameterName("pol_x","pol_y","pol_z",true,true);
  TargetNuclPolCmd->AvailableForStates(G4State_Idle);

  TargetElecPolCmd = new G4UIcmdWith3Vector("/hdgeant4/targetElectronPolarization",this);
  TargetElecPolCmd->SetGuidance("Control the primary target electron polarization.");
  TargetElecPolCmd->SetGuidance("Arguments are pol_x, pol_y, pol_z");
  TargetElecPolCmd->SetGuidance(" in fractions, normalized to a unit vector");
  TargetElecPolCmd->SetGuidance(" for 100% polarization.");
  TargetElecPolCmd->SetParameterName("pol_x","pol_y","pol_z",true,true);
  TargetElecPolCmd->AvailableForStates(G4State_Idle);
}

GlueXDetectorMessenger::~GlueXDetectorMessenger()
{
  delete FieldCmd;
  delete StepMaxCmd;
  delete detDir;
  delete hdgeant4Dir;
  delete OpenGeomCmd;
  delete RadiatorAnglesCmd;
  delete TargetNuclPolCmd;
  delete TargetElecPolCmd;
}

void GlueXDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if (command == FieldCmd) {
    myDetector->SetUniformField(FieldCmd->GetNewDoubleValue(newValue)/tesla);
  }
  else if (command == StepMaxCmd) {
    myDetector->SetMaxStep(StepMaxCmd->GetNewDoubleValue(newValue));
    std::cout << "MaxStep is now " << StepMaxCmd->GetCurrentValue() << std::endl;
  }
  else if (command == OpenGeomCmd) {
    G4GeometryManager::GetInstance()->OpenGeometry();
  }
  else if (command == CloseGeomCmd) {
    G4GeometryManager::GetInstance()->CloseGeometry();
  }
  else if (command == RadiatorAnglesCmd) {
    G4ThreeVector angles(RadiatorAnglesCmd->GetNew3VectorValue(newValue));
    std::map<CobremsGeneration*, int>::iterator igen;
    int assigned(0);
    for (igen = CobremsGeneration::CobremsGenerators.begin();
         igen != CobremsGeneration::CobremsGenerators.end();
         igen++)
    {
       igen->first->setTargetOrientation(angles[0], angles[1], angles[2]);
       std::cout << "radiator angles are now " 
                 << RadiatorAnglesCmd->GetCurrentValue() << std::endl;
       assigned += igen->second;
    }
    if (!assigned) {
       std::cerr << "CobremsGeneration has not been started yet, you must start"
                 << " up the photon beam generator before attempting to set"
                 << " the radiator orientation." << std::endl;
    }
  }
  else if (command == TargetNuclPolCmd) {
    G4ThreeVector polar(RadiatorAnglesCmd->GetNew3VectorValue(newValue));
    const GlueXPhysicsList *phy = dynamic_cast<const GlueXPhysicsList*>
       (G4RunManager::GetRunManager()->GetUserPhysicsList());
    GlueXBeamConversionProcess *con = (GlueXBeamConversionProcess*)
                                      phy->getBeamConversionProcess();
    if (con != 0 && polar.mag() <= 1) {
       double nucl[3];
       double elec[3];
       con->getTargetPolarization(nucl, elec);
       nucl[0] = polar[0];
       nucl[1] = polar[1];
       nucl[2] = polar[2];
       con->setTargetPolarization(nucl, elec);
    }
    else {
       std::cerr << "polarization vector normalization >1, assignment failed." 
                 << std::endl;
    }
    std::cout << "target nuclear polarization vector is now " 
              << TargetNuclPolCmd->GetCurrentValue() << std::endl;
  }
  else if (command == TargetElecPolCmd) {
    G4ThreeVector polar(RadiatorAnglesCmd->GetNew3VectorValue(newValue));
    const GlueXPhysicsList *phy = dynamic_cast<const GlueXPhysicsList*>
       (G4RunManager::GetRunManager()->GetUserPhysicsList());
    GlueXBeamConversionProcess *con = (GlueXBeamConversionProcess*)
                                      phy->getBeamConversionProcess();
    if (con != 0 && polar.mag() <= 1) {
       double nucl[3];
       double elec[3];
       con->getTargetPolarization(nucl, elec);
       elec[0] = polar[0];
       elec[1] = polar[1];
       elec[2] = polar[2];
       con->setTargetPolarization(nucl, elec);
    }
    else {
       std::cerr << "polarization vector normalization >1, assignment failed."
                 << std::endl;
    }
    std::cout << "target electron polarization vector is now " 
              << TargetElecPolCmd->GetCurrentValue() << std::endl;
  }
}

G4String GlueXDetectorMessenger::GetCurrentValue(G4UIcommand* command)
{
  G4String cv;
  if (command == StepMaxCmd) {
    cv = StepMaxCmd->ConvertToString(myDetector->GetMaxStep(cm));
  }
  else if (command == RadiatorAnglesCmd) {
    G4ThreeVector angles;
    std::map<CobremsGeneration*, int>::iterator igen;
    for (igen = CobremsGeneration::CobremsGenerators.begin();
         igen != CobremsGeneration::CobremsGenerators.end();
         igen++)
    {
       angles[0] = igen->first->getTargetThetax();
       angles[1] = igen->first->getTargetThetay();
       angles[2] = igen->first->getTargetThetaz();
       cv = RadiatorAnglesCmd->ConvertToString(angles);
    }
  }
  else if (command == TargetNuclPolCmd) {
    const GlueXPhysicsList *phy = dynamic_cast<const GlueXPhysicsList*>
                   (G4RunManager::GetRunManager()->GetUserPhysicsList());
    if (phy != 0) {
      GlueXBeamConversionProcess *con = (GlueXBeamConversionProcess*)
                                        phy->getBeamConversionProcess();
      double nucl[3];
      double elec[3];
      con->getTargetPolarization(nucl, elec);
      cv = TargetNuclPolCmd->ConvertToString(G4ThreeVector(nucl[0], nucl[1], nucl[2]));
    }
  }
  else if (command == TargetElecPolCmd) {
    const GlueXPhysicsList *phy = dynamic_cast<const GlueXPhysicsList*>
                   (G4RunManager::GetRunManager()->GetUserPhysicsList());
    if (phy != 0) {
      GlueXBeamConversionProcess *con = (GlueXBeamConversionProcess*)
                                        phy->getBeamConversionProcess();
      double nucl[3];
      double elec[3];
      con->getTargetPolarization(nucl, elec);
      cv = TargetElecPolCmd->ConvertToString(G4ThreeVector(elec[0], elec[1], elec[2]));
    }
  }
  return cv;
}
