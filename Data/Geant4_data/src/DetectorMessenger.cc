#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
: G4UImessenger(), fDetector(Det)
{
  fTestDir = new G4UIdirectory("/testhadr/");
  fTestDir->SetGuidance("Hadr01 detector control.");

  fMatCmd = new G4UIcmdWithAString("/testhadr/TargetMat",this);
  fMatCmd->SetGuidance("Select Material of the Target.");
  fMatCmd->SetParameterName("material",false);
  fMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fRCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/TargetRadius",this);
  fRCmd->SetGuidance("Set radius of the Target.");
  fRCmd->SetParameterName("radius",false);
  fRCmd->SetUnitCategory("Length");
  fRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fLCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/TargetLength",this);
  fLCmd->SetGuidance("Set length of the Target.");
  fLCmd->SetParameterName("length",false);
  fLCmd->SetUnitCategory("Length");
  fLCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

DetectorMessenger::~DetectorMessenger()
{
  delete fMatCmd;
  delete fRCmd;
  delete fLCmd;
  delete fTestDir;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == fMatCmd ) {
    fDetector->SetTargetMaterial(newValue);
  } else if( command == fRCmd ) {
    fDetector->SetTargetRadius(fRCmd->GetNewDoubleValue(newValue));
  } else if( command == fLCmd ) {
    fDetector->SetTargetLength(fLCmd->GetNewDoubleValue(newValue));
  }
}