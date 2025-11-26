#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "StackingAction.hh"
#include "SteppingAction.hh"

int main(int argc,char** argv) {
  G4RunManager* runManager = new G4RunManager;

  // 1. Geometry
  runManager->SetUserInitialization(new DetectorConstruction());

  // 2. Physics List (Your custom class)
  runManager->SetUserInitialization(new PhysicsList());

  // 3. User Actions
  runManager->SetUserAction(new PrimaryGeneratorAction());
  
  RunAction* runAction = new RunAction();
  runManager->SetUserAction(runAction);
  
  runManager->SetUserAction(new EventAction());
  runManager->SetUserAction(new StackingAction());
  runManager->SetUserAction(new SteppingAction(runAction));
  

  // 4. Visualization & UI
//  G4VisManager* visManager = new G4VisExecutive;
//  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  G4String command = "/control/execute ";
  G4String fileName = argv[1];
  UImanager->ApplyCommand(command+fileName);

//  delete visManager;
  delete runManager;
}