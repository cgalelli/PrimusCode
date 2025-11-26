#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fLogicTarget(nullptr), fLogicWorld(nullptr), 
  fTargetMaterial(nullptr), fWorldMaterial(nullptr),
  fTargetLength(200*m), fTargetRadius(10*m), // Default values
  fMessenger(nullptr)
{
  fMessenger = new DetectorMessenger(this);
  DefineMaterials();
}

DetectorConstruction::~DetectorConstruction()
{
  delete fMessenger;
}

void DetectorConstruction::DefineMaterials()
{
  G4NistManager* nist = G4NistManager::Instance();
  
  // --- Elements ---
  G4Element* H  = nist->FindOrBuildElement("H");
  G4Element* O  = nist->FindOrBuildElement("O");
  G4Element* Na = nist->FindOrBuildElement("Na");
  G4Element* Mg = nist->FindOrBuildElement("Mg");
  G4Element* Si = nist->FindOrBuildElement("Si");
  G4Element* Cl = nist->FindOrBuildElement("Cl");
  G4Element* Fe = nist->FindOrBuildElement("Fe");
  G4Element* Ca = nist->FindOrBuildElement("Ca");
  G4Element* C  = nist->FindOrBuildElement("C");

  // --- Default World Material (Vacuum) ---
  fWorldMaterial = nist->FindOrBuildMaterial("G4_Galactic");

  // --- HALITE (NaCl) ---
  // Density ~2.16 g/cm3
  G4Material* Halite = new G4Material("Halite", 2.16*CLHEP::g/CLHEP::cm3, 2);
  Halite->AddElement(Na, 0.5);
  Halite->AddElement(Cl, 0.5);

  // --- OLIVINE ((Mg,Fe)2SiO4) ---
  // Typical mantle xenolith composition: approx (Mg0.9 Fe0.1)2 SiO4
  // We define a generic Olivine with density ~3.32 g/cm3
  G4Material* Olivine = new G4Material("Olivine", 3.32*CLHEP::g/CLHEP::cm3, 4);
  Olivine->AddElement(Mg, 0.257); // Mg_1.8
  Olivine->AddElement(Fe, 0.029); // Fe_0.2
  Olivine->AddElement(Si, 0.571);
  Olivine->AddElement(O, 0.143);

  G4Material* StdRock = new G4Material("StdRock",2.65*CLHEP::g/CLHEP::cm3,4, kStateSolid );
  StdRock->AddElement(O,  52.*perCent);
  StdRock->AddElement(Ca, 27.*perCent);
  StdRock->AddElement(C,  12.*perCent);
  StdRock->AddElement(Mg,  9.*perCent);

  G4Material* Water = new G4Material("H2O",1.*CLHEP::g/CLHEP::cm3,2);
  Water->AddElement(H, 2);
  Water->AddElement(O, 1);

  // --- Standard Materials (Optional) ---
  nist->FindOrBuildMaterial("G4_SiO2");
  nist->FindOrBuildMaterial("G4_AIR");

  // Set default target (can be overridden by macro)
  fTargetMaterial = Halite;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // --- World ---
  // Make world larger than target
  G4double worldSize = 2.0 * std::max(fTargetLength, fTargetRadius) * 1.2; 
  
  G4Box* solidWorld = new G4Box("World", worldSize/2, worldSize/2, worldSize/2);
  fLogicWorld = new G4LogicalVolume(solidWorld, fWorldMaterial, "World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0);

  // --- Target (Cylinder) ---
  G4Tubs* solidTarget = new G4Tubs("Target", 0., fTargetRadius, fTargetLength/2, 0., 360.*deg);
  fLogicTarget = new G4LogicalVolume(solidTarget, fTargetMaterial, "Target");
  new G4PVPlacement(0, G4ThreeVector(), fLogicTarget, "Target", fLogicWorld, false, 0);

  return physWorld;
}

void DetectorConstruction::SetTargetMaterial(const G4String& matName)
{
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(matName);
  if (pttoMaterial && pttoMaterial != fTargetMaterial) {
     fTargetMaterial = pttoMaterial;
     if (fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
     G4RunManager::GetRunManager()->PhysicsHasBeenModified();
     G4cout << "### Target material set to " << matName << G4endl;
  } else if (pttoMaterial == fTargetMaterial){
     G4cout << "### Target material already is set to " << matName << G4endl;
  } else {
     G4cout << "### WARNING: Material " << matName << " not found!" << G4endl;
  }
}

void DetectorConstruction::SetTargetRadius(G4double val)
{
  fTargetRadius = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetTargetLength(G4double val)
{
  fTargetLength = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}