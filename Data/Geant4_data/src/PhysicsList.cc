#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

// Physics Constructors
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmExtraPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh" 
#include "G4StoppingPhysics.hh"
#include "G4IonPhysics.hh"

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  SetVerboseLevel(0);
  
  RegisterPhysics(new G4EmStandardPhysics_option4(0));
  RegisterPhysics(new G4EmExtraPhysics());

  RegisterPhysics(new G4DecayPhysics(0));
  RegisterPhysics(new G4RadioactiveDecayPhysics(0));

  RegisterPhysics(new G4HadronElasticPhysics(0));

  RegisterPhysics(new G4HadronPhysicsQGSP_BIC_HP(0)); 

  RegisterPhysics(new G4StoppingPhysics(0));

  RegisterPhysics(new G4IonPhysics(0));
}

PhysicsList::~PhysicsList()
{}

void PhysicsList::ConstructParticle()
{
  G4VModularPhysicsList::ConstructParticle();
}

void PhysicsList::ConstructProcess()
{
  G4VModularPhysicsList::ConstructProcess();
}

void PhysicsList::AddPhysicsList(const G4String& name)
{
  G4cout << "PhysicsList::AddPhysicsList: " << name << " (Ignored)" << G4endl;
}