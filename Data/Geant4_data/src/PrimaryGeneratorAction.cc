#include "PrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
  fParticleGun  = new G4ParticleGun(1);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,50.*CLHEP::m));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}