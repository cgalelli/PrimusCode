#include "StackingAction.hh"
#include "RunAction.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"

StackingAction::StackingAction(RunAction* runAction)
  : G4UserStackingAction(), fRunAction(runAction)
{}

StackingAction::~StackingAction() {}

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track* aTrack) {
    
    // 1. Filter: We only care about Secondaries (ParentID > 0)
    if (aTrack->GetParentID() == 0) return fUrgent;

    // 2. Get Particle Properties
    const G4ParticleDefinition* particle = aTrack->GetDefinition();
    G4String name = particle->GetParticleName();
    G4double energy = aTrack->GetKineticEnergy();
    
    // 3. Filter: Ignore electrons (they are numerous and not nuclear recoils)
    if (name == "e-") return fKill; // Kill electrons to speed up simulation

    // 4. Filter: Threshold Energy & Geometry
    // Your previous logic: energy > 1eV (approx) and inside Z range
    G4double minEnergy = 1.0 * eV; 
    G4double zPos = aTrack->GetPosition().z();

    if (energy > minEnergy && std::abs(zPos) < 5000*mm) {
        
        // 5. Filter: Select Nuclei (Charge > 1) OR Neutrons
        // Also ensure they are direct daughters of the primary (ParentID == 1) 
        // if you only want primary interaction recoils.
        if (aTrack->GetParentID() == 1 && (particle->GetPDGCharge() > 1 || name == "neutron" )) {
            
            std::ofstream* out = fRunAction->GetOutputFile();
            if (out && out->is_open()) {
                *out << name << "\t" 
                     << particle->GetPDGMass()/MeV << "\t" 
                     << energy/MeV << "\t"
                     << aTrack->GetParentID() << "\t"
                     << zPos/mm << "\t"
                     << aTrack->GetTrackID() << G4endl;
            }
        }
    }

    // Standard behavior: Track everything else urgently (unless you want to kill everything)
    return fUrgent;
}