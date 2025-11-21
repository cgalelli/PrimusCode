#include "EventAction.hh"
#include "G4Event.hh"
#include "EventActionMessenger.hh"

#include "G4UImanager.hh"
#include "G4ios.hh"

EventAction::EventAction():
  G4UserEventAction(),
  fEventMessenger(0),
  fUI(0),
  fSelectedEvents(),
  fPrintModulo(100),
  fSelected(0),
  fDebugStarted(false)
{
  fEventMessenger = new EventActionMessenger(this);
  fUI = G4UImanager::GetUIpointer();
}

EventAction::~EventAction()
{
  delete fEventMessenger;
}

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  // New event
  G4int nEvt = evt->GetEventID();

  if(fSelected>0) {
    for(G4int i=0; i<fSelected; ++i) {
      if(nEvt == fSelectedEvents[i]) {
        fUI->ApplyCommand("/random/saveThisEvent");
        fUI->ApplyCommand("/tracking/verbose  2");
        fDebugStarted = true;
        break;
      }
    }
  }

  // Initialize user actions
  if(G4int(nEvt/fPrintModulo)*fPrintModulo == nEvt) {
    G4cout << "EventAction: Event # "
           << nEvt << " started" << G4endl;
  }
}

void EventAction::EndOfEventAction(const G4Event* evt)
{
  if(fDebugStarted) {
    fUI->ApplyCommand("/tracking/verbose  0");
    fDebugStarted = false;
    G4cout << "EventAction: Event ended" << G4endl;
  }
}