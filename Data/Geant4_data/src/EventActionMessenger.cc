#include "EventActionMessenger.hh"
#include "EventAction.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

EventActionMessenger::EventActionMessenger(EventAction* evAct)
  : G4UImessenger(), fEventAction(evAct),
    fPrintCmd(0), fCmd(0)
{   
  fPrintCmd = new G4UIcmdWithAnInteger("/testhadr/PrintModulo",this);
  fPrintCmd->SetGuidance("Print events modulo n");
  fPrintCmd->SetParameterName("EventNb",false);
  fPrintCmd->SetRange("EventNb>0");
  fPrintCmd->AvailableForStates(G4State_PreInit,G4State_Idle);      

  fCmd = new G4UIcmdWithAnInteger("/testhadr/DebugEvent",this);
  fCmd->SetGuidance("D event to debug");
  fCmd->SetParameterName("fNb",false);
  fCmd->SetRange("fNb>0");
  fCmd->AvailableForStates(G4State_PreInit,G4State_Idle);      

}

EventActionMessenger::~EventActionMessenger()
{
  delete fPrintCmd;   
  delete fCmd;
}

void EventActionMessenger::SetNewValue(G4UIcommand* command,
                                       G4String newValue)
{ 
  if(command == fPrintCmd)
    {fEventAction->SetPrintModulo(fPrintCmd->GetNewIntValue(newValue));}
  if(command == fCmd)
    {fEventAction->AddEventToDebug(fCmd->GetNewIntValue(newValue));}
}