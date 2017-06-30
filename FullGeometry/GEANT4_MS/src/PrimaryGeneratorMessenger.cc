#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* gun)
    : primaryGeneratorAction(gun)
{
  cryDir = new G4UIdirectory("/cry/");
  cryDir->SetGuidance("CRY initialization");

  cryInputCmd = new G4UIcmdWithAString("/cry/input", this);
  cryInputCmd->SetGuidance("CRY input lines");
  cryInputCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  cryUpdateCmd = new G4UIcmdWithoutParameter("/cry/update", this);
  cryUpdateCmd->SetGuidance("Update CRY definition.");
  cryUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\".");
  cryUpdateCmd->SetGuidance("if you changed the CRY definition.");
  cryUpdateCmd->AvailableForStates(G4State_Idle);

  cryUseCmd = new G4UIcmdWithABool("/cry/useCry", this);
  cryUseCmd->SetGuidance("Enable/Disable CRY, 1 for using CRY, 0 for using general particle source.");

  cryDistanceZCmd = new G4UIcmdWithADoubleAndUnit("/cry/distanceZ", this);
  cryDistanceZCmd->SetGuidance("Set Z distance of the plane from which the cosmic rays are generated.");
  cryDistanceZCmd->SetParameterName("cryDistanceZCmd", false);
  cryDistanceZCmd->SetDefaultUnit("m");

  messengerInput = new std::string;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete cryDir;
  delete cryInputCmd;
  delete cryUpdateCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == cryInputCmd) {
    (*messengerInput).append(newValue);
    (*messengerInput).append(" ");
  }
  else if (command == cryUpdateCmd) {
    primaryGeneratorAction->updateCry(messengerInput);
    *messengerInput = "";
  }
  else if (command == cryUseCmd) {
    primaryGeneratorAction->setUseCry(cryUseCmd->GetNewBoolValue(newValue));
  }
  else if (command == cryDistanceZCmd) {
    primaryGeneratorAction->setCryDistanceZ(cryDistanceZCmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
