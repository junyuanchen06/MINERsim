#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UIcmdWithABool.hh>
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"


class PrimaryGeneratorAction;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
    ~PrimaryGeneratorMessenger();

    void SetNewValue(G4UIcommand*, G4String);

  private:
    PrimaryGeneratorAction* primaryGeneratorAction;
    G4UIdirectory* cryDir;
    G4UIcmdWithAString* cryInputCmd;
    G4UIcmdWithoutParameter* cryUpdateCmd;
    G4UIcmdWithABool* cryUseCmd;
    G4UIcmdWithADoubleAndUnit* cryDistanceZCmd;
    std::string* messengerInput;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
