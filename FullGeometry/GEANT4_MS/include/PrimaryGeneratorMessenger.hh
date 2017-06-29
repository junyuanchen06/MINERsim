#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"


class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;

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
    std::string* messengerInput;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

