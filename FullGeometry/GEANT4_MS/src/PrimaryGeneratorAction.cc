//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: PrimaryGeneratorAction.cc 67272 2013-02-13 12:05:47Z ihrivnac $
//

#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/Random.h>
#include <RunAction.hh>
#include "PrimaryGeneratorAction.hh"
#include "RNGWrapper.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(RunAction* run)
    : G4VUserPrimaryGeneratorAction(), generalParticleSource(0), useCry(false), cryDistanceZ(0 * m), runAction(run)
{
  generalParticleSource = new G4GeneralParticleSource();
  particleGun = new G4ParticleGun();
  particleTable = G4ParticleTable::GetParticleTable();
  cryParticles = new std::vector<CRYParticle*>;
  primaryGeneratorMessenger = new PrimaryGeneratorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete generalParticleSource;
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::updateCry(std::string *messengerInput)
{
  CRYSetup *crySetup = new CRYSetup(*messengerInput, getenv(CRY_DATA_PATH_ENVIRONMENT_VARIABLE));
  cryGenerator = new CRYGenerator(crySetup);

  RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(),&CLHEP::HepRandomEngine::flat);
  crySetup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::setUseCry(G4bool use)
{
  useCry = use;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::setCryDistanceZ(G4double distance)
{
  cryDistanceZ = distance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(useCry){
    G4String particleName;
    cryParticles->clear();
    cryGenerator->genEvent(cryParticles);
    runAction->setElapsedTime(cryGenerator->timeSimulated());

    for (unsigned j=0; j<cryParticles->size(); j++) {
      particleName = CRYUtils::partName((*cryParticles)[j]->id());
      particleGun->SetParticleDefinition(particleTable->FindParticle((*cryParticles)[j]->PDGid()));
      particleGun->SetParticleEnergy((*cryParticles)[j]->ke()*MeV);
      particleGun->SetParticlePosition(G4ThreeVector((*cryParticles)[j]->x()*m, (*cryParticles)[j]->y()*m, (*cryParticles)[j]->z()*m + cryDistanceZ));
      particleGun->SetParticleMomentumDirection(G4ThreeVector((*cryParticles)[j]->u(), (*cryParticles)[j]->v(), (*cryParticles)[j]->w()));
      particleGun->SetParticleTime((*cryParticles)[j]->t());
      particleGun->GeneratePrimaryVertex(anEvent);
      delete (*cryParticles)[j];
    }
  }
  else {
    generalParticleSource->GeneratePrimaryVertex(anEvent);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
