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

#include "EventAction.hh"

#include "TClonesArray.h"


#include "MinerHit.hh"
#include "BaseHit.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"

#include "RootIO.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
: G4UserEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  
  G4SDManager* fSDM = G4SDManager::GetSDMpointer();
  G4HCofThisEvent* HCofEvent = event->GetHCofThisEvent();
  MinerHitsCollection*  vacHits1 = (MinerHitsCollection*)(HCofEvent->GetHC(fSDM->GetCollectionID("ACVacHitsCollection1")));
  MinerHitsCollection*  vacHits2 = (MinerHitsCollection*)(HCofEvent->GetHC(fSDM->GetCollectionID("ACVacHitsCollection2")));
  MinerHitsCollection*  waterHits = (MinerHitsCollection*)(HCofEvent->GetHC(fSDM->GetCollectionID("ACWaterHitsCollection")));
  MinerHitsCollection*  graphiteHits = (MinerHitsCollection*)(HCofEvent->GetHC(fSDM->GetCollectionID("ACGraphiteHitsCollection")));
  MinerHitsCollection*  leadHits = (MinerHitsCollection*)(HCofEvent->GetHC(fSDM->GetCollectionID("ACLeadHitsCollection")));
  MinerHitsCollection*  polyHits = (MinerHitsCollection*)(HCofEvent->GetHC(fSDM->GetCollectionID("ACPolyHitsCollection")));
  MinerHitsCollection*  zipHits = (MinerHitsCollection*)(HCofEvent->GetHC(fSDM->GetCollectionID("ACZipHitsCollection")));
  MinerHitsCollection*  concreteHits = (MinerHitsCollection*)(HCofEvent->GetHC(fSDM->GetCollectionID("ACConcreteHitsCollection")));
  MinerHitsCollection*  scintHits = (MinerHitsCollection*)(HCofEvent->GetHC(fSDM->GetCollectionID("ACScintHitsCollection")));


  RootIO::GetInstance()->AddHits(waterHits,0);
  RootIO::GetInstance()->AddHits(graphiteHits,1);
  RootIO::GetInstance()->AddHits(leadHits,2);
  RootIO::GetInstance()->AddHits(polyHits,3);
  RootIO::GetInstance()->AddHits(vacHits1,4);
  RootIO::GetInstance()->AddHits(vacHits2,5);
  RootIO::GetInstance()->AddHits(zipHits,6);
  RootIO::GetInstance()->AddHits(concreteHits,7);
  RootIO::GetInstance()->AddHits(scintHits,8);

  RootIO::GetInstance()->Write();


}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
