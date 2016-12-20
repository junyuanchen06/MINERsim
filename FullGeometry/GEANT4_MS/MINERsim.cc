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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
*/
#include "G4MTRunManager.hh"
#include "G4RunManager.hh"

#include "G4UImanager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "MINERMaterials.hh"
#include "GeometryConstruction.hh"
#include "Shielding.hh"
#include "ActionInitialization.hh"
#include "RootIO.hh"

// stuff for importance sampling
#include "G4GeometryManager.hh"
#include "ImportanceGeometryConstruction.hh"
#include "MonitorGeometryConstruction.hh"
#include "G4ImportanceBiasing.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4GeometrySampler.hh"
#include "G4IStore.hh"
#include "G4VWeightWindowStore.hh"
#include "G4WeightWindowAlgorithm.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

int main(int argc,char** argv) {


  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  MINERMaterials::Instance();
  GeometryConstruction* detector = new GeometryConstruction;
  runManager->SetUserInitialization(detector);

  ImportanceGeometryConstruction *pdet = new ImportanceGeometryConstruction("ParallelBiasingWorld");
  detector->RegisterParallelWorld(pdet);

  MonitorGeometryConstruction *mondet = new MonitorGeometryConstruction("ParallelMonitoringWorld");
  detector->RegisterParallelWorld(mondet);

  G4GeometrySampler pgs(pdet->GetWorldVolume(),"gamma");
  pgs.SetParallel(true);

  G4GeometrySampler pgs2(pdet->GetWorldVolume(),"neutron");
  pgs2.SetParallel(true);

  G4VModularPhysicsList *physicsList = new Shielding;
  physicsList->RegisterPhysics(new G4ParallelWorldPhysics("ParallelBiasingWorld"));
  physicsList->RegisterPhysics(new G4ImportanceBiasing(&pgs,"ParallelBiasingWorld"));
  physicsList->RegisterPhysics(new G4ImportanceBiasing(&pgs2,"ParallelBiasingWorld"));
  physicsList->RegisterPhysics(new G4ParallelWorldPhysics("ParallelMonitoringWorld"));


  /*
  // add thermal neutron model
  G4HadronElasticProcess* theNeutronElasticProcess = new G4HadronElasticProcess;
  // Cross Section Data set
  G4NeutronHPElasticData* theHPElasticData = new G4NeutronHPElasticData;
  theNeutronElasticProcess->AddDataSet(theHPElasticData);
  G4NeutronHPThermalScatteringData* theHPThermalScatteringData = new G4NeutronHPThermalScatteringData;
  theNeutronElasticProcess->AddDataSet(theHPThermalScatteringData);
  // Models
  G4NeutronHPElastic* theNeutronElasticModel = new G4NeutronHPElastic;
  theNeutronElasticModel->SetMinEnergy(4.0*eV);
  theNeutronElasticProcess->RegisterMe(theNeutronElasticModel);
  G4NeutronHPThermalScattering* theNeutronThermalElasticModel = new G4NeutronHPThermalScattering;
  theNeutronThermalElasticModel->SetMaxEnergy(4.0*eV);
  theNeutronElasticProcess->RegisterMe(theNeutronThermalElasticModel);
  // Apply Processes to Process Manager of Neutron
  G4ProcessManager* pmanager = G4Neutron::Neutron()->GetProcessManager();
  pmanager->AddDiscreteProcess(theNeutronElasticProcess);
  */

  runManager->SetUserInitialization(physicsList);
  runManager->SetUserInitialization(new ActionInitialization);
  RootIO::GetInstance();

  runManager->Initialize();

  pdet->CreateImportanceStore();

  pgs2.PrepareImportanceSampling(G4IStore::GetInstance(pdet->GetName()),0);
  pgs2.Configure();


#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
    
  // get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);    
    }
  else          // interactive mode : define UI session
    {
#ifdef G4UI_USE
     G4UIExecutive* ui = new G4UIExecutive(argc, argv, "tcsh");
#ifdef G4VIS_USE
     UImanager->ApplyCommand("/control/execute vis.mac");
#endif             
     ui->SessionStart();
     delete ui;
#endif
     
#ifdef G4VIS_USE
     delete visManager;
#endif    
    }


  // job termination
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  pgs.ClearSampling();
  pgs2.ClearSampling();
  delete runManager;
  
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 
