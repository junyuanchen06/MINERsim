#include "globals.hh"
#include <sstream>

#include "MonitorGeometryConstruction.hh"
#include "MinerSD.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// For Primitive Scorers
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSNofCollision.hh"
#include "G4PSPopulation.hh"
#include "G4PSTrackCounter.hh"
#include "G4PSTrackLength.hh"

// for importance biasing
#include "G4IStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MonitorGeometryConstruction::
MonitorGeometryConstruction(G4String worldName) 
:G4VUserParallelWorld(worldName)
{
  //  Construct();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MonitorGeometryConstruction::~MonitorGeometryConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MonitorGeometryConstruction::Construct()
{  
  G4cout << " constructing parallel world " << G4endl;

  G4Material* dummyMat  = 0;

  //GetWorld methods create a clone of the mass world to the parallel world (!)
  fGhostWorld = GetWorld();
  G4LogicalVolume* worldLogical = fGhostWorld->GetLogicalVolume();


  G4double in = 2.54*cm;
  G4double worldY = 10*m;
  G4double worldX = 10*m;
  G4double worldZ = 6*m;
  G4double bioShield_thick = (19.625 + 49.6875) *in ; //1.7*m; // this includes SS plating which is why you'll see it subtracted out below
  G4double water_rad = 12.*4.5*in;  // from Sean's talk at TAMU CNS workshop
  G4double floor_to_TC_out = 48*in; // distance from ground floor to center of thermal column, ESTIMATE
  G4double floor_to_TC_pool = 60*in; // distance from pool floor to center of thermal column 
  G4double floor_thick_out = 0.5*m; // reasonable?
  G4double floor_thick_pool = (floor_thick_out-(floor_to_TC_pool-floor_to_TC_out)); // this includes SS plating which is why you'll see it subtracted out below
  G4double thermalC_height_big = 1.24*m; // larger part of thermal column, square face side
  G4double thermalC_height_small = 1.0*m; // smaller part of thermal column, square face side
  G4double thermalC_length_big = 49.6875*in;
  G4double thermalC_length_small = (19.625+8.875+25.8125)*in;
  G4double bioShield_flat_width = thermalC_height_big + 1*m; // this is the width of the "flat" part of the bio shield where the thermal column is ESTIMATE
  G4double aluminumPlate_thick = 0.5*in; // THIS IS A GUESS
  G4double SSPlate_thick = 0.5*in; // THIS IS A GUESS
  G4double poolC_to_TC = water_rad - (25.8125+8.875)*in; // this is the distance from the center of the circular part of the pool to the back of the thermal column
  G4double gb_depth = 26.5625*in; // depth of graphite block
  G4double gb_width = 29.5*in; // width of graphite block ESTIMATE
  G4double gb_height = 29.5*in; // height of graphite block
  G4double largeOverburden_length = 19.625*in; // this is the lenght of the part of the TC with the largest overburden
  G4double iceboxThick = 3./8.*in;  // ???
  G4double iceboxRadius = 6*in;
  G4double iceboxHeight = 18*in;
  G4double iceboxLead_thick = 2*in;
  G4double iceboxPoly_thick = 2*in;
  G4double muVeto_thick = 2*in;
  G4double backScint_thick = 8*in;
  G4double outerLead_thick = 2*in;
  G4double outerPoly_thick = 6*in;
  G4double outerPoly_back_thick = 4*in;
  G4double TCLead_thick = 8*in;
  G4double TCPoly_thick = thermalC_length_small-SSPlate_thick;
  G4double detThick = 1.33*in;
  G4double detSpace = 0.4*in;


  G4double layerThick = 0.001*mm;

  G4Tubs *insideIB1 = new G4Tubs("insideIB1", 0., iceboxRadius - 5*iceboxThick, iceboxHeight/2. - 5*iceboxThick,0,360*deg);
  G4Tubs *insideIB2 = new G4Tubs("insideIB1", 0., iceboxRadius - 5*iceboxThick - layerThick, iceboxHeight/2. - 5*iceboxThick - layerThick,0,360*deg);
  G4SubtractionSolid *insideIB = new G4SubtractionSolid("insideIB",insideIB1,insideIB2);
  G4LogicalVolume *insideIB_log = new G4LogicalVolume(insideIB, dummyMat, "insideIB_log");


  G4Box *FluxMonitor = new G4Box("FluxMonitor",layerThick/2.,thermalC_height_big/2.,thermalC_height_big/2.);

  G4LogicalVolume *atNeutronDet_log = new G4LogicalVolume(FluxMonitor, dummyMat, "atNeutronDet_log");
  G4LogicalVolume *atLead_log = new G4LogicalVolume(FluxMonitor, dummyMat, "atLead_log");
  G4LogicalVolume *atStartTC_log = new G4LogicalVolume(FluxMonitor, dummyMat, "atStartTC_log");
  G4LogicalVolume *atPoly1_log = new G4LogicalVolume(FluxMonitor, dummyMat, "atPoly1_log");
  G4LogicalVolume *atPoly2_log = new G4LogicalVolume(FluxMonitor, dummyMat, "atPoly2_log");
  G4LogicalVolume *atPoly3_log = new G4LogicalVolume(FluxMonitor, dummyMat, "atPoly3_log");
  G4LogicalVolume *atPoly4_log = new G4LogicalVolume(FluxMonitor, dummyMat, "atPoly4_log");
  G4LogicalVolume *atPoly5_log = new G4LogicalVolume(FluxMonitor, dummyMat, "atPoly5_log");
  G4LogicalVolume *atPoly6_log = new G4LogicalVolume(FluxMonitor, dummyMat, "atPoly6_log");
  G4LogicalVolume *atPoly7_log = new G4LogicalVolume(FluxMonitor, dummyMat, "atPoly7_log");
  G4LogicalVolume *atPoly8_log = new G4LogicalVolume(FluxMonitor, dummyMat, "atPoly8_log");
  G4LogicalVolume *atPoly9_log = new G4LogicalVolume(FluxMonitor, dummyMat, "atPoly9_log");


  G4ThreeVector   posShield(water_rad+bioShield_thick +SSPlate_thick  - thermalC_length_big + iceboxRadius+iceboxPoly_thick+iceboxLead_thick+TCLead_thick , 0, -(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  G4ThreeVector   atNeutronDet_place(posShield.x() - (iceboxRadius+iceboxPoly_thick+iceboxLead_thick),0.,posShield.z());
  G4ThreeVector   atLead_place(posShield.x() - (iceboxRadius+iceboxPoly_thick+iceboxLead_thick) - TCLead_thick,0.,posShield.z());
  G4ThreeVector   atStartTC_place(posShield.x() - (iceboxRadius+iceboxPoly_thick+iceboxLead_thick) - TCLead_thick - TCPoly_thick - SSPlate_thick,0.,posShield.z());

  G4ThreeVector   atPoly1_place(posShield.x() - (iceboxRadius+iceboxPoly_thick+iceboxLead_thick) - TCLead_thick - 0.9*TCPoly_thick - SSPlate_thick,0.,posShield.z());
  G4ThreeVector   atPoly2_place(posShield.x() - (iceboxRadius+iceboxPoly_thick+iceboxLead_thick) - TCLead_thick - 0.8*TCPoly_thick - SSPlate_thick,0.,posShield.z());
  G4ThreeVector   atPoly3_place(posShield.x() - (iceboxRadius+iceboxPoly_thick+iceboxLead_thick) - TCLead_thick - 0.7*TCPoly_thick - SSPlate_thick,0.,posShield.z());
  G4ThreeVector   atPoly4_place(posShield.x() - (iceboxRadius+iceboxPoly_thick+iceboxLead_thick) - TCLead_thick - 0.6*TCPoly_thick - SSPlate_thick,0.,posShield.z());
  G4ThreeVector   atPoly5_place(posShield.x() - (iceboxRadius+iceboxPoly_thick+iceboxLead_thick) - TCLead_thick - 0.5*TCPoly_thick - SSPlate_thick,0.,posShield.z());
  G4ThreeVector   atPoly6_place(posShield.x() - (iceboxRadius+iceboxPoly_thick+iceboxLead_thick) - TCLead_thick - 0.4*TCPoly_thick - SSPlate_thick,0.,posShield.z());
  G4ThreeVector   atPoly7_place(posShield.x() - (iceboxRadius+iceboxPoly_thick+iceboxLead_thick) - TCLead_thick - 0.3*TCPoly_thick - SSPlate_thick,0.,posShield.z());
  G4ThreeVector   atPoly8_place(posShield.x() - (iceboxRadius+iceboxPoly_thick+iceboxLead_thick) - TCLead_thick - 0.2*TCPoly_thick - SSPlate_thick,0.,posShield.z());
  G4ThreeVector   atPoly9_place(posShield.x() - (iceboxRadius+iceboxPoly_thick+iceboxLead_thick) - TCLead_thick - 0.1*TCPoly_thick - SSPlate_thick,0.,posShield.z());


  new G4PVPlacement(0,posShield,insideIB_log,"insideIB_phys",worldLogical,false,0,false);
  new G4PVPlacement(0,atNeutronDet_place,atNeutronDet_log,"atNeutronDet_phys",worldLogical,false,0,false);
  new G4PVPlacement(0,atLead_place,atLead_log,"atLead_phys",worldLogical,false,0,false);
  new G4PVPlacement(0,atStartTC_place,atStartTC_log,"atStartTC_phys",worldLogical,false,0,false);
  new G4PVPlacement(0,atPoly1_place,atPoly1_log,"atPoly1_phys",worldLogical,false,0,false);
  new G4PVPlacement(0,atPoly2_place,atPoly2_log,"atPoly2_phys",worldLogical,false,0,false);
  new G4PVPlacement(0,atPoly3_place,atPoly3_log,"atPoly3_phys",worldLogical,false,0,false);
  new G4PVPlacement(0,atPoly4_place,atPoly4_log,"atPoly4_phys",worldLogical,false,0,false);
  new G4PVPlacement(0,atPoly5_place,atPoly5_log,"atPoly5_phys",worldLogical,false,0,false);
  new G4PVPlacement(0,atPoly6_place,atPoly6_log,"atPoly6_phys",worldLogical,false,0,false);
  new G4PVPlacement(0,atPoly7_place,atPoly7_log,"atPoly7_phys",worldLogical,false,0,false);
  new G4PVPlacement(0,atPoly8_place,atPoly8_log,"atPoly8_phys",worldLogical,false,0,false);
  new G4PVPlacement(0,atPoly9_place,atPoly9_log,"atPoly9_phys",worldLogical,false,0,false);



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume &MonitorGeometryConstruction::
GetWorldVolumeAddress() const{
   return *fGhostWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *MonitorGeometryConstruction::GetWorldVolume() {
  return fGhostWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MonitorGeometryConstruction::SetSensitive(){

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MonitorGeometryConstruction::ConstructSD()
{

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  MinerSD* insideIB_Det = new MinerSD("/MINERsim/insideIB_Det","MS_insideIB_hits");
  SDman->AddNewDetector(insideIB_Det);
  SetSensitiveDetector("insideIB_log",insideIB_Det, true);

  MinerSD* atNeutronDet_Det = new MinerSD("/MINERsim/atNeutronDet_Det","MS_atNeutronDet_hits");
  SDman->AddNewDetector(atNeutronDet_Det);
  SetSensitiveDetector("atNeutronDet_log",atNeutronDet_Det, true);

  MinerSD* atLead_Det = new MinerSD("/MINERsim/atLead_Det","MS_atLead_hits");
  SDman->AddNewDetector(atLead_Det);
  SetSensitiveDetector("atLead_log",atLead_Det, true);

  MinerSD* atStartTC_Det = new MinerSD("/MINERsim/atStartTC_Det","MS_atStartTC_hits");
  SDman->AddNewDetector(atStartTC_Det);
  SetSensitiveDetector("atStartTC_log",atStartTC_Det, true);

  MinerSD* atPoly1_Det = new MinerSD("/MINERsim/atPoly1_Det","MS_atPoly1_hits");
  SDman->AddNewDetector(atPoly1_Det);
  SetSensitiveDetector("atPoly1_log",atPoly1_Det, true);

  MinerSD* atPoly2_Det = new MinerSD("/MINERsim/atPoly2_Det","MS_atPoly2_hits");
  SDman->AddNewDetector(atPoly2_Det);
  SetSensitiveDetector("atPoly2_log",atPoly2_Det, true);

  MinerSD* atPoly3_Det = new MinerSD("/MINERsim/atPoly3_Det","MS_atPoly3_hits");
  SDman->AddNewDetector(atPoly3_Det);
  SetSensitiveDetector("atPoly3_log",atPoly3_Det, true);

  MinerSD* atPoly4_Det = new MinerSD("/MINERsim/atPoly4_Det","MS_atPoly4_hits");
  SDman->AddNewDetector(atPoly4_Det);
  SetSensitiveDetector("atPoly4_log",atPoly4_Det, true);

  MinerSD* atPoly5_Det = new MinerSD("/MINERsim/atPoly5_Det","MS_atPoly5_hits");
  SDman->AddNewDetector(atPoly5_Det);
  SetSensitiveDetector("atPoly5_log",atPoly5_Det, true);

  MinerSD* atPoly6_Det = new MinerSD("/MINERsim/atPoly6_Det","MS_atPoly6_hits");
  SDman->AddNewDetector(atPoly6_Det);
  SetSensitiveDetector("atPoly6_log",atPoly6_Det, true);

  MinerSD* atPoly7_Det = new MinerSD("/MINERsim/atPoly7_Det","MS_atPoly7_hits");
  SDman->AddNewDetector(atPoly7_Det);
  SetSensitiveDetector("atPoly7_log",atPoly7_Det, true);

  MinerSD* atPoly8_Det = new MinerSD("/MINERsim/atPoly8_Det","MS_atPoly8_hits");
  SDman->AddNewDetector(atPoly8_Det);
  SetSensitiveDetector("atPoly8_log",atPoly8_Det, true);

  MinerSD* atPoly9_Det = new MinerSD("/MINERsim/atPoly9_Det","MS_atPoly9_hits");
  SDman->AddNewDetector(atPoly9_Det);
  SetSensitiveDetector("atPoly9_log",atPoly9_Det, true);



}

