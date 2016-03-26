#include "globals.hh"
#include <sstream>

#include "ImportanceGeometryConstruction.hh"

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

ImportanceGeometryConstruction::
ImportanceGeometryConstruction(G4String worldName) 
:G4VUserParallelWorld(worldName),fLogicalVolumeVector()
{
  //  Construct();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ImportanceGeometryConstruction::~ImportanceGeometryConstruction()
{
  fLogicalVolumeVector.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ImportanceGeometryConstruction::Construct()
{  
  G4cout << " constructing parallel world " << G4endl;

  G4Material* dummyMat  = 0;

  //GetWorld methods create a clone of the mass world to the parallel world (!)
  // via the transportation manager
  fGhostWorld = GetWorld();
  G4cout << " ImportanceGeometryConstruction:: ghostWorldName = " 
         << fGhostWorld->GetName() << G4endl;
  G4LogicalVolume* worldLogical = fGhostWorld->GetLogicalVolume();
  fLogicalVolumeVector.push_back(worldLogical);



  //  fPVolumeStore.AddPVolume(G4GeometryCell(*pWorldVolume, 0));
  fPVolumeStore.AddPVolume(G4GeometryCell(*fGhostWorld, 0));

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

  G4double water_thick = 3.86*m - gb_depth;

  G4double nRegions = 40.;

//  G4double sizeImpRegions = 2.*m;
  G4double sizeImpRegions = 1.*m;

  // this is the start point of our interesting region where we'll do importance sampling
  G4ThreeVector startPos(water_rad+bioShield_thick +SSPlate_thick - (water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big),0, -(worldZ/2.)+floor_thick_out+floor_to_TC_out);


  // define the part of the world region that we don't really care about (basically behind the reactor)
  G4Box * notInteresting  = new G4Box("notInteresting", (worldX/2. + startPos.x())/2. ,worldY/2.,worldZ/2.);
  G4LogicalVolume * notInteresting_log = new G4LogicalVolume(notInteresting,dummyMat,"notInteresting_log",0,0,0);

  G4Box * notInteresting2  = new G4Box("notInteresting2", (worldX - (water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big) - (worldX/2. + startPos.x()))/2.,worldY/2.,worldZ/2.);
  G4LogicalVolume * notInteresting2_log = new G4LogicalVolume(notInteresting2,dummyMat,"notInteresting2_log",0,0,0);



  G4ThreeVector interestingPlace(water_rad+bioShield_thick +SSPlate_thick - (water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/2.,0, -(worldZ/2.)+floor_thick_out+floor_to_TC_out);


  G4Box * Region = new G4Box("Region", (water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/(2.*nRegions), sizeImpRegions, sizeImpRegions);
  G4LogicalVolume * Region_1_log = new G4LogicalVolume(Region,dummyMat,"Region_1_log");
  G4LogicalVolume * Region_2_log = new G4LogicalVolume(Region,dummyMat,"Region_2_log");
  G4LogicalVolume * Region_3_log = new G4LogicalVolume(Region,dummyMat,"Region_3_log");
  G4LogicalVolume * Region_4_log = new G4LogicalVolume(Region,dummyMat,"Region_4_log");
  G4LogicalVolume * Region_5_log = new G4LogicalVolume(Region,dummyMat,"Region_5_log");
  G4LogicalVolume * Region_6_log = new G4LogicalVolume(Region,dummyMat,"Region_6_log");
  G4LogicalVolume * Region_7_log = new G4LogicalVolume(Region,dummyMat,"Region_7_log");
  G4LogicalVolume * Region_8_log = new G4LogicalVolume(Region,dummyMat,"Region_8_log");
  G4LogicalVolume * Region_9_log = new G4LogicalVolume(Region,dummyMat,"Region_9_log");
  G4LogicalVolume * Region_10_log = new G4LogicalVolume(Region,dummyMat,"Region_10_log");
  G4LogicalVolume * Region_11_log = new G4LogicalVolume(Region,dummyMat,"Region_11_log");
  G4LogicalVolume * Region_12_log = new G4LogicalVolume(Region,dummyMat,"Region_12_log");
  G4LogicalVolume * Region_13_log = new G4LogicalVolume(Region,dummyMat,"Region_13_log");
  G4LogicalVolume * Region_14_log = new G4LogicalVolume(Region,dummyMat,"Region_14_log");
  G4LogicalVolume * Region_15_log = new G4LogicalVolume(Region,dummyMat,"Region_15_log");
  G4LogicalVolume * Region_16_log = new G4LogicalVolume(Region,dummyMat,"Region_16_log");
  G4LogicalVolume * Region_17_log = new G4LogicalVolume(Region,dummyMat,"Region_17_log");
  G4LogicalVolume * Region_18_log = new G4LogicalVolume(Region,dummyMat,"Region_18_log");
  G4LogicalVolume * Region_19_log = new G4LogicalVolume(Region,dummyMat,"Region_19_log");
  G4LogicalVolume * Region_20_log = new G4LogicalVolume(Region,dummyMat,"Region_20_log");
  G4LogicalVolume * Region_21_log = new G4LogicalVolume(Region,dummyMat,"Region_21_log");
  G4LogicalVolume * Region_22_log = new G4LogicalVolume(Region,dummyMat,"Region_22_log");
  G4LogicalVolume * Region_23_log = new G4LogicalVolume(Region,dummyMat,"Region_23_log");
  G4LogicalVolume * Region_24_log = new G4LogicalVolume(Region,dummyMat,"Region_24_log");
  G4LogicalVolume * Region_25_log = new G4LogicalVolume(Region,dummyMat,"Region_25_log");
  G4LogicalVolume * Region_26_log = new G4LogicalVolume(Region,dummyMat,"Region_26_log");
  G4LogicalVolume * Region_27_log = new G4LogicalVolume(Region,dummyMat,"Region_27_log");
  G4LogicalVolume * Region_28_log = new G4LogicalVolume(Region,dummyMat,"Region_28_log");
  G4LogicalVolume * Region_29_log = new G4LogicalVolume(Region,dummyMat,"Region_29_log");
  G4LogicalVolume * Region_30_log = new G4LogicalVolume(Region,dummyMat,"Region_30_log");
  G4LogicalVolume * Region_31_log = new G4LogicalVolume(Region,dummyMat,"Region_31_log");
  G4LogicalVolume * Region_32_log = new G4LogicalVolume(Region,dummyMat,"Region_32_log");
  G4LogicalVolume * Region_33_log = new G4LogicalVolume(Region,dummyMat,"Region_33_log");
  G4LogicalVolume * Region_34_log = new G4LogicalVolume(Region,dummyMat,"Region_34_log");
  G4LogicalVolume * Region_35_log = new G4LogicalVolume(Region,dummyMat,"Region_35_log");
  G4LogicalVolume * Region_36_log = new G4LogicalVolume(Region,dummyMat,"Region_36_log");
  G4LogicalVolume * Region_37_log = new G4LogicalVolume(Region,dummyMat,"Region_37_log");
  G4LogicalVolume * Region_38_log = new G4LogicalVolume(Region,dummyMat,"Region_38_log");
  G4LogicalVolume * Region_39_log = new G4LogicalVolume(Region,dummyMat,"Region_39_log");
  G4LogicalVolume * Region_40_log = new G4LogicalVolume(Region,dummyMat,"Region_40_log");



  fLogicalVolumeVector.push_back(notInteresting_log);
  fLogicalVolumeVector.push_back(Region_1_log);
  fLogicalVolumeVector.push_back(Region_2_log);
  fLogicalVolumeVector.push_back(Region_3_log);
  fLogicalVolumeVector.push_back(Region_4_log);
  fLogicalVolumeVector.push_back(Region_5_log);
  fLogicalVolumeVector.push_back(Region_6_log);
  fLogicalVolumeVector.push_back(Region_7_log);
  fLogicalVolumeVector.push_back(Region_8_log);
  fLogicalVolumeVector.push_back(Region_9_log);
  fLogicalVolumeVector.push_back(Region_10_log);
  fLogicalVolumeVector.push_back(Region_11_log);
  fLogicalVolumeVector.push_back(Region_12_log);
  fLogicalVolumeVector.push_back(Region_13_log);
  fLogicalVolumeVector.push_back(Region_14_log);
  fLogicalVolumeVector.push_back(Region_15_log);
  fLogicalVolumeVector.push_back(Region_16_log);
  fLogicalVolumeVector.push_back(Region_17_log);
  fLogicalVolumeVector.push_back(Region_18_log);
  fLogicalVolumeVector.push_back(Region_19_log);
  fLogicalVolumeVector.push_back(Region_20_log);
  fLogicalVolumeVector.push_back(Region_21_log);
  fLogicalVolumeVector.push_back(Region_22_log);
  fLogicalVolumeVector.push_back(Region_23_log);
  fLogicalVolumeVector.push_back(Region_24_log);
  fLogicalVolumeVector.push_back(Region_25_log);
  fLogicalVolumeVector.push_back(Region_26_log);
  fLogicalVolumeVector.push_back(Region_27_log);
  fLogicalVolumeVector.push_back(Region_28_log);
  fLogicalVolumeVector.push_back(Region_29_log);
  fLogicalVolumeVector.push_back(Region_30_log);
  fLogicalVolumeVector.push_back(Region_31_log);
  fLogicalVolumeVector.push_back(Region_32_log);
  fLogicalVolumeVector.push_back(Region_33_log);
  fLogicalVolumeVector.push_back(Region_34_log);
  fLogicalVolumeVector.push_back(Region_35_log);
  fLogicalVolumeVector.push_back(Region_36_log);
  fLogicalVolumeVector.push_back(Region_37_log);
  fLogicalVolumeVector.push_back(Region_38_log);
  fLogicalVolumeVector.push_back(Region_39_log);
  fLogicalVolumeVector.push_back(Region_40_log);


  // Outer Regions
  G4Box * Out_RegionI = new G4Box("Out_RegionI", (water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/(2.*nRegions), worldY/2., worldZ/2.);
  G4ThreeVector outCut(0.,startPos.y(), startPos.z() );
  G4SubtractionSolid *Out_Region = new G4SubtractionSolid("Out_Region",Out_RegionI,Region,0,outCut);

  G4LogicalVolume * Out_Region_1_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_1_log");
  G4LogicalVolume * Out_Region_2_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_2_log");
  G4LogicalVolume * Out_Region_3_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_3_log");
  G4LogicalVolume * Out_Region_4_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_4_log");
  G4LogicalVolume * Out_Region_5_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_5_log");
  G4LogicalVolume * Out_Region_6_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_6_log");
  G4LogicalVolume * Out_Region_7_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_7_log");
  G4LogicalVolume * Out_Region_8_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_8_log");
  G4LogicalVolume * Out_Region_9_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_9_log");
  G4LogicalVolume * Out_Region_10_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_10_log");
  G4LogicalVolume * Out_Region_11_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_11_log");
  G4LogicalVolume * Out_Region_12_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_12_log");
  G4LogicalVolume * Out_Region_13_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_13_log");
  G4LogicalVolume * Out_Region_14_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_14_log");
  G4LogicalVolume * Out_Region_15_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_15_log");
  G4LogicalVolume * Out_Region_16_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_16_log");
  G4LogicalVolume * Out_Region_17_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_17_log");
  G4LogicalVolume * Out_Region_18_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_18_log");
  G4LogicalVolume * Out_Region_19_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_19_log");
  G4LogicalVolume * Out_Region_20_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_20_log");
  G4LogicalVolume * Out_Region_21_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_21_log");
  G4LogicalVolume * Out_Region_22_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_22_log");
  G4LogicalVolume * Out_Region_23_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_23_log");
  G4LogicalVolume * Out_Region_24_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_24_log");
  G4LogicalVolume * Out_Region_25_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_25_log");
  G4LogicalVolume * Out_Region_26_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_26_log");
  G4LogicalVolume * Out_Region_27_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_27_log");
  G4LogicalVolume * Out_Region_28_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_28_log");
  G4LogicalVolume * Out_Region_29_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_29_log");
  G4LogicalVolume * Out_Region_30_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_30_log");
  G4LogicalVolume * Out_Region_31_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_31_log");
  G4LogicalVolume * Out_Region_32_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_32_log");
  G4LogicalVolume * Out_Region_33_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_33_log");
  G4LogicalVolume * Out_Region_34_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_34_log");
  G4LogicalVolume * Out_Region_35_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_35_log");
  G4LogicalVolume * Out_Region_36_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_36_log");
  G4LogicalVolume * Out_Region_37_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_37_log");
  G4LogicalVolume * Out_Region_38_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_38_log");
  G4LogicalVolume * Out_Region_39_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_39_log");
  G4LogicalVolume * Out_Region_40_log = new G4LogicalVolume(Out_Region,dummyMat,"Out_Region_40_log");


  fLogicalVolumeVector.push_back(Out_Region_1_log);
  fLogicalVolumeVector.push_back(Out_Region_2_log);
  fLogicalVolumeVector.push_back(Out_Region_3_log);
  fLogicalVolumeVector.push_back(Out_Region_4_log);
  fLogicalVolumeVector.push_back(Out_Region_5_log);
  fLogicalVolumeVector.push_back(Out_Region_6_log);
  fLogicalVolumeVector.push_back(Out_Region_7_log);
  fLogicalVolumeVector.push_back(Out_Region_8_log);
  fLogicalVolumeVector.push_back(Out_Region_9_log);
  fLogicalVolumeVector.push_back(Out_Region_10_log);
  fLogicalVolumeVector.push_back(Out_Region_11_log);
  fLogicalVolumeVector.push_back(Out_Region_12_log);
  fLogicalVolumeVector.push_back(Out_Region_13_log);
  fLogicalVolumeVector.push_back(Out_Region_14_log);
  fLogicalVolumeVector.push_back(Out_Region_15_log);
  fLogicalVolumeVector.push_back(Out_Region_16_log);
  fLogicalVolumeVector.push_back(Out_Region_17_log);
  fLogicalVolumeVector.push_back(Out_Region_18_log);
  fLogicalVolumeVector.push_back(Out_Region_19_log);
  fLogicalVolumeVector.push_back(Out_Region_20_log);
  fLogicalVolumeVector.push_back(Out_Region_21_log);
  fLogicalVolumeVector.push_back(Out_Region_22_log);
  fLogicalVolumeVector.push_back(Out_Region_23_log);
  fLogicalVolumeVector.push_back(Out_Region_24_log);
  fLogicalVolumeVector.push_back(Out_Region_25_log);
  fLogicalVolumeVector.push_back(Out_Region_26_log);
  fLogicalVolumeVector.push_back(Out_Region_27_log);
  fLogicalVolumeVector.push_back(Out_Region_28_log);
  fLogicalVolumeVector.push_back(Out_Region_29_log);
  fLogicalVolumeVector.push_back(Out_Region_30_log);
  fLogicalVolumeVector.push_back(Out_Region_31_log);
  fLogicalVolumeVector.push_back(Out_Region_32_log);
  fLogicalVolumeVector.push_back(Out_Region_33_log);
  fLogicalVolumeVector.push_back(Out_Region_34_log);
  fLogicalVolumeVector.push_back(Out_Region_35_log);
  fLogicalVolumeVector.push_back(Out_Region_36_log);
  fLogicalVolumeVector.push_back(Out_Region_37_log);
  fLogicalVolumeVector.push_back(Out_Region_38_log);
  fLogicalVolumeVector.push_back(Out_Region_39_log);
  fLogicalVolumeVector.push_back(Out_Region_40_log);
  fLogicalVolumeVector.push_back(notInteresting2_log);


  G4ThreeVector   NotInteresting_place(-1*worldX/2. + 0.5*(worldX/2. + startPos.x()), 0.,0.);
  G4ThreeVector   NotInteresting2_place( worldX/2. - 0.5*((worldX - (water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big) - (worldX/2. + startPos.x()))), 0.,0.);

  G4ThreeVector   Region_1_place(startPos.x() + 0.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_2_place(startPos.x() + 1.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_3_place(startPos.x() + 2.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_4_place(startPos.x() + 3.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_5_place(startPos.x() + 4.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_6_place(startPos.x() + 5.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_7_place(startPos.x() + 6.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_8_place(startPos.x() + 7.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_9_place(startPos.x() + 8.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_10_place(startPos.x() + 9.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_11_place(startPos.x() + 10.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_12_place(startPos.x() + 11.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_13_place(startPos.x() + 12.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_14_place(startPos.x() + 13.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_15_place(startPos.x() + 14.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_16_place(startPos.x() + 15.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_17_place(startPos.x() + 16.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_18_place(startPos.x() + 17.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_19_place(startPos.x() + 18.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_20_place(startPos.x() + 19.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_21_place(startPos.x() + 20.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_22_place(startPos.x() + 21.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_23_place(startPos.x() + 22.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_24_place(startPos.x() + 23.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_25_place(startPos.x() + 24.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_26_place(startPos.x() + 25.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_27_place(startPos.x() + 26.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_28_place(startPos.x() + 27.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_29_place(startPos.x() + 28.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_30_place(startPos.x() + 29.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_31_place(startPos.x() + 30.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_32_place(startPos.x() + 31.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_33_place(startPos.x() + 32.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_34_place(startPos.x() + 33.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_35_place(startPos.x() + 34.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_36_place(startPos.x() + 35.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_37_place(startPos.x() + 36.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_38_place(startPos.x() + 37.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_39_place(startPos.x() + 38.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());
  G4ThreeVector   Region_40_place(startPos.x() + 39.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), startPos.y(), startPos.z());

  G4ThreeVector   Out_Region_1_place(startPos.x() + 0.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_2_place(startPos.x() + 1.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_3_place(startPos.x() + 2.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_4_place(startPos.x() + 3.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_5_place(startPos.x() + 4.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_6_place(startPos.x() + 5.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_7_place(startPos.x() + 6.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_8_place(startPos.x() + 7.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_9_place(startPos.x() + 8.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_10_place(startPos.x() + 9.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_11_place(startPos.x() + 10.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_12_place(startPos.x() + 11.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_13_place(startPos.x() + 12.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_14_place(startPos.x() + 13.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_15_place(startPos.x() + 14.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_16_place(startPos.x() + 15.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_17_place(startPos.x() + 16.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_18_place(startPos.x() + 17.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_19_place(startPos.x() + 18.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_20_place(startPos.x() + 19.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_21_place(startPos.x() + 20.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_22_place(startPos.x() + 21.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_23_place(startPos.x() + 22.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_24_place(startPos.x() + 23.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_25_place(startPos.x() + 24.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_26_place(startPos.x() + 25.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_27_place(startPos.x() + 26.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_28_place(startPos.x() + 27.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_29_place(startPos.x() + 28.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_30_place(startPos.x() + 29.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_31_place(startPos.x() + 30.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_32_place(startPos.x() + 31.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_33_place(startPos.x() + 32.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_34_place(startPos.x() + 33.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_35_place(startPos.x() + 34.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_36_place(startPos.x() + 35.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_37_place(startPos.x() + 36.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_38_place(startPos.x() + 37.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_39_place(startPos.x() + 38.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions),0.,0.);
  G4ThreeVector   Out_Region_40_place(startPos.x() + 39.5*((water_thick + gb_depth + SSPlate_thick + thermalC_length_small + thermalC_length_big)/nRegions), 0.,0.);


  G4VPhysicalVolume *pvol1 = new G4PVPlacement(0,NotInteresting_place,notInteresting_log,"notInteresting_phys",worldLogical,true,0);
  G4VPhysicalVolume *pvol2 = new G4PVPlacement(0,Region_1_place , Region_1_log,"Region_1_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol3 = new G4PVPlacement(0,Region_2_place , Region_2_log,"Region_2_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol4 = new G4PVPlacement(0,Region_3_place , Region_3_log,"Region_3_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol5 = new G4PVPlacement(0,Region_4_place , Region_4_log,"Region_4_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol6 = new G4PVPlacement(0,Region_5_place , Region_5_log,"Region_5_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol7 = new G4PVPlacement(0,Region_6_place , Region_6_log,"Region_6_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol8 = new G4PVPlacement(0,Region_7_place , Region_7_log,"Region_7_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol9 = new G4PVPlacement(0,Region_8_place , Region_8_log,"Region_8_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol10 = new G4PVPlacement(0,Region_9_place , Region_9_log,"Region_9_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol11 = new G4PVPlacement(0,Region_10_place , Region_10_log,"Region_10_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol12 = new G4PVPlacement(0,Region_11_place , Region_11_log,"Region_11_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol13 = new G4PVPlacement(0,Region_12_place , Region_12_log,"Region_12_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol14 = new G4PVPlacement(0,Region_13_place , Region_13_log,"Region_13_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol15 = new G4PVPlacement(0,Region_14_place , Region_14_log,"Region_14_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol16 = new G4PVPlacement(0,Region_15_place , Region_15_log,"Region_15_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol17 = new G4PVPlacement(0,Region_16_place , Region_16_log,"Region_16_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol18 = new G4PVPlacement(0,Region_17_place , Region_17_log,"Region_17_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol19 = new G4PVPlacement(0,Region_18_place , Region_18_log,"Region_18_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol20 = new G4PVPlacement(0,Region_19_place , Region_19_log,"Region_19_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol21 = new G4PVPlacement(0,Region_20_place , Region_20_log,"Region_20_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol22 = new G4PVPlacement(0,Region_21_place , Region_21_log,"Region_21_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol23 = new G4PVPlacement(0,Region_22_place , Region_22_log,"Region_22_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol24 = new G4PVPlacement(0,Region_23_place , Region_23_log,"Region_23_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol25 = new G4PVPlacement(0,Region_24_place , Region_24_log,"Region_24_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol26 = new G4PVPlacement(0,Region_25_place , Region_25_log,"Region_25_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol27 = new G4PVPlacement(0,Region_26_place , Region_26_log,"Region_26_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol28 = new G4PVPlacement(0,Region_27_place , Region_27_log,"Region_27_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol29 = new G4PVPlacement(0,Region_28_place , Region_28_log,"Region_28_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol30 = new G4PVPlacement(0,Region_29_place , Region_29_log,"Region_29_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol31 = new G4PVPlacement(0,Region_30_place , Region_30_log,"Region_30_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol32 = new G4PVPlacement(0,Region_31_place , Region_31_log,"Region_31_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol33 = new G4PVPlacement(0,Region_32_place , Region_32_log,"Region_32_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol34 = new G4PVPlacement(0,Region_33_place , Region_33_log,"Region_33_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol35 = new G4PVPlacement(0,Region_34_place , Region_34_log,"Region_34_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol36 = new G4PVPlacement(0,Region_35_place , Region_35_log,"Region_35_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol37 = new G4PVPlacement(0,Region_36_place , Region_36_log,"Region_36_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol38 = new G4PVPlacement(0,Region_37_place , Region_37_log,"Region_37_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol39 = new G4PVPlacement(0,Region_38_place , Region_38_log,"Region_38_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol40 = new G4PVPlacement(0,Region_39_place , Region_39_log,"Region_39_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol41 = new G4PVPlacement(0,Region_40_place , Region_40_log,"Region_40_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol42 = new G4PVPlacement(0,Out_Region_1_place , Out_Region_1_log,"Out_Region_1_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol43 = new G4PVPlacement(0,Out_Region_2_place , Out_Region_2_log,"Out_Region_2_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol44 = new G4PVPlacement(0,Out_Region_3_place , Out_Region_3_log,"Out_Region_3_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol45 = new G4PVPlacement(0,Out_Region_4_place , Out_Region_4_log,"Out_Region_4_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol46 = new G4PVPlacement(0,Out_Region_5_place , Out_Region_5_log,"Out_Region_5_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol47 = new G4PVPlacement(0,Out_Region_6_place , Out_Region_6_log,"Out_Region_6_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol48 = new G4PVPlacement(0,Out_Region_7_place , Out_Region_7_log,"Out_Region_7_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol49 = new G4PVPlacement(0,Out_Region_8_place , Out_Region_8_log,"Out_Region_8_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol50 = new G4PVPlacement(0,Out_Region_9_place , Out_Region_9_log,"Out_Region_9_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol51 = new G4PVPlacement(0,Out_Region_10_place , Out_Region_10_log,"Out_Region_10_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol52 = new G4PVPlacement(0,Out_Region_11_place , Out_Region_11_log,"Out_Region_11_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol53 = new G4PVPlacement(0,Out_Region_12_place , Out_Region_12_log,"Out_Region_12_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol54 = new G4PVPlacement(0,Out_Region_13_place , Out_Region_13_log,"Out_Region_13_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol55 = new G4PVPlacement(0,Out_Region_14_place , Out_Region_14_log,"Out_Region_14_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol56 = new G4PVPlacement(0,Out_Region_15_place , Out_Region_15_log,"Out_Region_15_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol57 = new G4PVPlacement(0,Out_Region_16_place , Out_Region_16_log,"Out_Region_16_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol58 = new G4PVPlacement(0,Out_Region_17_place , Out_Region_17_log,"Out_Region_17_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol59 = new G4PVPlacement(0,Out_Region_18_place , Out_Region_18_log,"Out_Region_18_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol60 = new G4PVPlacement(0,Out_Region_19_place , Out_Region_19_log,"Out_Region_19_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol61 = new G4PVPlacement(0,Out_Region_20_place , Out_Region_20_log,"Out_Region_20_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol62 = new G4PVPlacement(0,Out_Region_21_place , Out_Region_21_log,"Out_Region_21_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol63 = new G4PVPlacement(0,Out_Region_22_place , Out_Region_22_log,"Out_Region_22_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol64 = new G4PVPlacement(0,Out_Region_23_place , Out_Region_23_log,"Out_Region_23_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol65 = new G4PVPlacement(0,Out_Region_24_place , Out_Region_24_log,"Out_Region_24_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol66 = new G4PVPlacement(0,Out_Region_25_place , Out_Region_25_log,"Out_Region_25_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol67 = new G4PVPlacement(0,Out_Region_26_place , Out_Region_26_log,"Out_Region_26_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol68 = new G4PVPlacement(0,Out_Region_27_place , Out_Region_27_log,"Out_Region_27_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol69 = new G4PVPlacement(0,Out_Region_28_place , Out_Region_28_log,"Out_Region_28_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol70 = new G4PVPlacement(0,Out_Region_29_place , Out_Region_29_log,"Out_Region_29_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol71 = new G4PVPlacement(0,Out_Region_30_place , Out_Region_30_log,"Out_Region_30_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol72 = new G4PVPlacement(0,Out_Region_31_place , Out_Region_31_log,"Out_Region_31_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol73 = new G4PVPlacement(0,Out_Region_32_place , Out_Region_32_log,"Out_Region_32_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol74 = new G4PVPlacement(0,Out_Region_33_place , Out_Region_33_log,"Out_Region_33_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol75 = new G4PVPlacement(0,Out_Region_34_place , Out_Region_34_log,"Out_Region_34_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol76 = new G4PVPlacement(0,Out_Region_35_place , Out_Region_35_log,"Out_Region_35_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol77 = new G4PVPlacement(0,Out_Region_36_place , Out_Region_36_log,"Out_Region_36_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol78 = new G4PVPlacement(0,Out_Region_37_place , Out_Region_37_log,"Out_Region_37_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol79 = new G4PVPlacement(0,Out_Region_38_place , Out_Region_38_log,"Out_Region_38_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol80 = new G4PVPlacement(0,Out_Region_39_place , Out_Region_39_log,"Out_Region_39_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol81 = new G4PVPlacement(0,Out_Region_40_place , Out_Region_40_log,"Out_Region_40_phys" ,worldLogical,true,0);
  G4VPhysicalVolume *pvol82= new G4PVPlacement(0,NotInteresting2_place,notInteresting2_log,"notInteresting2_phys",worldLogical,true,0);


  G4GeometryCell cell1(*pvol1, 0);
  G4GeometryCell cell2(*pvol2, 0);
  G4GeometryCell cell3(*pvol3, 0);
  G4GeometryCell cell4(*pvol4, 0);
  G4GeometryCell cell5(*pvol5, 0);
  G4GeometryCell cell6(*pvol6, 0);
  G4GeometryCell cell7(*pvol7, 0);
  G4GeometryCell cell8(*pvol8, 0);
  G4GeometryCell cell9(*pvol9, 0);
  G4GeometryCell cell10(*pvol10, 0);
  G4GeometryCell cell11(*pvol11, 0);
  G4GeometryCell cell12(*pvol12, 0);
  G4GeometryCell cell13(*pvol13, 0);
  G4GeometryCell cell14(*pvol14, 0);
  G4GeometryCell cell15(*pvol15, 0);
  G4GeometryCell cell16(*pvol16, 0);
  G4GeometryCell cell17(*pvol17, 0);
  G4GeometryCell cell18(*pvol18, 0);
  G4GeometryCell cell19(*pvol19, 0);
  G4GeometryCell cell20(*pvol20, 0);
  G4GeometryCell cell21(*pvol21, 0);
  G4GeometryCell cell22(*pvol22, 0);
  G4GeometryCell cell23(*pvol23, 0);
  G4GeometryCell cell24(*pvol24, 0);
  G4GeometryCell cell25(*pvol25, 0);
  G4GeometryCell cell26(*pvol26, 0);
  G4GeometryCell cell27(*pvol27, 0);
  G4GeometryCell cell28(*pvol28, 0);
  G4GeometryCell cell29(*pvol29, 0);
  G4GeometryCell cell30(*pvol30, 0);
  G4GeometryCell cell31(*pvol31, 0);
  G4GeometryCell cell32(*pvol32, 0);
  G4GeometryCell cell33(*pvol33, 0);
  G4GeometryCell cell34(*pvol34, 0);
  G4GeometryCell cell35(*pvol35, 0);
  G4GeometryCell cell36(*pvol36, 0);
  G4GeometryCell cell37(*pvol37, 0);
  G4GeometryCell cell38(*pvol38, 0);
  G4GeometryCell cell39(*pvol39, 0);
  G4GeometryCell cell40(*pvol40, 0);
  G4GeometryCell cell41(*pvol41, 0);
  G4GeometryCell cell42(*pvol42, 0); 
  G4GeometryCell cell43(*pvol43, 0); 
  G4GeometryCell cell44(*pvol44, 0); 
  G4GeometryCell cell45(*pvol45, 0); 
  G4GeometryCell cell46(*pvol46, 0); 
  G4GeometryCell cell47(*pvol47, 0); 
  G4GeometryCell cell48(*pvol48, 0); 
  G4GeometryCell cell49(*pvol49, 0); 
  G4GeometryCell cell50(*pvol50, 0); 
  G4GeometryCell cell51(*pvol51, 0); 
  G4GeometryCell cell52(*pvol52, 0); 
  G4GeometryCell cell53(*pvol53, 0); 
  G4GeometryCell cell54(*pvol54, 0); 
  G4GeometryCell cell55(*pvol55, 0); 
  G4GeometryCell cell56(*pvol56, 0); 
  G4GeometryCell cell57(*pvol57, 0); 
  G4GeometryCell cell58(*pvol58, 0); 
  G4GeometryCell cell59(*pvol59, 0); 
  G4GeometryCell cell60(*pvol60, 0); 
  G4GeometryCell cell61(*pvol61, 0); 
  G4GeometryCell cell62(*pvol62, 0); 
  G4GeometryCell cell63(*pvol63, 0); 
  G4GeometryCell cell64(*pvol64, 0); 
  G4GeometryCell cell65(*pvol65, 0); 
  G4GeometryCell cell66(*pvol66, 0); 
  G4GeometryCell cell67(*pvol67, 0); 
  G4GeometryCell cell68(*pvol68, 0); 
  G4GeometryCell cell69(*pvol69, 0); 
  G4GeometryCell cell70(*pvol70, 0); 
  G4GeometryCell cell71(*pvol71, 0); 
  G4GeometryCell cell72(*pvol72, 0); 
  G4GeometryCell cell73(*pvol73, 0); 
  G4GeometryCell cell74(*pvol74, 0); 
  G4GeometryCell cell75(*pvol75, 0); 
  G4GeometryCell cell76(*pvol76, 0); 
  G4GeometryCell cell77(*pvol77, 0); 
  G4GeometryCell cell78(*pvol78, 0); 
  G4GeometryCell cell79(*pvol79, 0); 
  G4GeometryCell cell80(*pvol80, 0); 
  G4GeometryCell cell81(*pvol81, 0); 
  G4GeometryCell cell82(*pvol82, 0);




  fPVolumeStore.AddPVolume(cell1);
  fPVolumeStore.AddPVolume(cell2);
  fPVolumeStore.AddPVolume(cell3);
  fPVolumeStore.AddPVolume(cell4);
  fPVolumeStore.AddPVolume(cell5);
  fPVolumeStore.AddPVolume(cell6);
  fPVolumeStore.AddPVolume(cell7);
  fPVolumeStore.AddPVolume(cell8);
  fPVolumeStore.AddPVolume(cell9);
  fPVolumeStore.AddPVolume(cell10);
  fPVolumeStore.AddPVolume(cell11);
  fPVolumeStore.AddPVolume(cell12);
  fPVolumeStore.AddPVolume(cell13);
  fPVolumeStore.AddPVolume(cell14);
  fPVolumeStore.AddPVolume(cell15);
  fPVolumeStore.AddPVolume(cell16);
  fPVolumeStore.AddPVolume(cell17);
  fPVolumeStore.AddPVolume(cell18);
  fPVolumeStore.AddPVolume(cell19);
  fPVolumeStore.AddPVolume(cell20);
  fPVolumeStore.AddPVolume(cell21);
  fPVolumeStore.AddPVolume(cell22);
  fPVolumeStore.AddPVolume(cell23);
  fPVolumeStore.AddPVolume(cell24);
  fPVolumeStore.AddPVolume(cell25);
  fPVolumeStore.AddPVolume(cell26);
  fPVolumeStore.AddPVolume(cell27);
  fPVolumeStore.AddPVolume(cell28);
  fPVolumeStore.AddPVolume(cell29);
  fPVolumeStore.AddPVolume(cell30);
  fPVolumeStore.AddPVolume(cell31);
  fPVolumeStore.AddPVolume(cell32);
  fPVolumeStore.AddPVolume(cell33);
  fPVolumeStore.AddPVolume(cell34);
  fPVolumeStore.AddPVolume(cell35);
  fPVolumeStore.AddPVolume(cell36);
  fPVolumeStore.AddPVolume(cell37);
  fPVolumeStore.AddPVolume(cell38);
  fPVolumeStore.AddPVolume(cell39);
  fPVolumeStore.AddPVolume(cell40);
  fPVolumeStore.AddPVolume(cell41);
  fPVolumeStore.AddPVolume(cell42);
  fPVolumeStore.AddPVolume(cell43);
  fPVolumeStore.AddPVolume(cell44);
  fPVolumeStore.AddPVolume(cell45);
  fPVolumeStore.AddPVolume(cell46);
  fPVolumeStore.AddPVolume(cell47);
  fPVolumeStore.AddPVolume(cell48);
  fPVolumeStore.AddPVolume(cell49);
  fPVolumeStore.AddPVolume(cell50);
  fPVolumeStore.AddPVolume(cell51);
  fPVolumeStore.AddPVolume(cell52);
  fPVolumeStore.AddPVolume(cell53);
  fPVolumeStore.AddPVolume(cell54);
  fPVolumeStore.AddPVolume(cell55);
  fPVolumeStore.AddPVolume(cell56);
  fPVolumeStore.AddPVolume(cell57);
  fPVolumeStore.AddPVolume(cell58);
  fPVolumeStore.AddPVolume(cell59);
  fPVolumeStore.AddPVolume(cell60);
  fPVolumeStore.AddPVolume(cell61);
  fPVolumeStore.AddPVolume(cell62);
  fPVolumeStore.AddPVolume(cell63);
  fPVolumeStore.AddPVolume(cell64);
  fPVolumeStore.AddPVolume(cell65);
  fPVolumeStore.AddPVolume(cell66);
  fPVolumeStore.AddPVolume(cell67);
  fPVolumeStore.AddPVolume(cell68);
  fPVolumeStore.AddPVolume(cell69);
  fPVolumeStore.AddPVolume(cell70);
  fPVolumeStore.AddPVolume(cell71);
  fPVolumeStore.AddPVolume(cell72);
  fPVolumeStore.AddPVolume(cell73);
  fPVolumeStore.AddPVolume(cell74);
  fPVolumeStore.AddPVolume(cell75);
  fPVolumeStore.AddPVolume(cell76);
  fPVolumeStore.AddPVolume(cell77);
  fPVolumeStore.AddPVolume(cell78);
  fPVolumeStore.AddPVolume(cell79);
  fPVolumeStore.AddPVolume(cell80);
  fPVolumeStore.AddPVolume(cell81);
  fPVolumeStore.AddPVolume(cell82);




}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4VPhysicalVolume &ImportanceGeometryConstruction::
GetPhysicalVolumeByName(const G4String& name) const {
  return *fPVolumeStore.GetPVolume(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String ImportanceGeometryConstruction::ListPhysNamesAsG4String(){
  G4String names(fPVolumeStore.GetPNames());
  return names;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String ImportanceGeometryConstruction::GetCellName(G4int i) {
  if (i == 1) { return "notInteresting_phys"; }
  if (i == 2) { return "Region_1_phys"; }
  if (i == 3) { return "Region_2_phys"; }
  if (i == 4) { return "Region_3_phys"; }
  if (i == 5) { return "Region_4_phys"; }
  if (i == 6) { return "Region_5_phys"; }
  if (i == 7) { return "Region_6_phys"; }
  if (i == 8) { return "Region_7_phys"; }
  if (i == 9) { return "Region_8_phys"; }
  if (i == 10) { return "Region_9_phys"; }
  if (i == 11) { return "Region_10_phys"; }
  if (i == 12) { return "Region_11_phys"; }
  if (i == 13) { return "Region_12_phys"; }
  if (i == 14) { return "Region_13_phys"; }
  if (i == 15) { return "Region_14_phys"; }
  if (i == 16) { return "Region_15_phys"; }
  if (i == 17) { return "Region_16_phys"; }
  if (i == 18) { return "Region_17_phys"; }
  if (i == 19) { return "Region_18_phys"; }
  if (i == 20) { return "Region_19_phys"; }
  if (i == 21) { return "Region_20_phys"; }
  if (i == 22) { return "Region_21_phys"; }
  if (i == 23) { return "Region_22_phys"; }
  if (i == 24) { return "Region_23_phys"; }
  if (i == 25) { return "Region_24_phys"; }
  if (i == 26) { return "Region_25_phys"; }
  if (i == 27) { return "Region_26_phys"; }
  if (i == 28) { return "Region_27_phys"; }
  if (i == 29) { return "Region_28_phys"; }
  if (i == 30) { return "Region_29_phys"; }
  if (i == 31) { return "Region_30_phys"; }
  if (i == 32) { return "Region_31_phys"; }
  if (i == 33) { return "Region_32_phys"; }
  if (i == 34) { return "Region_33_phys"; }
  if (i == 35) { return "Region_34_phys"; }
  if (i == 36) { return "Region_35_phys"; }
  if (i == 37) { return "Region_36_phys"; }
  if (i == 38) { return "Region_37_phys"; }
  if (i == 39) { return "Region_38_phys"; }
  if (i == 40) { return "Region_39_phys"; }
  if (i == 41) { return "Region_40_phys"; }
  if (i == 42) { return "Out_Region_1_phys"; }
  if (i == 43) { return "Out_Region_2_phys"; }
  if (i == 44) { return "Out_Region_3_phys"; }
  if (i == 45) { return "Out_Region_4_phys"; }
  if (i == 46) { return "Out_Region_5_phys"; }
  if (i == 47) { return "Out_Region_6_phys"; }
  if (i == 48) { return "Out_Region_7_phys"; }
  if (i == 49) { return "Out_Region_8_phys"; }
  if (i == 50) { return "Out_Region_9_phys"; }
  if (i == 51) { return "Out_Region_10_phys"; }
  if (i == 52) { return "Out_Region_11_phys"; }
  if (i == 53) { return "Out_Region_12_phys"; }
  if (i == 54) { return "Out_Region_13_phys"; }
  if (i == 55) { return "Out_Region_14_phys"; }
  if (i == 56) { return "Out_Region_15_phys"; }
  if (i == 57) { return "Out_Region_16_phys"; }
  if (i == 58) { return "Out_Region_17_phys"; }
  if (i == 59) { return "Out_Region_18_phys"; }
  if (i == 60) { return "Out_Region_19_phys"; }
  if (i == 61) { return "Out_Region_20_phys"; }
  if (i == 62) { return "Out_Region_21_phys"; }
  if (i == 63) { return "Out_Region_22_phys"; }
  if (i == 64) { return "Out_Region_23_phys"; }
  if (i == 65) { return "Out_Region_24_phys"; }
  if (i == 66) { return "Out_Region_25_phys"; }
  if (i == 67) { return "Out_Region_26_phys"; }
  if (i == 68) { return "Out_Region_27_phys"; }
  if (i == 69) { return "Out_Region_28_phys"; }
  if (i == 70) { return "Out_Region_29_phys"; }
  if (i == 71) { return "Out_Region_30_phys"; }
  if (i == 72) { return "Out_Region_31_phys"; }
  if (i == 73) { return "Out_Region_32_phys"; }
  if (i == 74) { return "Out_Region_33_phys"; }
  if (i == 75) { return "Out_Region_34_phys"; }
  if (i == 76) { return "Out_Region_35_phys"; }
  if (i == 77) { return "Out_Region_36_phys"; }
  if (i == 78) { return "Out_Region_37_phys"; }
  if (i == 79) { return "Out_Region_38_phys"; }
  if (i == 80) { return "Out_Region_39_phys"; }
  if (i == 81) { return "Out_Region_40_phys"; }
  if (i == 82) { return "notInteresting2_phys"; }
  return "";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4GeometryCell ImportanceGeometryConstruction::GetGeometryCell(G4int i){
  G4String name(GetCellName(i));
  const G4VPhysicalVolume *p=0;
  p = fPVolumeStore.GetPVolume(name);
  if (p) {
    return G4GeometryCell(*p,0);
  }
  else {
    G4cout << "ImportanceGeometryConstruction::GetGeometryCell: " << G4endl
           << " couldn't get G4GeometryCell" << G4endl;
    return G4GeometryCell(*fGhostWorld,-2);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume &ImportanceGeometryConstruction::
GetWorldVolumeAddress() const{
   return *fGhostWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *ImportanceGeometryConstruction::GetWorldVolume() {
  return fGhostWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ImportanceGeometryConstruction::SetSensitive(){

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ImportanceGeometryConstruction::ConstructSD()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VIStore* ImportanceGeometryConstruction::CreateImportanceStore()
{


  G4cout << " ImportanceGeometryConstruction:: Creating Importance Store " 
         << G4endl;
  if (!fPVolumeStore.Size())
  {
    G4Exception("ImportanceGeometryConstruction::CreateImportanceStore"
               ,"example_0001",RunMustBeAborted
               ,"no physical volumes created yet!");
  }

  // creating and filling the importance store
  
  //  G4IStore *istore = new G4IStore(*fWorldVolume);
  
  G4IStore *istore = G4IStore::GetInstance(GetName());
  

  G4GeometryCell gWorldVolumeCell(GetWorldVolumeAddress(), 0);


  istore->AddImportanceGeometryCell(1, GetGeometryCell(1).GetPhysicalVolume());

  // the interesting region (thermal column, graphite, water between core and TC, etc.
  istore->AddImportanceGeometryCell(std::pow(2,1), GetGeometryCell(2).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,1), GetGeometryCell(3).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,1), GetGeometryCell(4).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,2), GetGeometryCell(5).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,2), GetGeometryCell(6).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,2), GetGeometryCell(7).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,3), GetGeometryCell(8).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,3), GetGeometryCell(9).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,3), GetGeometryCell(10).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,4), GetGeometryCell(11).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,4), GetGeometryCell(12).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,4), GetGeometryCell(13).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,5), GetGeometryCell(14).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,5), GetGeometryCell(15).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,5), GetGeometryCell(16).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,6), GetGeometryCell(17).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,6), GetGeometryCell(18).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,6), GetGeometryCell(19).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,8), GetGeometryCell(20).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,8), GetGeometryCell(21).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,9), GetGeometryCell(22).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,10), GetGeometryCell(23).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,11), GetGeometryCell(24).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,12), GetGeometryCell(25).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,13), GetGeometryCell(26).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,14), GetGeometryCell(27).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,15), GetGeometryCell(28).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,16), GetGeometryCell(29).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,17), GetGeometryCell(30).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,18), GetGeometryCell(31).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,19), GetGeometryCell(32).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,20), GetGeometryCell(33).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,21), GetGeometryCell(34).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,22), GetGeometryCell(35).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,22), GetGeometryCell(36).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,22), GetGeometryCell(37).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,22), GetGeometryCell(38).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,22), GetGeometryCell(39).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,22), GetGeometryCell(40).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,22), GetGeometryCell(41).GetPhysicalVolume());

  // these are regions outside the interesting region but touching them
  //  should be close but less than the inner part they touch in importance value
  // this helps reduce unwanted tracks but avoids transitions with huge importance value differences
  istore->AddImportanceGeometryCell(std::pow(2,1), GetGeometryCell(42).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,1), GetGeometryCell(43).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,1), GetGeometryCell(44).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,1), GetGeometryCell(45).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,1), GetGeometryCell(46).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,1), GetGeometryCell(47).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,1), GetGeometryCell(48).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,1), GetGeometryCell(49).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,1), GetGeometryCell(50).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,2), GetGeometryCell(51).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,2), GetGeometryCell(52).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,2), GetGeometryCell(53).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,3), GetGeometryCell(54).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,3), GetGeometryCell(55).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,3), GetGeometryCell(56).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,4), GetGeometryCell(57).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,4), GetGeometryCell(58).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,4), GetGeometryCell(59).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,5), GetGeometryCell(60).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,6), GetGeometryCell(61).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,7), GetGeometryCell(62).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,8), GetGeometryCell(63).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,9), GetGeometryCell(64).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,10), GetGeometryCell(65).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,11), GetGeometryCell(66).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,12), GetGeometryCell(67).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,13), GetGeometryCell(68).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,14), GetGeometryCell(69).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,15), GetGeometryCell(70).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,16), GetGeometryCell(71).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,17), GetGeometryCell(72).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,18), GetGeometryCell(73).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,19), GetGeometryCell(74).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,20), GetGeometryCell(75).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,20), GetGeometryCell(76).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,20), GetGeometryCell(77).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,20), GetGeometryCell(78).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,20), GetGeometryCell(79).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,20), GetGeometryCell(80).GetPhysicalVolume());
  istore->AddImportanceGeometryCell(std::pow(2,20), GetGeometryCell(81).GetPhysicalVolume());

  // finally, the region behind the TC
  istore->AddImportanceGeometryCell(std::pow(2,8), GetGeometryCell(82).GetPhysicalVolume());

  return istore;

}

