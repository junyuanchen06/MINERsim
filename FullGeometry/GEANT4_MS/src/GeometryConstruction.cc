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
// $Id: GeometryConstruction.cc 67272 2013-02-13 12:05:47Z ihrivnac $
//
/// \file eventgenerator/exgps/src/GeometryConstruction.cc
/// \brief Implementation of the GeometryConstruction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "GeometryConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "MinerSD.hh"
#include "MINERMaterials.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4SDManager.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GeometryConstruction::GeometryConstruction()
: G4VUserDetectorConstruction(),
  fUniverse_phys(0)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GeometryConstruction::~GeometryConstruction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GeometryConstruction::ConstructSDandField()
{ 

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  // Sensitive detectors
  MinerSD* aCraigSD = new MinerSD("/MINERsim/Craig", "MS_Craig_hits");
  SDman->AddNewDetector(aCraigSD);
  SetSensitiveDetector("Craig_log_insideCastle", aCraigSD, true);

  MinerSD* aCraigSD2 = new MinerSD("/MINERsim/Craig2", "MS_Craig_hits2");
  SDman->AddNewDetector(aCraigSD2);
  SetSensitiveDetector("Craig_log_outsideCastle", aCraigSD2, true);


}

G4AssemblyVolume* GeometryConstruction::ConstructBEGe(std::string name)
{

  // Craig's Ge detector
  // Need to put N2 into the detector still

  G4double detRad = 50.8*mm/2.; 
  G4double detHalfZ = (4.805/2.)*cm; 

  G4VSolid* target1 = new G4Tubs("Target1"+name,0.,detRad,detHalfZ,0,360*deg);
  G4VSolid* target2 = new G4Tubs("Target2"+name,0.,4.75*mm,15.0*mm,0,360*deg);
  G4ThreeVector targetOffset(0,0,-1*(detHalfZ - 15.0*mm));
  G4SubtractionSolid *target = new G4SubtractionSolid("Target"+name,target1,target2,0,targetOffset);
  G4LogicalVolume *fLogicDet = new G4LogicalVolume(target, mats->GetDetGe(),"Craig_log"+name);


  G4double cuCasingThick = 0.76*mm;
  //G4double cuCasingThick = 2.25*mm;
  G4double cuRingThick = 2.7*mm - 0.76*mm;
  //G4double cuRingThick = 3.5*mm;
  G4double cuRingHeight = 8.6*mm;


  G4VSolid* cuCasing1 = new G4Tubs("cuCasing1"+name,0.,detRad+cuCasingThick,35.5*mm/2. + detHalfZ + cuCasingThick/2.,0,360*deg);
  G4VSolid* cuCasing2 = new G4Tubs("cuCasing2"+name,0.,detRad+.001*mm,35.5*mm/2. + detHalfZ + cuCasingThick/2.,0,360*deg);
  G4ThreeVector casingOffset(0,0,3.0*mm);
  G4SubtractionSolid *cuCasing = new G4SubtractionSolid("cuCasting"+name,cuCasing1,cuCasing2,0,casingOffset);


  G4VSolid* cuRing = new G4Tubs("cuRing"+name,detRad+cuCasingThick,detRad+cuCasingThick+cuRingThick,cuRingHeight/2.,0,360*deg);
  G4VSolid* cuRingT = new G4Tubs("cuRingT"+name,detRad+cuCasingThick,detRad+cuCasingThick+cuRingThick,cuRingHeight/6.,0,360*deg);
  G4ThreeVector ringOffsetT(0,0,35.5*mm/2. + detHalfZ + cuCasingThick/2.-cuRingHeight/6.);
  G4ThreeVector ringOffset1(0,0,(35.5*mm/2. + detHalfZ + cuCasingThick/2.)/3.);
  G4ThreeVector ringOffset2(0,0,-1*(35.5*mm/2. + detHalfZ + cuCasingThick/2.)/3.);

  G4UnionSolid *cuDetHousing1 = new G4UnionSolid("cuDetHousing1"+name,cuCasing,cuRingT,0,ringOffsetT);
  G4UnionSolid *cuDetHousing2 = new G4UnionSolid("cuDetHousing2"+name,cuDetHousing1,cuRing,0,ringOffset1);
  G4UnionSolid *cuDetHousing = new G4UnionSolid("cuDetHousing"+name,cuDetHousing2,cuRing,0,ringOffset2);

  G4LogicalVolume *fLogicCasing = new G4LogicalVolume(cuDetHousing, mats->GetCopper(),"DetHousing_log"+name);
  G4ThreeVector cuHousingPos(0,0,-1*((35.5*mm/2. + detHalfZ + cuCasingThick/2.) - detHalfZ));




  G4double tubeAthick = (1./16.)*in;
  G4double tubeAlen = 4.38*in;
  G4double tubeAdia = 3.0*in;

  G4double tubeBthick = (1./8.)*in;
  G4double tubeBlen = 5.0*in;
  G4double tubeBdia = 3.13*in;

  G4double tubeCthick = (1./8.)*in;
  G4double tubeClen = 0.88*in;
  G4double tubeCdia = 1.13*in;

  G4double tubeDthick = (1./8.)*in;
  G4double tubeDlen = 13.75*in;
  G4double tubeDdia = 9.00*in;

  G4Tubs *outTubeA = new G4Tubs("outTubeA"+name, 0., tubeAdia/2., tubeAlen/2.,0,360*deg);
  G4Tubs *inTubeA = new G4Tubs("inTubeA"+name, 0., tubeAdia/2. - tubeAthick, tubeAlen/2. - tubeAthick/2.,0,360*deg);
  // right now the above defines the aluminum tube surrounding the detector, and has a 1/16" aluminum window, probably too thick?
  G4ThreeVector inTubeAOffset(0,0,-tubeAthick/2.);
  G4SubtractionSolid *detTubeA = new G4SubtractionSolid("detTubeA"+name,outTubeA,inTubeA,0,inTubeAOffset);

  G4Tubs *outTubeB = new G4Tubs("outTubeB"+name, 0., tubeBdia/2., tubeBlen/2.,0,360*deg);
  G4Tubs *inTubeB = new G4Tubs("inTubeB"+name, 0., tubeBdia/2. - tubeBthick, tubeBlen/2. - tubeBthick/2.,0,360*deg);
  G4ThreeVector inTubeBOffset(0,0,+tubeBthick/2.);
  G4SubtractionSolid *detTubeB = new G4SubtractionSolid("detTubeB"+name,outTubeB,inTubeB,0,inTubeBOffset);

  G4Tubs *detTubeC = new G4Tubs("detTubeC"+name, tubeCdia/2. - tubeCthick, tubeCdia/2., tubeClen/2.,0,360*deg);
  G4Tubs *inTubeC = new G4Tubs("inTubeC"+name, 0., tubeCdia/2. - tubeCthick, tubeClen/2.,0,360*deg);

  G4Tubs *outTubeD = new G4Tubs("outTubeD"+name, 0., tubeDdia/2., tubeDlen/2.,0,360*deg);
  G4Tubs *inTubeD = new G4Tubs("inTubeD"+name, 0., tubeDdia/2. - tubeDthick, tubeDlen/2. - tubeDthick,0,360*deg);
  G4SubtractionSolid *detTubeD = new G4SubtractionSolid("detTubeD"+name,outTubeD,inTubeD);


  G4AssemblyVolume *fullDet = new G4AssemblyVolume();
  G4ThreeVector tubeAPos(0,0,                 -1*(tubeAlen/2. - tubeAthick - detHalfZ) + 5*mm); //detector 5mm from front aluminum
  G4ThreeVector tubeBPos(0,0,                 (tubeAPos.z() - tubeAlen/2. - tubeBlen/2.));
  G4ThreeVector tubeCPos(0,0 - (tubeBdia/3.), (tubeAPos.z() - tubeAlen/2. - tubeBlen - tubeClen/2.));
  G4ThreeVector tubeDPos(0,0,                 (tubeAPos.z() - tubeAlen/2. - tubeBlen - tubeClen - tubeDlen/2.));

  G4ThreeVector tubeA_N2_Pos = tubeAPos+inTubeAOffset;
  G4ThreeVector tubeB_N2_Pos = tubeBPos+inTubeBOffset;

  G4Tubs *liqN2_TubeD = new G4Tubs("liqN2_TubeD"+name, 0., tubeDdia/2. - tubeDthick, tubeDlen/2. - tubeDthick,0,360*deg);

  // canberra lead shield
  G4Tubs *InnerCanberra = new G4Tubs("InnerCanberra"+name, 0., (279.*mm)/2.,203.*mm,0.,360*deg);
  G4Tubs *CopperCanberraI = new G4Tubs("CopperCanberraI"+name, 0., (279.*mm)/2. + 1.6*mm,203.*mm + 1.6*mm,0.,360*deg);
  G4Tubs *TinCanberraI = new G4Tubs("TinCanberraI"+name, 0., (279.*mm)/2. + 1.6*mm + 1.0*mm,203.*mm + 1.6*mm + 1.0*mm,0.,360*deg);
  G4Tubs *LeadCanberraI = new G4Tubs("LeadCanberraI"+name, 0., (279.*mm)/2. + 1.6*mm + 1.0*mm + 4.*in,203.*mm + 1.6*mm + 1.0*mm + 4.*in,0.,360*deg);
  G4Tubs *SteelCanberraI = new G4Tubs("SteelCanberraI"+name, 0., (279.*mm)/2. + 1.6*mm + 1.0*mm + 4.*in + 9.5*mm,203.*mm + 1.6*mm + 1.0*mm + 4.*in + 9.5*mm,0.,360*deg);

  G4SubtractionSolid *CopperCanberraI2 = new G4SubtractionSolid("CopperCanberraI2"+name,CopperCanberraI,InnerCanberra);
  G4SubtractionSolid *TinCanberraI2 = new G4SubtractionSolid("TinCanberraI2"+name,TinCanberraI,CopperCanberraI);
  G4SubtractionSolid *LeadCanberraI2 = new G4SubtractionSolid("LeadCanberraI2"+name,LeadCanberraI,TinCanberraI);
  G4SubtractionSolid *SteelCanberraI2 = new G4SubtractionSolid("SteelCanberraI2"+name,SteelCanberraI,LeadCanberraI);

  G4Tubs *bottomHoleCanberra = new G4Tubs("bottomHoleCanberra"+name,0.,(3.+(2./8.))*in/2.,(1.6*mm + 1.0*mm + 4.*in + 9.5*mm)/2.,0.,360*deg);
  G4ThreeVector bHC_Pos(0,0,-1*(1.6*mm + 1.0*mm + 4.*in + 9.5*mm)/2. - 203.*mm);

  G4SubtractionSolid *CopperCanberra = new G4SubtractionSolid("CopperCanberra"+name,CopperCanberraI2,bottomHoleCanberra,0,bHC_Pos);
  G4SubtractionSolid *TinCanberra = new G4SubtractionSolid("TinCanberra"+name,TinCanberraI2,bottomHoleCanberra,0,bHC_Pos);
  G4SubtractionSolid *LeadCanberra = new G4SubtractionSolid("LeadCanberra"+name,LeadCanberraI2,bottomHoleCanberra,0,bHC_Pos);
  G4SubtractionSolid *SteelCanberra = new G4SubtractionSolid("SteelCanberra"+name,SteelCanberraI2,bottomHoleCanberra,0,bHC_Pos);

  G4ThreeVector ShieldCanberra_Pos(0,0,645*mm/2. + tubeDPos.z() + tubeDlen/2.);


  G4LogicalVolume *fLogicDetTubeA = new G4LogicalVolume(detTubeA, mats->GetAluminum(),"detTubeA_log"+name);
  G4LogicalVolume *fLogicDetTubeB = new G4LogicalVolume(detTubeB, mats->GetAluminum(),"detTubeB_log"+name);
  G4LogicalVolume *fLogicDetTubeC = new G4LogicalVolume(detTubeC, mats->GetStainlessSteel(),"detTubeC_log"+name);
  G4LogicalVolume *fLogicDetTubeD = new G4LogicalVolume(detTubeD, mats->GetStainlessSteel(),"detTubeD_log"+name);
  G4LogicalVolume *fLogicDetLiN2D = new G4LogicalVolume(liqN2_TubeD,mats->GetLiN2(),"liqN2_TubeD_log"+name);
  //G4LogicalVolume *fLogicDetLead = new G4LogicalVolume(detLeadShield,mats->GetShieldLead(),"detLeadShield_log"+name);
  G4LogicalVolume *fLogicDetCSCu = new G4LogicalVolume(CopperCanberra,mats->GetCopper(),"CopperCanberra_log"+name);
//  G4LogicalVolume *fLogicDetCSTin = new G4LogicalVolume(TinCanberra,mats->GetTin(),"TinCanberra_log"+name);
  G4LogicalVolume *fLogicDetCSTin = new G4LogicalVolume(TinCanberra,mats->GetCadmium(),"TinCanberra_log"+name);  
  G4LogicalVolume *fLogicDetCSPb = new G4LogicalVolume(LeadCanberra,mats->GetShieldLead(),"LeadCanberra_log"+name);
  G4LogicalVolume *fLogicDetCSSteel = new G4LogicalVolume(SteelCanberra,mats->GetStainlessSteel(),"SteelCanberra_log"+name);

  fullDet->AddPlacedVolume(fLogicDetTubeA,tubeAPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetTubeB,tubeBPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetTubeC,tubeCPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetTubeD,tubeDPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetLiN2D,tubeDPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDet,zeroPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicCasing,cuHousingPos,zeroRot);
  //fullDet->AddPlacedVolume(fLogicDetLead,detLead_Pos,zeroRot);
  //fullDet->AddPlacedVolume(fLogicDetCSCu,ShieldCanberra_Pos,zeroRot);
  //fullDet->AddPlacedVolume(fLogicDetCSTin,ShieldCanberra_Pos,zeroRot);
  //fullDet->AddPlacedVolume(fLogicDetCSPb,ShieldCanberra_Pos,zeroRot);
  //fullDet->AddPlacedVolume(fLogicDetCSSteel,ShieldCanberra_Pos,zeroRot);

  // visualization stuff
  fLogicDetTubeA->SetVisAttributes(aVisAttAlTubes);
  fLogicDetTubeB->SetVisAttributes(aVisAttAlTubes);
  fLogicDetTubeC->SetVisAttributes(aVisAttStainless);
  fLogicDetTubeD->SetVisAttributes(aVisAttStainless);
  fLogicDetLiN2D->SetVisAttributes(aVisAttLiqN2);
  fLogicDet->SetVisAttributes(aVisAttWater);
  fLogicDetCSCu->SetVisAttributes(aVisAttCu);
  fLogicDetCSTin->SetVisAttributes(aVisAttAlTubes);
  fLogicDetCSPb->SetVisAttributes(aVisAttLead);
  fLogicDetCSSteel->SetVisAttributes(aVisAttStainless);


  return fullDet;
}

G4VPhysicalVolume* GeometryConstruction::Construct()
{

  // some useful stuff
  //G4double in = 2.54*cm;
  //G4bool fCheckOverlaps = true;
  //const MINERMaterials * mats = MINERMaterials::GetInstance();


  // ok, let's make a real attempt at the thermal column geometry.  I'm going to attempt 
  // to paramertize this as much as possible to facilitate quick changes.  We're going to set 
  // the origin of the coordinate system at the center of the "circular" part of the pool 
  // for now to make defining the geometry easier.

  // z is up
  // x is cental axis of thermal column
  // y is across face of thermal column


  // mixed units suck but them's whats we gots

  G4double worldY = 10*m;
  G4double worldX = 10*m;
  G4double worldZ = 6*m;

  // define the world volume
  G4Box * universe_s  = new G4Box("universe_s", worldX/2.,worldY/2.,worldZ/2.);
  G4LogicalVolume * universe_log = new G4LogicalVolume(universe_s,mats->GetAir(),"universe_L",0,0,0);


  G4double bioShield_thick = (19.625 + 49.6875) *in ; //1.7*m; // this includes SS plating which is why you'll see it subtracted out below
  G4double water_rad = 12.*4.5*in;  // from Sean's talk at TAMU CNS workshop
  G4double floor_to_TC_out = 48*in; // distance from ground floor to center of thermal column, ESTIMATE
  G4double floor_to_TC_pool = 60*in; // distance from pool floor to center of thermal column 
  G4double floor_thick_out = 0.5*m; // reasonable?
  G4double floor_thick_pool = (floor_thick_out-(floor_to_TC_pool-floor_to_TC_out)); // this includes SS plating which is why you'll see it subtracted out below
  G4double thermalC_height_big = 1.24*m; // larger part of thermal column, square face side
//  G4double thermalC_height_small = 1.0*m; // smaller part of thermal column, square face side
  G4double thermalC_height_small = 39.5*in; // smaller part of thermal column, square face side
  G4double thermalC_width_small = 41*in; // smaller part of thermal column, square face side
  G4double thermalC_length_big = 49.6875*in;

//  G4double thermalC_length_small = (19.625+8.875+25.8125)*in;
  G4double thermalC_length_small = (19.625+8.875+25.8125+1.)*in;

  G4double bioShield_flat_width = thermalC_height_big + 1*m; // this is the width of the "flat" part of the bio shield where the thermal column is ESTIMATE
  G4double aluminumPlate_thick = 0.5*in; // THIS IS A GUESS
  G4double SSPlate_thick = 0.5*in; // THIS IS A GUESS
  G4double poolC_to_TC = water_rad - (25.8125+8.875)*in; // this is the distance from the center of the circular part of the pool to the back of the thermal column
  G4double gb_depth = 26.5625*in; // depth of graphite block
  G4double gb_width = 29.5*in; // width of graphite block ESTIMATE
  G4double gb_height = 29.5*in; // height of graphite block
  G4double largeOverburden_length = 19.625*in; // this is the length of the part of the TC with the largest overburden
  G4double TC_water_length = 25.8125*in; // this is the length of the part of the TC that extends into the water (ignoring the curvature of the wall)
                                         //  and also = the length of the extension that's aluminum instead of SS


  G4cout << water_rad+bioShield_thick-thermalC_length_big-thermalC_length_small-SSPlate_thick << 0. << -(worldZ/2.)+floor_thick_out+floor_to_TC_out << G4endl;
  //G4cout << water_rad+bioShield_thick-thermalC_length_big-((8.875+25.8125)*in)/2. << "  " <<  0. << "  " << -(worldZ/2.)+floor_thick_out+floor_to_TC_out + 0.5*thermalC_height_small << G4endl;
  //G4cout << water_rad+bioShield_thick-thermalC_length_big-((8.875+25.8125)*in)/2. << "  " <<  0. << "  " << -(worldZ/2.)+floor_thick_out+floor_to_TC_out - 0.5*thermalC_height_small << G4endl;
  //G4cout << water_rad+bioShield_thick-thermalC_length_big-((8.875+25.8125)*in)/2. << "  " <<  0.5*thermalC_height_small  << "  " << -(worldZ/2.)+floor_thick_out+floor_to_TC_out << G4endl;

  //G4cout << poolC_to_TC - 3.86*m << "  " << 0. << "  " << -(worldZ/2.)+floor_thick_out+floor_to_TC_out << G4endl;

  // floor
  // fill the entire world volume x and y with floor
  // we'll subtract out the difference in pool and outer floor heights to get floor in the pool
  G4Box *floor_out = new G4Box("floor_out",worldX/2.,worldY/2.,floor_thick_out/2.);
  G4Box *floor_rec_pool = new G4Box("floor_rec_pool",worldX/4.,water_rad + SSPlate_thick,(floor_to_TC_pool-floor_to_TC_out + SSPlate_thick)/2.);
  G4Tubs *floor_tubs_pool = new G4Tubs("floor_tubs_pool",0.,water_rad + SSPlate_thick,(floor_to_TC_pool-floor_to_TC_out + SSPlate_thick)/2.,-pi/2.,pi);
  G4ThreeVector poolFloorTransX(-worldX/4.,0,0);
  G4UnionSolid *floorI = new G4UnionSolid("floorI",floor_tubs_pool,floor_rec_pool,0,poolFloorTransX);
  G4ThreeVector poolFloorTransZ(0,0,(floor_thick_out - (floor_to_TC_pool-floor_to_TC_out) - SSPlate_thick)/2.); // check this
  G4SubtractionSolid *floor = new G4SubtractionSolid("floor",floor_out,floorI,0,poolFloorTransZ);
  G4LogicalVolume *floor_log = new G4LogicalVolume(floor,mats->GetHDConcrete() ,"floor_log");
  G4ThreeVector floorPlacement(0,0,-(worldZ/2.)+(floor_thick_out/2.));


  // graphite block in the ppol
  G4Box *graphite_block = new G4Box("graphite_block", gb_depth/2.,gb_width/2.,gb_height/2.);
  G4ThreeVector graphite_Placement(water_rad+bioShield_thick-thermalC_length_big-thermalC_length_small- (gb_depth/2.),0., -(worldZ/2.) + floor_thick_out + floor_to_TC_out);
  G4LogicalVolume *graphite_log = new G4LogicalVolume(graphite_block,mats->GetGraphite(),"graphite_lot");

  //water 
  // all the way up for now
  // two parts: semicircle and rectangular
  G4Box *pool_rec = new G4Box("pool_rec",worldX/4.,water_rad,(worldZ/2.) - floor_thick_pool/2.);
  G4Tubs *pool_tubs = new G4Tubs("pool_tubs",0.,water_rad,(worldZ/2.) - floor_thick_pool/2.,-pi/2.,pi);
  G4UnionSolid *poolI = new G4UnionSolid("poolI",pool_tubs,pool_rec,0,poolFloorTransX);
  G4Box *TC_pool_cutout = new G4Box("TC_pool_cutout",(thermalC_length_small)/2.,(thermalC_width_small)/2. + SSPlate_thick,(thermalC_height_small)/2. + SSPlate_thick);
  G4ThreeVector TC_in_poolPlacement(water_rad+bioShield_thick-thermalC_length_big-(thermalC_length_small/2.), 0, -((worldZ/2.) - floor_thick_pool/2.) + floor_to_TC_pool); // move to bottom of floor then up to proper height
  G4ThreeVector GB_in_poolPlacement(water_rad+bioShield_thick-thermalC_length_big-thermalC_length_small-gb_depth/2., 0, -((worldZ/2.) - floor_thick_pool/2.) + floor_to_TC_pool);
  G4SubtractionSolid *poolI2 = new G4SubtractionSolid("poolI2",poolI,TC_pool_cutout,0,TC_in_poolPlacement);
  G4SubtractionSolid *pool = new G4SubtractionSolid("pool",poolI2,graphite_block,0,GB_in_poolPlacement);
  G4ThreeVector poolPlacement(0,0,floor_thick_pool/2.);
  G4LogicalVolume *pool_log = new G4LogicalVolume(pool,mats->GetWater(),"pool_log");
  // could add support legs, probably a third order effect


  //bioShield (aka the HDConcrete wall between us and the water
  // all the way up for now
  // three parts: semicircle and two boxes
  // have to leave room for the SS plates
  // Need to add upper level "shelf" still
  // The front part with the TC entrance is actually flat, need to shave it off so the door can be flush
  G4Box *bioShield_rec = new G4Box("bioShield_rec",worldX/4.,(bioShield_thick-SSPlate_thick)/2.,(worldZ/2.)-(floor_thick_out/2.) );
  G4Tubs *bioShield_tubs = new G4Tubs("bioShield_tubs",water_rad+SSPlate_thick,water_rad+bioShield_thick,(worldZ/2.)-(floor_thick_out/2.),-pi/2.,pi);
  G4ThreeVector bioShield1Trans(-worldX/4.,water_rad+((bioShield_thick-SSPlate_thick)/2.)+SSPlate_thick,0.);
  G4ThreeVector bioShield2Trans(-worldX/4.,-1*(water_rad+((bioShield_thick-SSPlate_thick)/2.)+SSPlate_thick),0.);
  G4UnionSolid *bioShieldI1 = new G4UnionSolid("bioShieldI1",bioShield_tubs,bioShield_rec,0,bioShield1Trans);
  G4UnionSolid *bioShieldI2 = new G4UnionSolid("bioShieldI2",bioShieldI1,bioShield_rec,0,bioShield2Trans);
  G4Box *TC_BS_cutout1 = new G4Box("TC_BS_cutout1",(thermalC_length_small)/2.,(thermalC_width_small)/2. + SSPlate_thick,(thermalC_height_small)/2. + SSPlate_thick);
  G4Box *TC_BS_cutout2 = new G4Box("TC_BS_cutout2",(thermalC_length_big+aluminumPlate_thick)/2.,(thermalC_height_big)/2. + aluminumPlate_thick,(thermalC_height_big)/2. + aluminumPlate_thick);
  //there's still some HDC left because of the curvature of the wall so cut that out as well
  G4Box *TC_BS_cutout3 = new G4Box("TC_BS_cutout3",(thermalC_length_big+aluminumPlate_thick)/2.,(thermalC_height_big)/2. + aluminumPlate_thick,(thermalC_height_big)/2. +  aluminumPlate_thick);
  G4ThreeVector TCsmall_in_BSPlacement(water_rad+bioShield_thick-thermalC_length_big-(thermalC_length_small/2.),0,-((worldZ/2.)-(floor_thick_out/2.))+floor_to_TC_out);
  G4ThreeVector TCbig_in_BSPlacement(water_rad+bioShield_thick-(thermalC_length_big/2.),0,-((worldZ/2.)-(floor_thick_out/2.))+floor_to_TC_out);
  G4ThreeVector TCleftover_in_BSPlacement(water_rad+bioShield_thick,0,-((worldZ/2.)-(floor_thick_out/2.))+floor_to_TC_out);
  G4SubtractionSolid * TC_in_BS1 = new G4SubtractionSolid("TC_in_BS1",bioShieldI2,TC_BS_cutout1,0,TCsmall_in_BSPlacement);
  G4SubtractionSolid * TC_in_BS2 = new G4SubtractionSolid("TC_in_BS2",TC_in_BS1,TC_BS_cutout2,0,TCbig_in_BSPlacement);
  G4SubtractionSolid * bioShieldI3 = new G4SubtractionSolid("bioShieldI3",TC_in_BS2,TC_BS_cutout3,0,TCleftover_in_BSPlacement);
  G4double sagitta = (water_rad + bioShield_thick) - sqrt(pow((water_rad + bioShield_thick),2) - pow((bioShield_flat_width/2.),2));
  G4Box *bioShield_front_add = new G4Box("bioShield_front_add",sagitta/2.,bioShield_flat_width/2.,(worldZ/2.)-(floor_thick_out/2.) );
  G4ThreeVector TCbig_in_BSFlatPlacement(0,0,-((worldZ/2.)-(floor_thick_out/2.))+floor_to_TC_out);
  G4SubtractionSolid * bioShield_front_sub = new G4SubtractionSolid("bioShield_front_sub",bioShield_front_add,TC_BS_cutout2,0,TCbig_in_BSFlatPlacement);
  G4ThreeVector Flat_BSPlacement(water_rad + bioShield_thick - sagitta/2.,0,0);
  G4UnionSolid *bioShield = new G4UnionSolid("bioShield",bioShieldI3,bioShield_front_sub,0,Flat_BSPlacement);
  G4ThreeVector bioShieldPlacement(0,0,(floor_thick_out)/2.);
  G4LogicalVolume *bioShield_log = new G4LogicalVolume(bioShield,mats->GetHDConcrete(),"bioShield_log");


  // SS pool liner
  // five parts: semicircle and two boxes for sides
  // semicircle and box for bottom
  G4Box *poolSS_rec = new G4Box("poolSS_rec",worldX/4.,SSPlate_thick/2.,(worldZ/2.) - floor_thick_pool/2. );
  G4Tubs *poolSS_tubs = new G4Tubs("poolSS_tubs",water_rad,water_rad+SSPlate_thick,(worldZ/2.) - floor_thick_pool/2.,-pi/2.,pi);
  G4ThreeVector poolSS1Trans(-worldX/4.,water_rad+SSPlate_thick/2.,0.);
  G4ThreeVector poolSS2Trans(-worldX/4.,-1*(water_rad+SSPlate_thick/2.),0.);
  G4UnionSolid *poolSSI1 = new G4UnionSolid("poolSSI1",poolSS_tubs,poolSS_rec,0,poolSS1Trans);
  G4UnionSolid *poolSSI2 = new G4UnionSolid("poolSSI2",poolSSI1,poolSS_rec,0,poolSS2Trans);
  G4ThreeVector TC_in_SSPlacement(water_rad+bioShield_thick-thermalC_length_big-(thermalC_length_small/2), 0, -((worldZ/2.) - floor_thick_pool/2.) + floor_to_TC_pool); // move to bottom of floor then up to proper height
  G4SubtractionSolid * poolSSI3 = new G4SubtractionSolid("poolSSI3",poolSSI2,TC_pool_cutout,0,TC_in_SSPlacement);
  G4Tubs *poolSS_bot_tubs = new G4Tubs("poolSS_bot_tubs",0,water_rad+SSPlate_thick,SSPlate_thick/2.,-pi/2.,pi);
  G4Box *poolSS_bot_rec = new G4Box("poolSS_bot_rec",worldX/4.,(water_rad+SSPlate_thick)/2.,SSPlate_thick/2.);
  G4ThreeVector poolSS3Trans(-worldX/4.,0.,-((worldZ/2.) - floor_thick_pool/2.)-SSPlate_thick/2.);
  G4ThreeVector poolSS4Trans(0.,0.,-((worldZ/2.) - floor_thick_pool/2.)-SSPlate_thick/2.);
  G4UnionSolid *poolSSI4 = new G4UnionSolid("poolSSI4",poolSSI3,poolSS_bot_rec,0,poolSS3Trans);
  G4UnionSolid *poolSS = new G4UnionSolid("poolSS",poolSSI4,poolSS_bot_tubs,0,poolSS4Trans);
  G4ThreeVector poolSSPlacement(0,0,floor_thick_pool/2.);
  G4LogicalVolume *poolSS_log = new G4LogicalVolume(poolSS,mats->GetStainlessSteel(),"poolSS_log");

  // SS TC liner 
  // SS only in the small box part
  // there's some SS that goes into the HDC that I'm not including yet
  // Fixed March 9: TC extension is actually aluminum
  G4Box *TC_SS_box = new G4Box("TC_SS_box", (thermalC_length_small-TC_water_length)/2.,(thermalC_width_small)/2. + SSPlate_thick,(thermalC_height_small)/2. + SSPlate_thick);
  G4Box *TC_SS_cutout = new G4Box("TC_SS_cutout", thermalC_length_small/2.,thermalC_width_small/2.,thermalC_height_small/2.);
  G4ThreeVector TC_SS_placement(water_rad+bioShield_thick-thermalC_length_big- (thermalC_length_small-TC_water_length)/2.,0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  G4SubtractionSolid *TC_SS = new G4SubtractionSolid("TC_SS",TC_SS_box,TC_SS_cutout);

  G4Box *TC_Al_extI = new G4Box("TC_Al_extI", TC_water_length/2.,(thermalC_width_small)/2. + SSPlate_thick,(thermalC_height_small)/2. + SSPlate_thick);
  G4Box *TC_Al_cut = new G4Box("TC_Al_cut", TC_water_length/2.,(thermalC_width_small)/2. + SSPlate_thick,(thermalC_height_small)/2. + SSPlate_thick);
  G4ThreeVector TC_Alext_cutout_Trans(SSPlate_thick,0.,0.);
  G4SubtractionSolid *TC_Al_ext = new G4SubtractionSolid("TC_Al_ext",TC_Al_extI,TC_Al_cut,0,TC_Alext_cutout_Trans);
  G4ThreeVector TC_Alext_placement(water_rad+bioShield_thick-thermalC_length_big - (thermalC_length_small-TC_water_length) - (TC_water_length/2.),0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);

  G4LogicalVolume *TC_SS_log =  new G4LogicalVolume(TC_SS,mats->GetStainlessSteel(),"TC_SS_log");
  G4LogicalVolume *TC_Alext_log =  new G4LogicalVolume(TC_Al_ext,mats->GetAluminum(),"TC_Alext_log");

  // Al TC liner 
  // Al only in the big box part
  G4Box *TC_Al_box = new G4Box("TC_Al_box", thermalC_length_big/2.,(thermalC_height_big)/2. + aluminumPlate_thick,(thermalC_height_big)/2. + aluminumPlate_thick);
  G4Box *TC_Al_cutout = new G4Box("TC_Al_cutout", thermalC_length_big/2.,thermalC_height_big/2.,thermalC_height_big/2.);
  G4ThreeVector TC_Al_cutout_Trans(aluminumPlate_thick,0.,0.);
  G4SubtractionSolid *TC_AlI = new G4SubtractionSolid("TC_AlI",TC_Al_box,TC_Al_cutout,0,TC_Al_cutout_Trans);
  G4ThreeVector TC_Al_cutout_Trans2(-aluminumPlate_thick,0.,0.);
  G4SubtractionSolid *TC_Al = new G4SubtractionSolid("TC_Al",TC_AlI,TC_SS_cutout,0,TC_Al_cutout_Trans2);
  G4ThreeVector TC_Al_placement(water_rad+bioShield_thick-thermalC_length_big/2.,0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  G4LogicalVolume *TC_Al_log =  new G4LogicalVolume(TC_Al,mats->GetAluminum(),"TC_Al_log");

  // the movable door
  // keep it simple for now
  // 3/8" boral liner on front NOT IN YET
  // 4" lead
  // 2ft - 4" HDConcrete
  // 1/16" steel on back
  G4Box *doorLiner = new G4Box( "doorLiner",((3./8.)*in)/2., ((6*12 + 8)*in)/2., (6.*12.*in)/2.);
  G4Box *doorLead = new G4Box("doorLead", (4.*in)/2., ((6*12 + 8)*in)/2., (6.*12.*in)/2.);
  G4Box *doorConcrete = new G4Box("doorConcrete", (20.*in)/2., ((6*12 + 8)*in)/2., (6.*12.*in)/2.);
  G4Box *doorSteel = new G4Box("doorSteel", ((1./16.)*in)/2., ((6*12 + 8)*in)/2., (6.*12.*in)/2.);
  G4ThreeVector doorLinerPos(((3./8.)*in)/2.,0.,0.);
  G4ThreeVector doorLeadPos(((3./8.)*in) + (4.*in)/2.,0.,0.);
  G4ThreeVector doorConcretePos(((3./8.)*in) + (4.*in) + (20.*in)/2.,0.,0.);
  G4ThreeVector doorSteelPos(((3./8.)*in) + (4.*in) + (20.*in) + ((1./16.)*in)/2.,0.,0.);
  G4LogicalVolume *door_Liner_log = new G4LogicalVolume(doorLiner,mats->GetBoraflex(),"door_Liner_log");
  G4LogicalVolume *door_Lead_log = new G4LogicalVolume(doorLead,mats->GetShieldLead() ,"door_Lead_log");
  G4LogicalVolume *door_Concrete_log = new G4LogicalVolume(doorConcrete,mats->GetHDConcrete() ,"door_Concrete_log");
  G4LogicalVolume *door_Steel_log = new G4LogicalVolume(doorSteel,mats->GetStainlessSteel() ,"door_Steel_log");


  G4AssemblyVolume *door = new G4AssemblyVolume();
  //G4RotationMatrix *zeroRot = new G4RotationMatrix;
  //G4ThreeVector zeroPos;
  door->AddPlacedVolume(door_Liner_log,doorLinerPos,zeroRot);
  door->AddPlacedVolume(door_Lead_log,doorLeadPos,zeroRot);
  door->AddPlacedVolume(door_Concrete_log,doorConcretePos,zeroRot);
  door->AddPlacedVolume(door_Steel_log,doorSteelPos,zeroRot);
  G4ThreeVector posDoor(water_rad+bioShield_thick,0.,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  // move door back for now
  //G4ThreeVector posDoor(water_rad+bioShield_thick + 3.*12.*in,0.,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);


  //G4AssemblyVolume *waterBrick = new G4AssemblyVolume();
  // Define a water brick (does not yet include side divets, corners, top groove, lid)
  G4double brick_x = 22.86*cm;
  G4double brick_y = 45.72*cm;
  G4double brick_z = 15.24*cm;
  G4double brick_thickness = 0.238*cm;

  G4Box *Outer = new G4Box("Outer", brick_x/2, brick_y/2, brick_z/2);
  G4Box * Inner = new G4Box("Inner", (brick_x-2*brick_thickness)/2, (brick_y-2*brick_thickness)/2,(brick_z-2*brick_thickness)/2);
  G4SubtractionSolid *hollowBox = new G4SubtractionSolid("Hollow Box",Outer, Inner);

  G4Box * inner_water = new G4Box("Inner Water",  (brick_x-2*brick_thickness)/2, (brick_y-2*brick_thickness)/2,(brick_z-2*brick_thickness)/2);
  G4Tubs * hole_out = new G4Tubs("hole_out", 0, 2*cm + brick_thickness, brick_z/2, 0, 360);
  G4Tubs * poly_hole = new G4Tubs("poly_hole", 2*cm, 2*cm + brick_thickness, brick_z/2, 0, 360);


  G4ThreeVector zTrans(0,8*cm,0*cm);
  G4ThreeVector NzTrans(0,-8*cm,0*cm);

  G4SubtractionSolid *Shell = new G4SubtractionSolid("Shell", hollowBox, hole_out,0, zTrans);
  G4SubtractionSolid *Shell2 = new G4SubtractionSolid("Shell2", Shell, hole_out,0, NzTrans);
  G4UnionSolid *Shell3 = new G4UnionSolid("Shell3",Shell2,poly_hole,0,zTrans);
  G4UnionSolid *Shell4 = new G4UnionSolid("Shell4",Shell3,poly_hole,0,NzTrans);


  G4SubtractionSolid *filling = new G4SubtractionSolid("Filling", inner_water,  hole_out,0, zTrans);
  G4SubtractionSolid *filling2 = new G4SubtractionSolid("Filling2", filling,  hole_out,0, NzTrans);

  G4LogicalVolume * logic_shell = new G4LogicalVolume(Shell4, mats->GetHDPoly(), "Shell");
  G4LogicalVolume * logic_filling = new G4LogicalVolume(filling2, mats->GetWater(), "Filling");
  //waterBrick->AddPlacedVolume(logic_shell,zeroPos,zeroRot);
  //waterBrick->AddPlacedVolume(logic_filling,zeroPos,zeroRot);


  // ************************************************** //
  // ************************************************** //
  // ***********   WaterBrick Setup V1  *************** //
  // ************************************************** //
  // ************************************************** //

  G4double wgGap = 0.0*cm; // potential gap between all bricks
  G4AssemblyVolume *WBsubstruct1 = new G4AssemblyVolume();
  G4ThreeVector   brick1_place(0.5*(brick_x + wgGap),0.,0.);
  G4ThreeVector   brick2_place(-0.5*(brick_x + wgGap),0.,0.);
  G4ThreeVector   brick3_place(0.,(0.5*brick_y + 0.5*brick_x + 0.5*wgGap),0.);
  G4ThreeVector   brick4_place(0.,-1*(0.5*brick_y + 0.5*brick_x + 0.5*wgGap),0.);
//  G4ThreeVector   brick3_place(0.,0.5*(brick_y + wgGap) + 0.5*(brick_x + wgGap),0.);
//  G4ThreeVector   brick4_place(0.,-1*(0.5*(brick_y + wgGap) + 0.5*(brick_x + wgGap)),0.);
  G4RotationMatrix* RotBrick = new G4RotationMatrix;
  RotBrick->rotateZ(pi/2.);
  RotBrick->rotateX(0.);
  RotBrick->rotateY(0.);
  WBsubstruct1->AddPlacedVolume(logic_shell,brick1_place,zeroRot);
  WBsubstruct1->AddPlacedVolume(logic_shell,brick2_place,zeroRot);
  WBsubstruct1->AddPlacedVolume(logic_shell,brick3_place,RotBrick);
  WBsubstruct1->AddPlacedVolume(logic_shell,brick4_place,RotBrick);
  WBsubstruct1->AddPlacedVolume(logic_filling,brick1_place,zeroRot);
  WBsubstruct1->AddPlacedVolume(logic_filling,brick2_place,zeroRot);
  WBsubstruct1->AddPlacedVolume(logic_filling,brick3_place,RotBrick);
  WBsubstruct1->AddPlacedVolume(logic_filling,brick4_place,RotBrick);

  G4AssemblyVolume *WBsubstruct2 = new G4AssemblyVolume();
  G4ThreeVector   brick5_place(0.5*(brick_x + wgGap),0.5*(brick_y + wgGap),0.);
  G4ThreeVector   brick6_place(-0.5*(brick_x + wgGap),(0.5*(brick_y + wgGap)),0.);
  G4ThreeVector   brick7_place(0.5*(brick_x + wgGap),-0.5*(brick_y + wgGap),0.);
  G4ThreeVector   brick8_place(-0.5*(brick_x + wgGap),-1*(0.5*(brick_y + wgGap)),0.);
  WBsubstruct2->AddPlacedVolume(logic_shell,brick5_place,zeroRot);
  WBsubstruct2->AddPlacedVolume(logic_shell,brick6_place,zeroRot);
  WBsubstruct2->AddPlacedVolume(logic_shell,brick7_place,zeroRot);
  WBsubstruct2->AddPlacedVolume(logic_shell,brick8_place,zeroRot);
  WBsubstruct2->AddPlacedVolume(logic_filling,brick5_place,zeroRot);
  WBsubstruct2->AddPlacedVolume(logic_filling,brick6_place,zeroRot);
  WBsubstruct2->AddPlacedVolume(logic_filling,brick7_place,zeroRot);
  WBsubstruct2->AddPlacedVolume(logic_filling,brick8_place,zeroRot);


  // Lead Wall
  G4double leadWall_x = 4.*in;
  G4double leadWall_y = thermalC_width_small;
  G4double leadWall_z = thermalC_height_small;

  G4Box *leadWall = new G4Box("leadWall", leadWall_x/2, leadWall_y/2, leadWall_z/2);
  G4LogicalVolume * leadWall_log = new G4LogicalVolume(leadWall, mats->GetShieldLead(), "leadWall_log");

  G4Box *leadWall3 = new G4Box("leadWall3", leadWall_x/2, leadWall_y/2, leadWall_z/2);
  G4LogicalVolume * leadWall3_log = new G4LogicalVolume(leadWall, mats->GetShieldLead(), "leadWall3_log");


  G4double leadWall2_x = 4.*in;
  G4double leadWall2_y = thermalC_height_big;
  G4double leadWall2_z = thermalC_height_big;

  G4Box *leadWall2 = new G4Box("leadWall2", leadWall2_x/2, leadWall2_y/2, leadWall2_z/2);
  G4LogicalVolume * leadWall2_log = new G4LogicalVolume(leadWall2, mats->GetShieldLead(), "leadWall2_log");


  // poly sheets
  G4double polySheet1_y = 1.0*in;
  G4double polySheet1_z = 38.25*in;
  G4double polySheet1_x = 18.0*in;

  G4double polySheet2_y = 36.0*in;
  G4double polySheet2_z = 1.0*in;
  G4double polySheet2_x = 18.0*in;

  G4double polySheet3_y = 40.5*in;
  G4double polySheet3_z = 1.0*in;
  G4double polySheet3_x = 18.0*in;

  G4Box *poly1 = new G4Box("poly1", polySheet1_x/2, polySheet1_y/2, polySheet1_z/2);
  G4Box *poly2 = new G4Box("poly2", polySheet2_x/2, polySheet2_y/2, polySheet2_z/2);
  G4Box *poly3 = new G4Box("poly3", polySheet3_x/2, polySheet3_y/2, polySheet3_z/2);

  // neoprene/rubber sheets
  G4double rubberSheet1_y = 1.0*in;
  G4double rubberSheet1_z = 38.25*in;
  G4double rubberSheet1_x = 18.0*in;

  G4double rubberSheet2_y = 36.0*in;
  G4double rubberSheet2_z = 0.5*in;
  G4double rubberSheet2_x = 18.0*in;

  G4double rubberSheet4_y = 2.0*in;
  G4double rubberSheet4_z = 0.25*in;
  G4double rubberSheet4_x = 18.0*in;

  G4double rubberSheet5_y = 0.5*in;
  G4double rubberSheet5_z = 1.0*in;
  G4double rubberSheet5_x = 18.0*in;

  G4Box *rubber1 = new G4Box("rubber1", rubberSheet1_x/2, rubberSheet1_y/2, rubberSheet1_z/2);
  G4Box *rubber2 = new G4Box("rubber2", rubberSheet2_x/2, rubberSheet2_y/2, rubberSheet2_z/2);
  G4Box *rubber4 = new G4Box("rubber4", rubberSheet4_x/2, rubberSheet4_y/2, rubberSheet4_z/2);
  G4Box *rubber5 = new G4Box("rubber5", rubberSheet5_x/2, rubberSheet5_y/2, rubberSheet5_z/2);


  // Layer 1 logical volumes
  G4LogicalVolume * poly_left_1_layer1 = new G4LogicalVolume(poly1, mats->GetHDPoly(), "poly_left_1_layer1");
  G4LogicalVolume * poly_left_2_layer1 = new G4LogicalVolume(poly1, mats->GetHDPoly(), "poly_left_2_layer1");
  G4LogicalVolume * poly_right_1_layer1 = new G4LogicalVolume(poly1, mats->GetHDPoly(), "poly_right_1_layer1");
  G4LogicalVolume * poly_right_2_layer1 = new G4LogicalVolume(poly1, mats->GetHDPoly(), "poly_right_2_layer1");
  G4LogicalVolume * poly_top_1_layer1 = new G4LogicalVolume(poly2, mats->GetHDPoly(), "poly_top_1_layer1");
  G4LogicalVolume * poly_top_2_layer1 = new G4LogicalVolume(poly2, mats->GetHDPoly(), "poly_top_2_layer1");
  G4LogicalVolume * poly_top_3_layer1 = new G4LogicalVolume(poly3, mats->GetHDPoly(), "poly_top_3_layer1");

  G4LogicalVolume * rub_left_1_layer1 = new G4LogicalVolume(rubber1, mats->GetNeopreneBlend(), "rubber_left_1_layer1");
  G4LogicalVolume * rub_right_1_layer1 = new G4LogicalVolume(rubber5, mats->GetNeopreneBlend(), "rubber_right_1_layer1");
  G4LogicalVolume * rub_top_1_layer1 = new G4LogicalVolume(rubber4, mats->GetNeopreneBlend(), "rubber_top_1_layer1");
  G4LogicalVolume * rub_top_2_layer1 = new G4LogicalVolume(rubber4, mats->GetNeopreneBlend(), "rubber_top_2_layer1");
  G4LogicalVolume * rub_top_3_layer1 = new G4LogicalVolume(rubber2, mats->GetNeopreneBlend(), "rubber_top_3_layer1");


  // Layer 2 logical volumes
  G4LogicalVolume * poly_left_1_layer2 = new G4LogicalVolume(poly1, mats->GetHDPoly(), "poly_left_1_layer2");
  G4LogicalVolume * poly_left_2_layer2 = new G4LogicalVolume(poly1, mats->GetHDPoly(), "poly_left_2_layer2");
  G4LogicalVolume * poly_right_1_layer2 = new G4LogicalVolume(poly1, mats->GetHDPoly(), "poly_right_1_layer2");
  G4LogicalVolume * poly_right_2_layer2 = new G4LogicalVolume(poly1, mats->GetHDPoly(), "poly_right_2_layer2");
  G4LogicalVolume * poly_bottom_1_layer2 = new G4LogicalVolume(poly2, mats->GetHDPoly(), "poly_bottom_1_layer2");
  G4LogicalVolume * poly_top_1_layer2 = new G4LogicalVolume(poly2, mats->GetHDPoly(), "poly_top_1_layer2");
  G4LogicalVolume * poly_top_2_layer2 = new G4LogicalVolume(poly3, mats->GetHDPoly(), "poly_top_2_layer2");
  
  G4LogicalVolume * rub_right_1_layer2 = new G4LogicalVolume(rubber1, mats->GetNeopreneBlend(), "rubber_right_1_layer2");
  G4LogicalVolume * rub_left_1_layer2 = new G4LogicalVolume(rubber5, mats->GetNeopreneBlend(), "rubber_left_1_layer2");
  G4LogicalVolume * rub_top_1_layer2 = new G4LogicalVolume(rubber4, mats->GetNeopreneBlend(), "rubber_top_1_layer2");
  G4LogicalVolume * rub_top_2_layer2 = new G4LogicalVolume(rubber4, mats->GetNeopreneBlend(), "rubber_top_2_layer2");
  G4LogicalVolume * rub_top_3_layer2 = new G4LogicalVolume(rubber2, mats->GetNeopreneBlend(), "rubber_top_3_layer2");

  // Layer 3 logical volumes
  G4LogicalVolume * poly_left_1_layer3 = new G4LogicalVolume(poly1, mats->GetHDPoly(), "poly_left_1_layer3");
  G4LogicalVolume * poly_left_2_layer3 = new G4LogicalVolume(poly1, mats->GetHDPoly(), "poly_left_2_layer3");
  G4LogicalVolume * poly_right_1_layer3 = new G4LogicalVolume(poly1, mats->GetHDPoly(), "poly_right_1_layer3");
  G4LogicalVolume * poly_right_2_layer3 = new G4LogicalVolume(poly1, mats->GetHDPoly(), "poly_right_2_layer3");
  G4LogicalVolume * poly_bottom_1_layer3 = new G4LogicalVolume(poly2, mats->GetHDPoly(), "poly_bottom_1_layer3");
  G4LogicalVolume * poly_bottom_2_layer3 = new G4LogicalVolume(poly2, mats->GetHDPoly(), "poly_bottom_2_layer3");
  G4LogicalVolume * poly_top_1_layer3 = new G4LogicalVolume(poly3, mats->GetHDPoly(), "poly_top_1_layer3");

  G4LogicalVolume * rub_left_1_layer3 = new G4LogicalVolume(rubber1, mats->GetNeopreneBlend(), "rubber_left_1_layer3");
  G4LogicalVolume * rub_right_1_layer3 = new G4LogicalVolume(rubber5, mats->GetNeopreneBlend(), "rubber_right_1_layer3");
  G4LogicalVolume * rub_top_1_layer3 = new G4LogicalVolume(rubber4, mats->GetNeopreneBlend(), "rubber_top_1_layer3");
  G4LogicalVolume * rub_top_2_layer3 = new G4LogicalVolume(rubber4, mats->GetNeopreneBlend(), "rubber_top_2_layer3");
  G4LogicalVolume * rub_top_3_layer3 = new G4LogicalVolume(rubber2, mats->GetNeopreneBlend(), "rubber_top_3_layer3");  

  // placements for everything relative to (0,0,0)
  G4ThreeVector   posLayer1(water_rad+bioShield_thick +SSPlate_thick  - thermalC_length_big - thermalC_length_small + 9.*in,0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  G4ThreeVector   posLayer2(water_rad+bioShield_thick +SSPlate_thick  - thermalC_length_big - thermalC_length_small + 27.*in,0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  G4ThreeVector   posLayer3(water_rad+bioShield_thick +SSPlate_thick  - thermalC_length_big - thermalC_length_small + 53.*in,0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  G4ThreeVector   posLayerLead(water_rad+bioShield_thick +SSPlate_thick  - thermalC_length_big - thermalC_length_small + 38.*in,0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  G4ThreeVector   posLayerLead2(water_rad+bioShield_thick +SSPlate_thick  - thermalC_length_big - thermalC_length_small + 64.*in,0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  G4ThreeVector   posLayerLead3(water_rad+bioShield_thick +SSPlate_thick  - thermalC_length_big - thermalC_length_small + 42.*in,0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);

  // Layer 1 positions
  G4ThreeVector   posBrick1_1(posLayer1.x(),0.5*in ,posLayer1.z() + -1.75*in - 2.5*brick_z );
  G4ThreeVector   posBrick2_1(posLayer1.x(),0.5*in ,posLayer1.z() + -1.75*in - 1.5*brick_z );
  G4ThreeVector   posBrick3_1(posLayer1.x(),0.5*in ,posLayer1.z() + -1.75*in - 0.5*brick_z );
  G4ThreeVector   posBrick4_1(posLayer1.x(),0.5*in ,posLayer1.z() + -1.75*in + 0.5*brick_z );
  G4ThreeVector   posBrick5_1(posLayer1.x(),0.5*in ,posLayer1.z() + -1.75*in + 1.5*brick_z );
  G4ThreeVector   posBrick6_1(posLayer1.x(),0.5*in ,posLayer1.z() + -1.75*in + 2.5*brick_z );
  G4ThreeVector   posPL1_1(posLayer1.x(), 0.5*in - 2.5*in - 18.*in,posLayer1.z() + -0.625*in );
  G4ThreeVector   posPL2_1(posLayer1.x(), 0.5*in - 1.5*in - 18.*in,posLayer1.z() + -0.625*in );
  G4ThreeVector   posPR1_1(posLayer1.x(), 0.5*in + 0.5*in + 18.*in,posLayer1.z() + -0.625*in);
  G4ThreeVector   posPR2_1(posLayer1.x(), 0.5*in + 1.5*in + 18.*in,posLayer1.z() + -0.625*in);
  G4ThreeVector   posPT1_1(posLayer1.x(), 0.5*in,posLayer1.z() + 19.75*in - 2.5*in);
  G4ThreeVector   posPT2_1(posLayer1.x(), 0.5*in,posLayer1.z() + 19.75*in - 1.5*in);
  G4ThreeVector   posPT3_1(posLayer1.x(), -0.25*in,posLayer1.z() + 19.75*in - 0.5*in);
  G4ThreeVector   posRL1_1(posLayer1.x(), 0.5*in - 18.0*in - 0.5*in,posLayer1.z() + -0.625*in );
  G4ThreeVector   posRR1_1(posLayer1.x(), 20.5*in - 0.25*in,posLayer1.z() + 19.75*in - 0.5*in);
  G4ThreeVector   posRT1_1(posLayer1.x(), 20.5*in - 1.0*in ,posLayer1.z() + 19.75*in - 1.125*in);
  G4ThreeVector   posRT2_1(posLayer1.x(), -20.5*in + 1.0*in ,posLayer1.z() + 19.75*in - 1.125*in);
  G4ThreeVector   posRT3_1(posLayer1.x(), 0.5*in,posLayer1.z() + 19.75*in - 3.25*in);

  // Layer 2 positions
  G4ThreeVector   posBrick1_2(posLayer2.x(),-0.5*in ,posLayer1.z() + -0.75*in - 2.5*brick_z );
  G4ThreeVector   posBrick2_2(posLayer2.x(),-0.5*in ,posLayer1.z() + -0.75*in - 1.5*brick_z );
  G4ThreeVector   posBrick3_2(posLayer2.x(),-0.5*in ,posLayer1.z() + -0.75*in - 0.5*brick_z );
  G4ThreeVector   posBrick4_2(posLayer2.x(),-0.5*in ,posLayer1.z() + -0.75*in + 0.5*brick_z );
  G4ThreeVector   posBrick5_2(posLayer2.x(),-0.5*in ,posLayer1.z() + -0.75*in + 1.5*brick_z );
  G4ThreeVector   posBrick6_2(posLayer2.x(),-0.5*in ,posLayer1.z() + -0.75*in + 2.5*brick_z );
  G4ThreeVector   posPL1_2(posLayer2.x(), -0.5*in - 1.5*in - 18.*in,posLayer1.z() + -0.625*in);
  G4ThreeVector   posPL2_2(posLayer2.x(), -0.5*in - 0.5*in - 18.*in,posLayer1.z() + -0.625*in);
  G4ThreeVector   posPR1_2(posLayer2.x(), -0.5*in + 1.5*in + 18.*in,posLayer1.z() + -0.625*in);
  G4ThreeVector   posPR2_2(posLayer2.x(), -0.5*in + 2.5*in + 18.*in,posLayer1.z() + -0.625*in);
  G4ThreeVector   posPB1_2(posLayer2.x(), -0.5*in,posLayer1.z() + -19.75*in + 0.5*in);
  G4ThreeVector   posPT1_2(posLayer2.x(), -0.5*in,posLayer1.z() + 19.75*in - 1.5*in);
  G4ThreeVector   posPT2_2(posLayer2.x(), 0.25*in,posLayer1.z() + 19.75*in - 0.5*in);
  G4ThreeVector   posRL1_2(posLayer2.x(), -20.5*in + 0.25*in,posLayer1.z() + 19.75*in - 0.5*in);
  G4ThreeVector   posRR1_2(posLayer2.x(), -0.5*in + 18.0*in + 0.5*in,posLayer1.z() + -0.625*in);
  G4ThreeVector   posRT1_2(posLayer2.x(), 20.5*in - 1.0*in ,posLayer1.z() + 19.75*in - 1.125*in);
  G4ThreeVector   posRT2_2(posLayer2.x(), 20.5*in - 1.0*in ,posLayer1.z() + 19.75*in - 1.125*in);
  G4ThreeVector   posRT3_2(posLayer2.x(), -0.5*in,posLayer1.z() + 19.75*in - 2.25*in);

  // Layer 3 positions
  G4ThreeVector   posBrick1_3(posLayer3.x(),0.5*in ,posLayer1.z() + .25*in - 2.5*brick_z );
  G4ThreeVector   posBrick2_3(posLayer3.x(),0.5*in ,posLayer1.z() + .25*in - 1.5*brick_z );
  G4ThreeVector   posBrick3_3(posLayer3.x(),0.5*in ,posLayer1.z() + .25*in - 0.5*brick_z );
  G4ThreeVector   posBrick4_3(posLayer3.x(),0.5*in ,posLayer1.z() + .25*in + 0.5*brick_z );
  G4ThreeVector   posBrick5_3(posLayer3.x(),0.5*in ,posLayer1.z() + .25*in + 1.5*brick_z );
  G4ThreeVector   posBrick6_3(posLayer3.x(),0.5*in ,posLayer1.z() + .25*in + 2.5*brick_z );
  G4ThreeVector   posPL1_3(posLayer3.x(), 0.5*in - 2.5*in - 18.*in,posLayer1.z() + -0.625*in );
  G4ThreeVector   posPL2_3(posLayer3.x(), 0.5*in - 1.5*in - 18.*in,posLayer1.z() + -0.625*in );
  G4ThreeVector   posPR1_3(posLayer3.x(), 0.5*in + 0.5*in + 18.*in,posLayer1.z() + -0.625*in);
  G4ThreeVector   posPR2_3(posLayer3.x(), 0.5*in + 1.5*in + 18.*in,posLayer1.z() + -0.625*in);
  G4ThreeVector   posPB1_3(posLayer3.x(), 0.5*in,posLayer1.z() + -19.75*in + 1.5*in);
  G4ThreeVector   posPB2_3(posLayer3.x(), 0.5*in,posLayer1.z() + -19.75*in + 0.5*in);
  G4ThreeVector   posPT1_3(posLayer3.x(), -0.25*in,posLayer1.z() + 19.75*in - 0.5*in);
  G4ThreeVector   posRL1_3(posLayer3.x(), 0.5*in - 18.0*in - 0.5*in,posLayer1.z() + -0.625*in );
  G4ThreeVector   posRR1_3(posLayer3.x(), 20.5*in - 0.25*in,posLayer1.z() + 19.75*in - 0.5*in);
  G4ThreeVector   posRT1_3(posLayer3.x(), 20.5*in - 1.0*in ,posLayer1.z() + 19.75*in - 1.125*in);
  G4ThreeVector   posRT2_3(posLayer3.x(), -20.5*in + 1.0*in ,posLayer1.z() + 19.75*in - 1.125*in);
  G4ThreeVector   posRT3_3(posLayer3.x(), 0.5*in,posLayer1.z() + 19.75*in - 1.25*in);


  // Lead Castle

  G4double leadCastle_x = 24.25*in;
  G4double leadCastle_y = 24.25*in;
  G4double leadCastle_z = 22.0*in;

  G4Box *leadCastleA = new G4Box("leadCastleA", leadCastle_x/2, leadCastle_y/2, leadCastle_z/2);
  G4Box *leadCastleB = new G4Box("leadCastleB", 6.75*in/2., 16.25*in/2, 14.0*in/2);
  G4Box *leadCastleC = new G4Box("leadCastleC", 8.*in/2., 4.*in/2, 4.*in/2);
  G4Box *leadCastleD = new G4Box("leadCastleD", 9.5*in/2., 16.25*in/2, 14.0*in/2);

  G4ThreeVector castleOffsetB((12.125-(6.75/2.))*in,0,0);
  G4ThreeVector castleOffsetC((12.125-6.75-2.)*in,0,0);
  G4ThreeVector castleOffsetD((-12.125+4.+4.75)*in,0,0);

  G4SubtractionSolid *castle1 = new G4SubtractionSolid("castle1",leadCastleA,leadCastleB,0,castleOffsetB);
  G4SubtractionSolid *castle2 = new G4SubtractionSolid("castle2",castle1,leadCastleC,0,castleOffsetC);
  G4SubtractionSolid *castle = new G4SubtractionSolid("castle",castle2,leadCastleD,0,castleOffsetD);

  G4LogicalVolume * fLogicCastle = new G4LogicalVolume(castle, mats->GetShieldLead(), "fLogicCastle");

  G4ThreeVector   posLayerCastle(water_rad+bioShield_thick +SSPlate_thick  - thermalC_length_big - thermalC_length_small + 78.125*in,-5.5*in,-(worldZ/2.)+floor_thick_out+floor_to_TC_out-3.5*in);


  // Paraffin shield
  G4double parOuterDia = 11.*in;
  G4double parInnerDia = 3.131*in;
  G4double parInnerDepth = 5.*in;
  G4double parOuterDepth = 9.5*in;

  G4VSolid* paraffin1 = new G4Tubs("paraffin1",0.,parOuterDia/2.,parOuterDepth/2.,0,360*deg);
  G4VSolid* paraffin2 = new G4Tubs("paraffin2",0.,parInnerDia/2.,parInnerDepth/2.,0,360*deg);

  G4RotationMatrix* RotParaffin = new G4RotationMatrix;
  RotParaffin->rotateZ(0.);
  RotParaffin->rotateX(0.);
  RotParaffin->rotateY(pi/2.);

  G4ThreeVector paraffinOffset(0,0,-2.25*in);
  G4SubtractionSolid *paraffin = new G4SubtractionSolid("paraffin",paraffin1,paraffin2,0,paraffinOffset);

  G4LogicalVolume * fLogicParaffin = new G4LogicalVolume(paraffin, mats->GetBParaffin(), "fLogicParaffin");

  G4ThreeVector posLayerParaffin(posLayerCastle.x()+castleOffsetD.x(),posLayerCastle.y(),posLayerCastle.z());


  // Craig's Ge detector
  //G4double detHalfZ = 16.891*mm;
  G4double detHalfZ = (4.805/2.)*cm;
  G4AssemblyVolume *fullDet = ConstructBEGe("_insideCastle");
  G4AssemblyVolume *fullDet2 = ConstructBEGe("_outsideCastle");

  G4RotationMatrix* RotDet = new G4RotationMatrix;
  RotDet->rotateZ(-pi/4.);
  RotDet->rotateX(0.);
  RotDet->rotateY(-pi/2.);

  G4ThreeVector posDet = G4ThreeVector(water_rad+bioShield_thick - thermalC_length_big - thermalC_length_small + 66.*in + 8.5*in + 0.5*detHalfZ + 5*mm + 1.1*in,-5.5*in,-(worldZ/2.)+floor_thick_out+floor_to_TC_out-3.5*in);

  G4ThreeVector posDet2 = G4ThreeVector(water_rad+bioShield_thick - thermalC_length_big - thermalC_length_small + 66.*in + 0.5*detHalfZ + 5*mm + 1.1*in,-5.5*in,-(worldZ/2.)+floor_thick_out+floor_to_TC_out+13.5*in);


  // collecting placements here to make it easy to turn stuff off
  fUniverse_phys = new G4PVPlacement(0,G4ThreeVector(),"universe_P",universe_log,0,false,0);
  new G4PVPlacement(0,floorPlacement,floor_log,"floor_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,graphite_Placement,graphite_log,"graphite_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,poolPlacement,pool_log,"pool_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,bioShieldPlacement,bioShield_log,"bioShield_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,poolSSPlacement,poolSS_log,"poolSS_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TC_SS_placement,TC_SS_log,"TC_SS_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TC_Alext_placement,TC_Alext_log,"TC_Alext_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TC_Al_placement,TC_Al_log,"TC_Al_phys",universe_log,false,0,fCheckOverlaps);
  door->MakeImprint(universe_log,posDoor,zeroRot,0,fCheckOverlaps);
  fullDet->MakeImprint(universe_log,posDet,RotDet,0,fCheckOverlaps);
  fullDet2->MakeImprint(universe_log,posDet2,RotDet,0,fCheckOverlaps);

  //shielding->MakeImprint(universe_log,posShield,zeroRot,0,fCheckOverlaps);

  // Lead Wall Place
  new G4PVPlacement(0,posLayerLead,leadWall_log,"leadWall_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posLayerLead2,leadWall2_log,"leadWall2_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posLayerLead3,leadWall3_log,"leadWall3_phys",universe_log,false,0,fCheckOverlaps);

  // Shield Layer 1 Place
  WBsubstruct1->MakeImprint(universe_log,posBrick1_1,zeroRot,0,fCheckOverlaps);
  WBsubstruct2->MakeImprint(universe_log,posBrick2_1,zeroRot,0,fCheckOverlaps);
  WBsubstruct1->MakeImprint(universe_log,posBrick3_1,zeroRot,0,fCheckOverlaps);
  WBsubstruct2->MakeImprint(universe_log,posBrick4_1,zeroRot,0,fCheckOverlaps);
  WBsubstruct1->MakeImprint(universe_log,posBrick5_1,zeroRot,0,fCheckOverlaps);
  WBsubstruct2->MakeImprint(universe_log,posBrick6_1,zeroRot,0,fCheckOverlaps);
  new G4PVPlacement(0,posPL1_1,poly_left_1_layer1,"poly_left_1_layer1_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPL2_1,poly_left_2_layer1,"poly_left_2_layer1_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPR1_1,poly_right_1_layer1,"poly_right_1_layer1_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPR2_1,poly_right_2_layer1,"poly_right_2_layer1_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPT1_1,poly_top_1_layer1,"poly_top_1_layer1_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPT2_1,poly_top_2_layer1,"poly_top_2_layer1_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPT3_1,poly_top_3_layer1,"poly_top_3_layer1_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posRL1_1,rub_left_1_layer1,"rub_left_1_layer1_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posRR1_1,rub_right_1_layer1,"rub_right_1_layer1_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posRT1_1,rub_top_1_layer1,"rub_top_1_layer1_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posRT2_1,rub_top_2_layer1,"rub_top_2_layer1_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posRT3_1,rub_top_3_layer1,"rub_top_3_layer1_phys",universe_log,false,0,fCheckOverlaps);

  // Shield Layer 2 Place
  WBsubstruct2->MakeImprint(universe_log,posBrick1_2,zeroRot,0,fCheckOverlaps);
  WBsubstruct1->MakeImprint(universe_log,posBrick2_2,zeroRot,0,fCheckOverlaps);
  WBsubstruct2->MakeImprint(universe_log,posBrick3_2,zeroRot,0,fCheckOverlaps);
  WBsubstruct1->MakeImprint(universe_log,posBrick4_2,zeroRot,0,fCheckOverlaps);
  WBsubstruct2->MakeImprint(universe_log,posBrick5_2,zeroRot,0,fCheckOverlaps);
  WBsubstruct1->MakeImprint(universe_log,posBrick6_2,zeroRot,0,fCheckOverlaps);
  new G4PVPlacement(0,posPL1_2,poly_left_1_layer2,"poly_left_1_layer2_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPL2_2,poly_left_2_layer2,"poly_left_2_layer2_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPR1_2,poly_right_1_layer2,"poly_right_1_layer2_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPR2_2,poly_right_2_layer2,"poly_right_2_layer2_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPB1_2,poly_bottom_1_layer2,"poly_bottom_1_layer2_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPT1_2,poly_top_1_layer2,"poly_top_1_layer2_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPT2_2,poly_top_2_layer2,"poly_top_2_layer2_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posRL1_2,rub_left_1_layer2,"rub_right_1_layer2_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posRR1_2,rub_right_1_layer2,"rub_left_1_layer2_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posRT1_2,rub_top_1_layer2,"rub_top_1_layer2_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posRT2_2,rub_top_2_layer2,"rub_top_2_layer2_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posRT3_2,rub_top_3_layer2,"rub_top_3_layer2_phys",universe_log,false,0,fCheckOverlaps);


  // Shield Layer 3 Place
  WBsubstruct1->MakeImprint(universe_log,posBrick1_3,zeroRot,0,fCheckOverlaps);
  WBsubstruct2->MakeImprint(universe_log,posBrick2_3,zeroRot,0,fCheckOverlaps);
  WBsubstruct1->MakeImprint(universe_log,posBrick3_3,zeroRot,0,fCheckOverlaps);
  WBsubstruct2->MakeImprint(universe_log,posBrick4_3,zeroRot,0,fCheckOverlaps);
  WBsubstruct1->MakeImprint(universe_log,posBrick5_3,zeroRot,0,fCheckOverlaps);
  WBsubstruct2->MakeImprint(universe_log,posBrick6_3,zeroRot,0,fCheckOverlaps);
  new G4PVPlacement(0,posPL1_3,poly_left_1_layer3,"poly_left_1_layer3_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPL2_3,poly_left_2_layer3,"poly_left_2_layer3_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPR1_3,poly_right_1_layer3,"poly_right_1_layer3_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPR2_3,poly_right_2_layer3,"poly_right_2_layer3_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPB1_3,poly_bottom_1_layer3,"poly_bottom_1_layer3_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPB2_3,poly_bottom_2_layer3,"poly_bottom_2_layer3_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posPT1_3,poly_top_1_layer3,"poly_top_1_layer3_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posRL1_3,rub_left_1_layer3,"rub_left_1_layer3_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posRR1_3,rub_right_1_layer3,"rub_right_1_layer3_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posRT1_3,rub_top_1_layer3,"rub_top_1_layer3_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posRT2_3,rub_top_2_layer3,"rub_top_2_layer3_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posRT3_3,rub_top_3_layer3,"rub_top_3_layer3_phys",universe_log,false,0,fCheckOverlaps);


  // Lead Castle Place
  new G4PVPlacement(0,posLayerCastle,fLogicCastle,"PbCastle_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(RotParaffin,posLayerParaffin,fLogicParaffin,"Paraffin_phys",universe_log,false,0,fCheckOverlaps);


//--------- Visualization attributes -------------------------------
  universe_log->SetVisAttributes(G4VisAttributes::Invisible);


  //aVisAttAlTubes->SetForceWireframe(true);
  aVisAttWater->SetForceWireframe(true);
  aVisAttBioShield->SetForceWireframe(true);
  aVisAttHDConcrete->SetForceWireframe(true);
  aVisAttStainless->SetForceWireframe(true);
  aVisAttSSLiner->SetForceWireframe(true);
  //aVisAttLiqN2->SetForceWireframe(true);

  pool_log->SetVisAttributes(aVisAttWater);
  floor_log->SetVisAttributes(aVisAttHDConcrete);
  bioShield_log->SetVisAttributes(aVisAttBioShield);
  poolSS_log->SetVisAttributes(aVisAttSSLiner);

  graphite_log->SetVisAttributes(aVisAttGraphite);

  TC_SS_log->SetVisAttributes(aVisAttSSLiner);
  TC_Al_log->SetVisAttributes(aVisAttSSLiner);

  door_Liner_log->SetVisAttributes(aVisAttHDConcrete);
  door_Lead_log->SetVisAttributes(aVisAttLead);
  door_Concrete_log->SetVisAttributes(aVisAttHDConcrete);
  door_Steel_log->SetVisAttributes(aVisAttStainless);


  logic_shell->SetVisAttributes(aVisAttWhite);
  logic_shell->SetVisAttributes(aVisAttWater);


  poly_left_1_layer1->SetVisAttributes(aVisAttWhite);
  poly_left_2_layer1->SetVisAttributes(aVisAttWhite);
  poly_right_1_layer1->SetVisAttributes(aVisAttWhite);
  poly_right_2_layer1->SetVisAttributes(aVisAttWhite);
  poly_top_1_layer1->SetVisAttributes(aVisAttWhite);
  poly_top_2_layer1->SetVisAttributes(aVisAttWhite);
  poly_top_3_layer1->SetVisAttributes(aVisAttWhite);

  rub_left_1_layer1->SetVisAttributes(aVisAttGraphite);
  rub_right_1_layer1->SetVisAttributes(aVisAttGraphite);
  rub_top_1_layer1->SetVisAttributes(aVisAttGraphite);
  rub_top_2_layer1->SetVisAttributes(aVisAttGraphite);
  rub_top_3_layer1->SetVisAttributes(aVisAttGraphite);


  poly_left_1_layer2->SetVisAttributes(aVisAttWhite);
  poly_left_2_layer2->SetVisAttributes(aVisAttWhite);
  poly_right_1_layer2->SetVisAttributes(aVisAttWhite);
  poly_right_2_layer2->SetVisAttributes(aVisAttWhite);
  poly_bottom_1_layer2->SetVisAttributes(aVisAttWhite);
  poly_top_1_layer2->SetVisAttributes(aVisAttWhite);
  poly_top_2_layer2->SetVisAttributes(aVisAttWhite);
 
  rub_right_1_layer2->SetVisAttributes(aVisAttGraphite);
  rub_left_1_layer2->SetVisAttributes(aVisAttGraphite);
  rub_top_1_layer2->SetVisAttributes(aVisAttGraphite);
  rub_top_2_layer2->SetVisAttributes(aVisAttGraphite);
  rub_top_3_layer2->SetVisAttributes(aVisAttGraphite);

  poly_left_1_layer3->SetVisAttributes(aVisAttWhite);
  poly_left_2_layer3->SetVisAttributes(aVisAttWhite);
  poly_right_1_layer3->SetVisAttributes(aVisAttWhite);
  poly_right_2_layer3->SetVisAttributes(aVisAttWhite);
  poly_bottom_1_layer3->SetVisAttributes(aVisAttWhite);
  poly_bottom_2_layer3->SetVisAttributes(aVisAttWhite);
  poly_top_1_layer3->SetVisAttributes(aVisAttWhite);

  rub_left_1_layer3->SetVisAttributes(aVisAttGraphite);
  rub_right_1_layer3->SetVisAttributes(aVisAttGraphite);
  rub_top_1_layer3->SetVisAttributes(aVisAttGraphite);
  rub_top_2_layer3->SetVisAttributes(aVisAttGraphite);
  rub_top_3_layer3->SetVisAttributes(aVisAttGraphite);



  return fUniverse_phys;
}


G4VPhysicalVolume *GeometryConstruction::GetWorldVolume() {
   return fUniverse_phys;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
