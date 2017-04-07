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
  SetSensitiveDetector("Craig_log", aCraigSD, true);

  /*
  MinerSD* muVetoTopSD = new MinerSD("/MINERsim/muVetoTop","MS_muVetoTop_hits");
  SDman->AddNewDetector(muVetoTopSD);
  SetSensitiveDetector("muVetoTop_log",muVetoTopSD, true);

  MinerSD* muVetoBottomSD = new MinerSD("/MINERsim/muVetoBottom","MS_muVetoBottom_hits");
  SDman->AddNewDetector(muVetoBottomSD);
  SetSensitiveDetector("muVetoBottom_log",muVetoBottomSD, true);

  MinerSD* backScintSD = new MinerSD("/MINERsim/backScint","MS_backScint_hits");
  SDman->AddNewDetector(backScintSD);
  SetSensitiveDetector("backScint_log",backScintSD, true);

  MinerSD* det1SD = new MinerSD("/MINERsim/det1","MS_det1_hits");
  SDman->AddNewDetector(det1SD);
  SetSensitiveDetector("det1_log",det1SD, true);

  MinerSD* det2SD = new MinerSD("/MINERsim/det2","MS_det2_hits");
  SDman->AddNewDetector(det2SD);
  SetSensitiveDetector("det2_log",det2SD, true);

  MinerSD* det3SD = new MinerSD("/MINERsim/det3","MS_det3_hits");
  SDman->AddNewDetector(det3SD);
  SetSensitiveDetector("det3_log",det3SD, true);

  MinerSD* det4SD = new MinerSD("/MINERsim/det4","MS_det4_hits");
  SDman->AddNewDetector(det4SD);
  SetSensitiveDetector("det4_log",det4SD, true);

  MinerSD* det5SD = new MinerSD("/MINERsim/det5","MS_det5_hits");
  SDman->AddNewDetector(det5SD);
  SetSensitiveDetector("det5_log",det5SD, true);

  MinerSD* det6SD = new MinerSD("/MINERsim/det6","MS_det6_hits");
  SDman->AddNewDetector(det6SD);
  SetSensitiveDetector("det6_log",det6SD, true);

  MinerSD* det7SD = new MinerSD("/MINERsim/det7","MS_det7_hits");
  SDman->AddNewDetector(det7SD);
  SetSensitiveDetector("det7_log",det7SD, true);

  MinerSD* det8SD = new MinerSD("/MINERsim/det8","MS_det8_hits");
  SDman->AddNewDetector(det8SD);
  SetSensitiveDetector("det8_log",det8SD, true);
  */
}


G4VPhysicalVolume* GeometryConstruction::Construct()
{

  // some useful stuff
  G4double in = 2.54*cm;
  G4bool fCheckOverlaps = true;
  const MINERMaterials * mats = MINERMaterials::GetInstance();


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
  G4double thermalC_length_small = (19.625+8.875+25.8125)*in;
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


  G4cout << water_rad+bioShield_thick-thermalC_length_big-thermalC_length_small-SSPlate_thick-gb_depth << 0. << -(worldZ/2.)+floor_thick_out+floor_to_TC_out << G4endl;
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
  G4RotationMatrix *zeroRot = new G4RotationMatrix;
  G4ThreeVector zeroPos;
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
  G4ThreeVector   brick3_place(0.,0.5*(brick_y + wgGap) + 0.5*(brick_x + wgGap),0.);
  G4ThreeVector   brick4_place(0.,-1*(0.5*(brick_y + wgGap) + 0.5*(brick_x + wgGap)),0.);
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
  G4ThreeVector   brick6_place(-0.5*(brick_x + wgGap),-1*(0.5*(brick_y + wgGap)),0.);
  G4ThreeVector   brick7_place(0.5*(brick_x + wgGap),0.5*(brick_y + wgGap),0.);
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
  G4ThreeVector   posLayer3(water_rad+bioShield_thick +SSPlate_thick  - thermalC_length_big - thermalC_length_small + 47.*in,0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  G4ThreeVector   posLayerLead(water_rad+bioShield_thick +SSPlate_thick  - thermalC_length_big - thermalC_length_small + 37.*in,0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  G4ThreeVector   posLayerLead2(water_rad+bioShield_thick +SSPlate_thick  - thermalC_length_big - thermalC_length_small + 49.*in,0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);

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

/*

  // ************************************************** //
  // ************************************************** //
  // ***********    Shielding setup #1  *************** //
  // ************************************************** //
  // ************************************************** //

  // TC small part completely full of poly
  // 8" lead barrier after
  // 2" poly surounding icebox
  // 2" lead after that
  // 6" poly top/bottom, 4" back
  // 2" muon veto top and bottom, 8" nova-like detector back

  // going to treat this as a single object 
  // because it's easier to define locally then place 
  // it in the thermal column after
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


  // reactor poly shielding, completely filling the small part for now
  G4Box *poly_TC = new G4Box("poly_TC", TCPoly_thick/2.,thermalC_width_small/2.,thermalC_height_small/2.);
  G4LogicalVolume *polyTC_log = new G4LogicalVolume(poly_TC,mats->GetBoratedPoly05(),"polyTC_log");

  // reactor Front lead shielding
  G4Box *lead_TC = new G4Box("lead_TC",TCLead_thick/2.,thermalC_height_big/2.,thermalC_height_big/2.);
  G4LogicalVolume *leadTC_log = new G4LogicalVolume(lead_TC,mats->GetShieldLead(),"leadTC_log");


  // icebox oversimplified for now, three concentric cans with a can thickness of space between each
  G4Tubs *icebox1_out = new G4Tubs("icebox1_out", 0., iceboxRadius, iceboxHeight/2.,0,360*deg);
  G4Tubs *icebox1_in = new G4Tubs("icebox1_in", 0., iceboxRadius - iceboxThick, iceboxHeight/2. - iceboxThick,0,360*deg);
  G4SubtractionSolid *icebox1 = new G4SubtractionSolid("icebox1",icebox1_out,icebox1_in);
  G4LogicalVolume *icebox1_log = new G4LogicalVolume(icebox1,mats->GetCopper(),"icebox1_log");

  G4Tubs *icebox2_out = new G4Tubs("icebox2_out", 0., iceboxRadius - 2*iceboxThick, iceboxHeight/2. - 2*iceboxThick,0,360*deg);
  G4Tubs *icebox2_in = new G4Tubs("icebox2_in", 0., iceboxRadius - 3*iceboxThick, iceboxHeight/2. - 3*iceboxThick,0,360*deg);
  G4SubtractionSolid *icebox2 = new G4SubtractionSolid("icebox2",icebox2_out,icebox2_in);
  G4LogicalVolume *icebox2_log = new G4LogicalVolume(icebox2,mats->GetCopper(),"icebox2_log");

  G4Tubs *icebox3_out = new G4Tubs("icebox3_out", 0., iceboxRadius - 4*iceboxThick, iceboxHeight/2. - 4*iceboxThick,0,360*deg);
  G4Tubs *icebox3_in = new G4Tubs("icebox3_in", 0., iceboxRadius - 5*iceboxThick, iceboxHeight/2. - 5*iceboxThick,0,360*deg);
  G4SubtractionSolid *icebox3 = new G4SubtractionSolid("icebox3",icebox3_out,icebox3_in);
  G4LogicalVolume *icebox3_log = new G4LogicalVolume(icebox3,mats->GetCopper(),"icebox3_log");

  // Let's put some detectors in the icebox, they can just float for now
  G4Tubs *det = new G4Tubs("det",0.,2*in,detThick/2.,0,360*deg);
  G4LogicalVolume *det1_log = new G4LogicalVolume(det, mats->GetDetGe(),"det1_log");
  G4LogicalVolume *det2_log = new G4LogicalVolume(det, mats->GetDetGe(),"det2_log");
  G4LogicalVolume *det3_log = new G4LogicalVolume(det, mats->GetDetGe(),"det3_log");
  G4LogicalVolume *det4_log = new G4LogicalVolume(det, mats->GetDetGe(),"det4_log");
  G4LogicalVolume *det5_log = new G4LogicalVolume(det, mats->GetDetSi(),"det5_log");
  G4LogicalVolume *det6_log = new G4LogicalVolume(det, mats->GetDetSi(),"det6_log");
  G4LogicalVolume *det7_log = new G4LogicalVolume(det, mats->GetDetSi(),"det7_log");
  G4LogicalVolume *det8_log = new G4LogicalVolume(det, mats->GetDetSi(),"det8_log");
  // should have slightly more than 3" for spacing of detectors

  // poly around the icebox
  G4Box *poly_IBshield_out = new G4Box("poly_IBshield_out",iceboxRadius+iceboxPoly_thick,iceboxRadius+iceboxPoly_thick,iceboxHeight/2. + iceboxPoly_thick);
  G4Box *poly_IBshield_in = new G4Box("poly_IBshield_in",iceboxRadius,iceboxRadius,iceboxHeight/2.);
  G4SubtractionSolid *poly_IBshield = new G4SubtractionSolid("poly_IBshield",poly_IBshield_out,poly_IBshield_in);
  G4LogicalVolume *polyIBshield_log = new G4LogicalVolume(poly_IBshield,mats->GetBoratedPoly05(),"poly_IBshield_log");

  // lead around the icebox
  G4Box *lead_IBshield_out = new G4Box("lead_IBshield_out",iceboxRadius+iceboxPoly_thick+iceboxLead_thick,iceboxRadius+iceboxPoly_thick+iceboxLead_thick,iceboxHeight/2. +iceboxPoly_thick+ iceboxLead_thick);
  G4Box *lead_IBshield_in = new G4Box("lead_IBshield_in",iceboxRadius+iceboxPoly_thick,iceboxRadius+iceboxPoly_thick,iceboxHeight/2. + iceboxPoly_thick);
  G4SubtractionSolid *lead_IBshield = new G4SubtractionSolid("lead_IBshield",lead_IBshield_out,lead_IBshield_in);
  G4LogicalVolume *leadIBshield_log = new G4LogicalVolume(lead_IBshield,mats->GetShieldLead(),"lead_IBshield_log");

  // poly around the IB lead shielding
  G4Box *poly_outershield_out = new G4Box("poly_outershield_out",iceboxRadius+iceboxPoly_thick+iceboxLead_thick,iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick,iceboxHeight/2. + iceboxPoly_thick+iceboxLead_thick+outerPoly_thick);
  //current design has only 4" on doorside and no extra on backside
  G4Box *poly_outershield_in = new G4Box("poly_outershield_in",iceboxRadius+iceboxPoly_thick+iceboxLead_thick,iceboxRadius+iceboxPoly_thick+iceboxLead_thick,iceboxHeight/2.+iceboxPoly_thick+iceboxLead_thick);
  G4SubtractionSolid *poly_outershield_shell = new G4SubtractionSolid("poly_outershield_shell",poly_outershield_out,poly_outershield_in);
  G4ThreeVector  poly_outer_back_place(iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_back_thick/2.,0,0);  
  G4Box *poly_outershield_back = new G4Box("poly_outershield_back",outerPoly_back_thick/2.,iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick,iceboxHeight/2.+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick);
  G4UnionSolid *poly_outershield = new G4UnionSolid("poly_outershield",poly_outershield_shell,poly_outershield_back,0,poly_outer_back_place);
  G4LogicalVolume *polyoutershield_log = new G4LogicalVolume(poly_outershield,mats->GetBoratedPoly05(),"poly_outershield_log");

  //scintillator veto, just top and bottom, this means there will be some extra space on the sides
  G4Box *muVeto = new G4Box("muVeto",(thermalC_length_big-aluminumPlate_thick-TCLead_thick)/2.,iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick,muVeto_thick/2.);
  G4LogicalVolume *muVetoTop_log = new G4LogicalVolume(muVeto,mats->GetAcrylic() ,"muVetoTop_log");
  G4LogicalVolume *muVetoBottom_log = new G4LogicalVolume(muVeto,mats->GetAcrylic() ,"muVetoBottom_log");
  //back scintillator
  G4Box *backScint = new G4Box("backScint",backScint_thick/2.,iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick,iceboxHeight/2.+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick);
  G4LogicalVolume *backScint_log = new G4LogicalVolume(backScint,mats->GetAcrylic() ,"backScint_log");


  // Lead next to the liner on outside of thermal column
  G4Box *lead_outershield_in = new G4Box("lead_outershield_in",2*m,thermalC_height_big/2. -outerLead_thick,thermalC_height_big/2. -outerLead_thick);
  G4Box *lead_outershield_out = new G4Box("lead_outershield_out",(thermalC_length_big-aluminumPlate_thick-TCLead_thick)/2.,thermalC_height_big/2.,thermalC_height_big/2.);
  G4SubtractionSolid *lead_outershield = new G4SubtractionSolid("lead_outershield",lead_outershield_out,lead_outershield_in);
  G4LogicalVolume *leadoutershield_log = new G4LogicalVolume(lead_outershield,mats->GetShieldLead(),"lead_outershield_log");

  // ok, let's connect everything into an assembly
  G4AssemblyVolume *shielding = new G4AssemblyVolume();
  // placements for everything relative to (0,0,0) of the assembly
  G4ThreeVector   posShield(water_rad+bioShield_thick +SSPlate_thick  - thermalC_length_big + iceboxRadius+iceboxPoly_thick+iceboxLead_thick+TCLead_thick , 0, -(worldZ/2.)+floor_thick_out+floor_to_TC_out);

  G4ThreeVector   polyTC_place(-1*(TCPoly_thick/2.  + SSPlate_thick+ TCLead_thick + iceboxRadius+iceboxPoly_thick+iceboxLead_thick) ,0.,0.);
  G4ThreeVector   leadTC_place(-1*(0.5*TCLead_thick +  iceboxRadius+iceboxPoly_thick+iceboxLead_thick),0.,0.);
  G4ThreeVector   det1_place(0.,0.,(detThick/2. + detSpace/2. +  3.*detThick + 3.*detSpace));
  G4ThreeVector   det2_place(0.,0.,(detThick/2. + detSpace/2. +  2.*detThick + 2.*detSpace));
  G4ThreeVector   det3_place(0.,0.,(detThick/2. + detSpace/2. +  1.*detThick + 1.*detSpace));
  G4ThreeVector   det4_place(0.,0.,(detThick/2. + detSpace/2. +  0.*detThick + 0.*detSpace));
  G4ThreeVector   det5_place(0.,0.,-1 * (detThick/2. + detSpace/2. +  0.*detThick + 0.*detSpace));
  G4ThreeVector   det6_place(0.,0.,-1 * (detThick/2. + detSpace/2. +  1.*detThick + 1.*detSpace));
  G4ThreeVector   det7_place(0.,0.,-1 * (detThick/2. + detSpace/2. +  2.*detThick + 2.*detSpace));
  G4ThreeVector   det8_place(0.,0.,-1 * (detThick/2. + detSpace/2. +  3.*detThick + 3.*detSpace));
  G4ThreeVector   can1_place(0.,0.,0.);
  G4ThreeVector   can2_place(0.,0.,0.);
  G4ThreeVector   can3_place(0.,0.,0.);
  G4ThreeVector   polyIB_place(0.,0.,0.);
  G4ThreeVector   leadIB_place(0.,0.,0.);
  G4ThreeVector   polyOut_place(0.,0.,0.);
  G4ThreeVector   leadOut_place((thermalC_length_big-aluminumPlate_thick-TCLead_thick)/2. - (iceboxRadius+iceboxPoly_thick+iceboxLead_thick ),0.,0.); 
  G4ThreeVector   muVetoTop_place(leadOut_place.x(),0,(iceboxHeight/2. +iceboxPoly_thick+ iceboxLead_thick + outerPoly_thick+muVeto_thick/2.));
  G4ThreeVector   muVetoBot_place(leadOut_place.x(),0,-1*(iceboxHeight/2. +iceboxPoly_thick+ iceboxLead_thick + outerPoly_thick+muVeto_thick/2.));
  G4ThreeVector   backScint_place(iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_back_thick+backScint_thick/2.,0,0);
  // Detector assemble!
  shielding->AddPlacedVolume(icebox1_log,can1_place,zeroRot);
  shielding->AddPlacedVolume(icebox2_log,can2_place,zeroRot);
  shielding->AddPlacedVolume(icebox3_log,can3_place,zeroRot);
  shielding->AddPlacedVolume(polyTC_log,polyTC_place,zeroRot);
  shielding->AddPlacedVolume(leadTC_log,leadTC_place,zeroRot);
  shielding->AddPlacedVolume(det1_log,det1_place,zeroRot);
  shielding->AddPlacedVolume(det2_log,det2_place,zeroRot);
  shielding->AddPlacedVolume(det3_log,det3_place,zeroRot);
  shielding->AddPlacedVolume(det4_log,det4_place,zeroRot);
  shielding->AddPlacedVolume(det5_log,det5_place,zeroRot);
  shielding->AddPlacedVolume(det6_log,det6_place,zeroRot);
  shielding->AddPlacedVolume(det7_log,det7_place,zeroRot);
  shielding->AddPlacedVolume(det8_log,det8_place,zeroRot);
  shielding->AddPlacedVolume(polyIBshield_log,polyIB_place,zeroRot);
  shielding->AddPlacedVolume(leadIBshield_log,leadIB_place,zeroRot);
  shielding->AddPlacedVolume(polyoutershield_log,polyOut_place,zeroRot);
  shielding->AddPlacedVolume(leadoutershield_log,leadOut_place,zeroRot);
  shielding->AddPlacedVolume(muVetoTop_log,muVetoTop_place,zeroRot);
  shielding->AddPlacedVolume(muVetoBottom_log,muVetoBot_place,zeroRot);
  shielding->AddPlacedVolume(backScint_log,backScint_place,zeroRot);

*/

  // Craig's Ge detector
  // Need to put N2 into the detector still

  G4double detRad = 50.8*mm; // 4" diameter
  G4double detHalfZ = 16.891*mm; // 1.33" thick 

  G4VSolid* target = new G4Tubs("Target",0.,2.5*cm,(4.805/2.)*cm,0,360*deg);
  G4RotationMatrix* RotDet = new G4RotationMatrix;
  RotDet->rotateZ(-pi/4.);
  RotDet->rotateX(0.);
  RotDet->rotateY(-pi/2.);
  G4LogicalVolume *fLogicDet = new G4LogicalVolume(target, mats->GetDetGe(),"Craig_log");

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

  G4Tubs *outTubeA = new G4Tubs("outTubeA", 0., tubeAdia/2., tubeAlen/2.,0,360*deg);
  G4Tubs *inTubeA = new G4Tubs("inTubeA", 0., tubeAdia/2. - tubeAthick, tubeAlen/2. - tubeAthick/2.,0,360*deg);
  // right now the above defines the aluminum tube surrounding the detector, and has a 1/16" aluminum window, probably too thick?
  G4ThreeVector inTubeAOffset(0,0,-tubeAthick/2.);
  G4SubtractionSolid *detTubeA = new G4SubtractionSolid("detTubeA",outTubeA,inTubeA,0,inTubeAOffset);

  G4Tubs *outTubeB = new G4Tubs("outTubeB", 0., tubeBdia/2., tubeBlen/2.,0,360*deg);
  G4Tubs *inTubeB = new G4Tubs("inTubeB", 0., tubeBdia/2. - tubeBthick, tubeBlen/2. - tubeBthick/2.,0,360*deg);
  G4ThreeVector inTubeBOffset(0,0,+tubeBthick/2.);
  G4SubtractionSolid *detTubeB = new G4SubtractionSolid("detTubeB",outTubeB,inTubeB,0,inTubeBOffset);

  G4Tubs *detTubeC = new G4Tubs("detTubeC", tubeCdia/2. - tubeCthick, tubeCdia/2., tubeClen/2.,0,360*deg);
  G4Tubs *inTubeC = new G4Tubs("inTubeC", 0., tubeCdia/2. - tubeCthick, tubeClen/2.,0,360*deg);

  G4Tubs *outTubeD = new G4Tubs("outTubeD", 0., tubeDdia/2., tubeDlen/2.,0,360*deg);
  G4Tubs *inTubeD = new G4Tubs("inTubeD", 0., tubeDdia/2. - tubeDthick, tubeDlen/2. - tubeDthick,0,360*deg);
  G4SubtractionSolid *detTubeD = new G4SubtractionSolid("detTubeD",outTubeD,inTubeD);


  G4AssemblyVolume *fullDet = new G4AssemblyVolume();
  G4ThreeVector tubeAPos(0,0,                 -1*(tubeAlen/2. - tubeAthick - (4.805/2.)*cm) + (1./16.)*in); //detector 1/16" from front aluminum
  G4ThreeVector tubeBPos(0,0,                 (tubeAPos.z() - tubeAlen/2. - tubeBlen/2.)); 
  G4ThreeVector tubeCPos(0,0 - (tubeBdia/3.), (tubeAPos.z() - tubeAlen/2. - tubeBlen - tubeClen/2.));
  G4ThreeVector tubeDPos(0,0,                 (tubeAPos.z() - tubeAlen/2. - tubeBlen - tubeClen - tubeDlen/2.)); 

  G4ThreeVector tubeA_N2_Pos = tubeAPos+inTubeAOffset;
  G4ThreeVector tubeB_N2_Pos = tubeBPos+inTubeBOffset;

  G4Tubs *liqN2_inTubeA = new G4Tubs("liqN2_inTubeA", 0., tubeAdia/2. - tubeAthick, tubeAlen/2. - tubeAthick/2.,0,360*deg);
  G4SubtractionSolid *liqN2_TubeA = new G4SubtractionSolid("liqN2_TubeA",liqN2_inTubeA,target,0,-1*tubeA_N2_Pos);
  G4Tubs *liqN2_TubeB = new G4Tubs("liqN2_TubeB", 0., tubeBdia/2. - tubeBthick, tubeBlen/2. - tubeBthick/2.,0,360*deg);
  G4Tubs *liqN2_TubeC = new G4Tubs("liqN2_TubeC", 0., tubeCdia/2. - tubeCthick, tubeClen/2.,0,360*deg);
  G4Tubs *liqN2_TubeD = new G4Tubs("liqN2_TubeD", 0., tubeDdia/2. - tubeDthick, tubeDlen/2. - tubeDthick,0,360*deg);


  // canberra lead shield
  G4Tubs *InnerCanberra = new G4Tubs("InnerCanberra", 0., (279.*mm)/2.,203.*mm,0.,360*deg);
  G4Tubs *CopperCanberraI = new G4Tubs("CopperCanberraI", 0., (279.*mm)/2. + 1.6*mm,203.*mm + 1.6*mm,0.,360*deg);
  G4Tubs *TinCanberraI = new G4Tubs("TinCanberraI", 0., (279.*mm)/2. + 1.6*mm + 1.0*mm,203.*mm + 1.6*mm + 1.0*mm,0.,360*deg);
  G4Tubs *LeadCanberraI = new G4Tubs("LeadCanberraI", 0., (279.*mm)/2. + 1.6*mm + 1.0*mm + 4.*in,203.*mm + 1.6*mm + 1.0*mm + 4.*in,0.,360*deg);
  G4Tubs *SteelCanberraI = new G4Tubs("SteelCanberraI", 0., (279.*mm)/2. + 1.6*mm + 1.0*mm + 4.*in + 9.5*mm,203.*mm + 1.6*mm + 1.0*mm + 4.*in + 9.5*mm,0.,360*deg);

  G4SubtractionSolid *CopperCanberraI2 = new G4SubtractionSolid("CopperCanberraI2",CopperCanberraI,InnerCanberra);
  G4SubtractionSolid *TinCanberraI2 = new G4SubtractionSolid("TinCanberraI2",TinCanberraI,CopperCanberraI);
  G4SubtractionSolid *LeadCanberraI2 = new G4SubtractionSolid("LeadCanberraI2",LeadCanberraI,TinCanberraI);
  G4SubtractionSolid *SteelCanberraI2 = new G4SubtractionSolid("SteelCanberraI2",SteelCanberraI,LeadCanberraI);

  G4Tubs *bottomHoleCanberra = new G4Tubs("bottomHoleCanberra",0.,(3.+(2./8.))*in/2.,(1.6*mm + 1.0*mm + 4.*in + 9.5*mm)/2.,0.,360*deg);
  G4ThreeVector bHC_Pos(0,0,-1*(1.6*mm + 1.0*mm + 4.*in + 9.5*mm)/2. - 203.*mm);

  G4SubtractionSolid *CopperCanberra = new G4SubtractionSolid("CopperCanberra",CopperCanberraI2,bottomHoleCanberra,0,bHC_Pos);
  G4SubtractionSolid *TinCanberra = new G4SubtractionSolid("TinCanberra",TinCanberraI2,bottomHoleCanberra,0,bHC_Pos);
  G4SubtractionSolid *LeadCanberra = new G4SubtractionSolid("LeadCanberra",LeadCanberraI2,bottomHoleCanberra,0,bHC_Pos);
  G4SubtractionSolid *SteelCanberra = new G4SubtractionSolid("SteelCanberra",SteelCanberraI2,bottomHoleCanberra,0,bHC_Pos);

  G4ThreeVector ShieldCanberra_Pos(0,0,645*mm/2. + tubeDPos.z() + tubeDlen/2.);


  G4LogicalVolume *fLogicDetTubeA = new G4LogicalVolume(detTubeA, mats->GetAluminum(),"detTubeA_log");
  G4LogicalVolume *fLogicDetTubeB = new G4LogicalVolume(detTubeB, mats->GetAluminum(),"detTubeB_log");
  G4LogicalVolume *fLogicDetTubeC = new G4LogicalVolume(detTubeC, mats->GetStainlessSteel(),"detTubeC_log");
  G4LogicalVolume *fLogicDetTubeD = new G4LogicalVolume(detTubeD, mats->GetStainlessSteel(),"detTubeD_log");
  G4LogicalVolume *fLogicDetLiN2A = new G4LogicalVolume(liqN2_TubeA,mats->GetLiN2(),"liqN2_TubeA_log");
  G4LogicalVolume *fLogicDetLiN2B = new G4LogicalVolume(liqN2_TubeB,mats->GetLiN2(),"liqN2_TubeB_log");
  G4LogicalVolume *fLogicDetLiN2C = new G4LogicalVolume(liqN2_TubeC,mats->GetLiN2(),"liqN2_TubeC_log");
  G4LogicalVolume *fLogicDetLiN2D = new G4LogicalVolume(liqN2_TubeD,mats->GetLiN2(),"liqN2_TubeD_log");
  //G4LogicalVolume *fLogicDetLead = new G4LogicalVolume(detLeadShield,mats->GetShieldLead(),"detLeadShield_log");
  G4LogicalVolume *fLogicDetCSCu = new G4LogicalVolume(CopperCanberra,mats->GetCopper(),"CopperCanberra_log");
//  G4LogicalVolume *fLogicDetCSTin = new G4LogicalVolume(TinCanberra,mats->GetTin(),"TinCanberra_log");
  G4LogicalVolume *fLogicDetCSTin = new G4LogicalVolume(TinCanberra,mats->GetCadmium(),"TinCanberra_log");  
  G4LogicalVolume *fLogicDetCSPb = new G4LogicalVolume(LeadCanberra,mats->GetShieldLead(),"LeadCanberra_log");
  G4LogicalVolume *fLogicDetCSSteel = new G4LogicalVolume(SteelCanberra,mats->GetStainlessSteel(),"SteelCanberra_log");

  fullDet->AddPlacedVolume(fLogicDetTubeA,tubeAPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetTubeB,tubeBPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetTubeC,tubeCPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetTubeD,tubeDPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetLiN2A,tubeA_N2_Pos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetLiN2B,tubeB_N2_Pos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetLiN2C,tubeCPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetLiN2D,tubeDPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDet,zeroPos,zeroRot);
  //fullDet->AddPlacedVolume(fLogicDetLead,detLead_Pos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetCSCu,ShieldCanberra_Pos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetCSTin,ShieldCanberra_Pos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetCSPb,ShieldCanberra_Pos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetCSSteel,ShieldCanberra_Pos,zeroRot);

//  G4ThreeVector posDet = G4ThreeVector(water_rad+bioShield_thick - 7.5*in - tubeDlen - tubeClen - tubeBlen - tubeAlen + 0.5*detHalfZ,0.,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);  // 3 ft outside thermal column
  G4ThreeVector posDet = G4ThreeVector(water_rad+bioShield_thick - 6.*in - tubeDlen - tubeClen - tubeBlen - tubeAlen + 0.5*detHalfZ + (1./16.)*in,0.,-(worldZ/2.)+floor_thick_out+floor_to_TC_out-thermalC_height_big/2.+12.*in);


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
  //shielding->MakeImprint(universe_log,posShield,zeroRot,0,fCheckOverlaps);

  // Lead Wall Place
  new G4PVPlacement(0,posLayerLead,leadWall_log,"leadWall_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,posLayerLead2,leadWall2_log,"leadWall2_phys",universe_log,false,0,fCheckOverlaps);

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


//--------- Visualization attributes -------------------------------
  universe_log->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* aVisAttWater = new G4VisAttributes(G4Colour(0.0,0,1.0));
  G4VisAttributes* aVisAttGraphite = new G4VisAttributes(G4Colour(.2,0.2,.2));
  G4VisAttributes* aVisAttLead = new G4VisAttributes(G4Colour(.1,0.1,.1));
  G4VisAttributes* aVisAttPoly = new G4VisAttributes(G4Colour(31./255.,191./255.,127./255.));
  G4VisAttributes* aVisAttAlTubes = new G4VisAttributes(G4Colour(.5,.5,.5));
  G4VisAttributes* aVisAttStainless = new G4VisAttributes(G4Colour(.5,.5,.5));
  G4VisAttributes* aVisAttHDConcrete = new G4VisAttributes(G4Colour(0.3,0.3,.3));
  G4VisAttributes* aVisAttBioShield = new G4VisAttributes(G4Colour(0.55,0.55,.55));
  G4VisAttributes* aVisAttSSLiner = new G4VisAttributes(G4Colour(0.4,0.4,.4));
  G4VisAttributes* aVisAttLiqN2 = new G4VisAttributes(G4Colour(0.0,0.0,.8));
  G4VisAttributes* aVisAttDet = new G4VisAttributes(G4Colour(145./255.,151./255.,171./255.));
  G4VisAttributes* aVisAttCu = new G4VisAttributes(G4Colour(252./255.,187./255.,34./255.));
  G4VisAttributes* aVisAttScint = new G4VisAttributes(G4Colour(245./255.,241./255.,12./255.));
  G4VisAttributes* aVisAttWhite = new G4VisAttributes(G4Colour::White());

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

  fLogicDetTubeA->SetVisAttributes(aVisAttAlTubes);
  fLogicDetTubeB->SetVisAttributes(aVisAttAlTubes);
  fLogicDetTubeC->SetVisAttributes(aVisAttStainless);
  fLogicDetTubeD->SetVisAttributes(aVisAttStainless);
  fLogicDetLiN2A->SetVisAttributes(aVisAttLiqN2);
  fLogicDetLiN2B->SetVisAttributes(aVisAttLiqN2);
  fLogicDetLiN2C->SetVisAttributes(aVisAttLiqN2);
  fLogicDetLiN2D->SetVisAttributes(aVisAttLiqN2);
  fLogicDet->SetVisAttributes(aVisAttWater);

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


  //icebox1_log->SetVisAttributes(aVisAttCu);
  //icebox2_log->SetVisAttributes(aVisAttCu);
  //icebox3_log->SetVisAttributes(aVisAttCu);
  //polyTC_log->SetVisAttributes(aVisAttPoly);
  //leadTC_log->SetVisAttributes(aVisAttLead);
  //det1_log->SetVisAttributes(aVisAttDet);
  //det2_log->SetVisAttributes(aVisAttDet);
  //det3_log->SetVisAttributes(aVisAttDet);
  //det4_log->SetVisAttributes(aVisAttDet);
  //det5_log->SetVisAttributes(aVisAttDet);
  //det6_log->SetVisAttributes(aVisAttDet);
  //det7_log->SetVisAttributes(aVisAttDet);
  //det8_log->SetVisAttributes(aVisAttDet);
  //polyIBshield_log->SetVisAttributes(aVisAttPoly);
  //leadIBshield_log->SetVisAttributes(aVisAttLead);
  //polyoutershield_log->SetVisAttributes(aVisAttPoly);
  //leadoutershield_log->SetVisAttributes(aVisAttLead);
  //muVetoTop_log->SetVisAttributes(aVisAttScint);
  //muVetoBottom_log->SetVisAttributes(aVisAttScint);
  //backScint_log->SetVisAttributes(aVisAttScint);

  fLogicDetCSCu->SetVisAttributes(aVisAttCu);
  fLogicDetCSTin->SetVisAttributes(aVisAttAlTubes);
  fLogicDetCSPb->SetVisAttributes(aVisAttLead);
  fLogicDetCSSteel->SetVisAttributes(aVisAttStainless);

  return fUniverse_phys;
}


G4VPhysicalVolume *GeometryConstruction::GetWorldVolume() {
   return fUniverse_phys;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
