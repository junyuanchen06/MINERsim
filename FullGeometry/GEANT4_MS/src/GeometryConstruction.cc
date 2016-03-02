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

}


G4VPhysicalVolume* GeometryConstruction::Construct()
{

  // some useful stuff
  G4NistManager* nistManager = G4NistManager::Instance();
  G4double in = 2.54*cm;
  G4bool fCheckOverlaps = true;



  // material definitions go here:

  G4Material* vac = new G4Material("Vacuum", 1., 1.01*g/mole, universe_mean_density,
                                 kStateGas,0.000017*kelvin,1.e-19*pascal);


  //Define Stainless Steel
  G4Element* C  = nistManager->FindOrBuildElement("C");
  G4Element* Si = nistManager->FindOrBuildElement("Si");
  G4Element* Cr = nistManager->FindOrBuildElement("Cr");
  G4Element* Mn = nistManager->FindOrBuildElement("Mn");
  G4Element* Fe = nistManager->FindOrBuildElement("Fe");
  G4Element* Ni = nistManager->FindOrBuildElement("Ni");
  G4Material* StainlessSteel = new G4Material("StainlessSteel", 8.06*g/cm3, 6);
  StainlessSteel->AddElement(C, 0.001);
  StainlessSteel->AddElement(Si, 0.007);
  StainlessSteel->AddElement(Cr, 0.18);
  StainlessSteel->AddElement(Mn, 0.01);
  StainlessSteel->AddElement(Fe, 0.712);
  StainlessSteel->AddElement(Ni, 0.09);

  G4Material* liquidNitrogen = new G4Material("LiquidNitrogen",7.,14.01*g/mole,0.808*g/cm3, kStateLiquid,77.*kelvin);
  G4Material* air  = nistManager->FindOrBuildMaterial("G4_AIR");


  // this is temporary until we have better understanding of HD concrete composition
  //G4double concreteDensity = 3.53*g/cm3;
  //G4Material* HDConcrete = nistManager->BuildMaterialWithNewDensity("HDConcrete","G4_CONCRETE",concreteDensity);

  //G4_POLYETHYLENE  (C_2H_4)_N-Polyethylene  0.94 g/cm3,  
  //             1     0.143711
  //             6     0.856289
  G4Material *boron   = nistManager->FindOrBuildMaterial("G4_B");
  G4Material *copper   = nistManager->FindOrBuildMaterial("G4_Cu");

  G4Material *det_Ge = nistManager->FindOrBuildMaterial("G4_Ge");
  G4Material *det_Si = nistManager->FindOrBuildMaterial("G4_Si");

  G4Material *shieldLead = nistManager->FindOrBuildMaterial("G4_Pb");     
  G4Material *HDPE = nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");  
  //http://www.shieldwerx.com
  // 5% borated polyethylene = SWX203
  G4Material *boratedpoly05 = new G4Material("BoratedPoly05",1.06*g/cm3, 2); 
  boratedpoly05->AddMaterial(boron,0.05);
  boratedpoly05->AddMaterial(HDPE,0.95);

  // 30% borated polyethylene = SWX210
  G4Material *boratedpoly30 = new G4Material("BoratedPoly30",1.19*g/cm3, 2);
  boratedpoly30->AddMaterial(boron,0.3);
  boratedpoly30->AddMaterial(HDPE,0.70);


  // stand in for boraflex
  G4Material *Boraflex = new G4Material("Boraflex",1.64*g/cm3, 2);
  Boraflex->AddMaterial(boron,0.276);
  Boraflex->AddMaterial(HDPE,0.724);

  // HD Concrete in NSC
  G4Element *H1 = new G4Element("H1","H1", 1, 1.0078*g/mole);
  G4Element *O16 = new G4Element("O16","O16", 8, 15.9949*g/mole);
  G4Element *Mg24 = new G4Element("Mg24","MG24", 12, 23.985*g/mole);
  G4Element *Mg25 = new G4Element("Mg25","Mg25", 12, 24.9858*g/mole);
  G4Element *Mg26 = new G4Element("Mg26","Mg26", 12, 25.9826*g/mole);
  G4Element *Al27 = new G4Element("Al27","Al27", 13, 26.9815*g/mole);
  G4Element *Si28 = new G4Element("Si28","Si28", 14, 27.9769*g/mole);
  G4Element *Si29 = new G4Element("Si29","Si29", 14, 28.9765*g/mole);
  G4Element *Si30 = new G4Element("Si30","Si30", 14, 29.9738*g/mole);
  G4Element *S32 = new G4Element("S32","S32", 16, 31.9721*g/mole);
  G4Element *S33 = new G4Element("S33","S33", 16, 32.9745*g/mole);
  G4Element *S34 = new G4Element("S34","S34", 16, 33.9679*g/mole);
  G4Element *Ca40 = new G4Element("Ca40","Ca40", 20, 39.9626*g/mole);
  G4Element *Ca42 = new G4Element("Ca42","Ca42", 20, 41.9586*g/mole);
  G4Element *Ca43 = new G4Element("Ca43","Ca43", 20, 42.9588*g/mole);
  G4Element *Ca44 = new G4Element("Ca44","Ca44", 20, 43.9555*g/mole);
  G4Element *Ca46 = new G4Element("Ca46","Ca46", 20, 45.9537*g/mole);
  G4Element *Ca48 = new G4Element("Ca48","Ca48", 20, 47.9525*g/mole);
  G4Element *Ti46 = new G4Element("Ti46","Ti46", 22, 45.9526*g/mole);
  G4Element *Ti47 = new G4Element("Ti47","Ti47", 22, 46.9518*g/mole);
  G4Element *Ti48 = new G4Element("Ti48","Ti48", 22, 47.9479*g/mole);
  G4Element *Ti49 = new G4Element("Ti49","Ti49", 22, 48.9479*g/mole);
  G4Element *Ti50 = new G4Element("Ti50","Ti50", 22, 49.9448*g/mole);
  G4Element *V51 = new G4Element("V51","V51", 23, 50.9440*g/mole);
  G4Element *Cr50 = new G4Element("Cr50","Cr50", 24, 49.9460*g/mole);
  G4Element *Cr52 = new G4Element("Cr52","Cr52", 24, 51.9405*g/mole);
  G4Element *Cr53 = new G4Element("Cr53","Cr53", 24, 52.9407*g/mole);
  G4Element *Cr54 = new G4Element("Cr54","Cr54", 24, 53.9389*g/mole);
  G4Element *Mn55 = new G4Element("Mn55","Mn55", 25, 54.9380*g/mole);
  G4Element *Fe54 = new G4Element("Fe54","Fe54", 26, 53.9396*g/mole);
  G4Element *Fe56 = new G4Element("Fe56","Fe56", 26, 55.9308*g/mole);
  G4Element *Fe57 = new G4Element("Fe57","Fe57", 26, 56.9354*g/mole);
  G4Element *Fe58 = new G4Element("Fe58","Fe58", 26, 57.9333*g/mole);


  G4Material *HDConcrete = new G4Material("HDConcrete",3.53*g/cm3, 33);
  HDConcrete->AddElement(H1,.003113);
  HDConcrete->AddElement(O16,.330504);
  HDConcrete->AddElement(Mg24,.007279);
  HDConcrete->AddElement(Mg25,.00096);
  HDConcrete->AddElement(Mg26,.001099);
  HDConcrete->AddElement(Al27,.023468);
  HDConcrete->AddElement(Si28,.023657);
  HDConcrete->AddElement(Si29,.001241);
  HDConcrete->AddElement(Si30,.000852);
  HDConcrete->AddElement(S32,.001341);
  HDConcrete->AddElement(S33,.000011);
  HDConcrete->AddElement(S34,.000063);
  HDConcrete->AddElement(Ca40,.068668);
  HDConcrete->AddElement(Ca42,.000481);
  HDConcrete->AddElement(Ca43,.000103);
  HDConcrete->AddElement(Ca44,.00161);
  HDConcrete->AddElement(Ca46,.000003);
  HDConcrete->AddElement(Ca48,.000159);
  HDConcrete->AddElement(Ti46,.004171);
  HDConcrete->AddElement(Ti47,.003889);
  HDConcrete->AddElement(Ti48,.040153);
  HDConcrete->AddElement(Ti49,.003055);
  HDConcrete->AddElement(Ti50,.00306);
  HDConcrete->AddElement(V51,.003113);
  HDConcrete->AddElement(Cr50,.000071);
  HDConcrete->AddElement(Cr52,.001421);
  HDConcrete->AddElement(Cr53,.000164);
  HDConcrete->AddElement(Cr54,.000041);
  HDConcrete->AddElement(Mn55,.001981);
  HDConcrete->AddElement(Fe54,.027027);
  HDConcrete->AddElement(Fe56,.435692);
  HDConcrete->AddElement(Fe57,.010154);
  HDConcrete->AddElement(Fe58,.001378);


  // this is roughly NE213 liquid scint.
  G4double densityS = 0.874*g/cm3;
  G4Material* NE213 = new G4Material("NE213", densityS,2);

  G4double aH = 1.01*g/mole;
  G4Element* elH  = new G4Element("Hydrogen","H" ,1., aH);
  //might be able to replace with xylene?
  G4double aC = 12.01*g/mole;
  G4Element* elC  = new G4Element("Carbon"  ,"C" ,6., aC);
  NE213->AddElement(elC,0.5479);
  NE213->AddElement(elH,0.4521);

  



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
  G4LogicalVolume * universe_log = new G4LogicalVolume(universe_s,air,"universe_L",0,0,0);


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

  G4cout << water_rad+bioShield_thick-thermalC_length_big-thermalC_length_small-SSPlate_thick-gb_depth << 0. << -(worldZ/2.)+floor_thick_out+floor_to_TC_out << G4endl;

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
  G4LogicalVolume *floor_log = new G4LogicalVolume(floor,HDConcrete ,"floor_log");
  G4ThreeVector floorPlacement(0,0,-(worldZ/2.)+(floor_thick_out/2.));


  // graphite block in the ppol
  G4Box *graphite_block = new G4Box("graphite_block", gb_depth/2.,gb_width/2.,gb_height/2.);
  G4ThreeVector graphite_Placement(water_rad+bioShield_thick-thermalC_length_big-thermalC_length_small- (gb_depth/2.),0., -(worldZ/2.) + floor_thick_out + floor_to_TC_out);
  G4LogicalVolume *graphite_log = new G4LogicalVolume(graphite_block,nistManager->FindOrBuildMaterial("G4_GRAPHITE"),"graphite_lot");

  //water 
  // all the way up for now
  // two parts: semicircle and rectangular
  G4Box *pool_rec = new G4Box("pool_rec",worldX/4.,water_rad,(worldZ/2.) - floor_thick_pool/2.);
  G4Tubs *pool_tubs = new G4Tubs("pool_tubs",0.,water_rad,(worldZ/2.) - floor_thick_pool/2.,-pi/2.,pi);
  G4UnionSolid *poolI = new G4UnionSolid("poolI",pool_tubs,pool_rec,0,poolFloorTransX);
  G4Box *TC_pool_cutout = new G4Box("TC_pool_cutout",(thermalC_length_small)/2.,(thermalC_height_small)/2. + SSPlate_thick,(thermalC_height_small)/2. + SSPlate_thick);
  G4ThreeVector TC_in_poolPlacement(water_rad+bioShield_thick-thermalC_length_big-(thermalC_length_small/2.), 0, -((worldZ/2.) - floor_thick_pool/2.) + floor_to_TC_pool); // move to bottom of floor then up to proper height
  G4ThreeVector GB_in_poolPlacement(water_rad+bioShield_thick-thermalC_length_big-thermalC_length_small-gb_depth/2., 0, -((worldZ/2.) - floor_thick_pool/2.) + floor_to_TC_pool);
  G4SubtractionSolid *poolI2 = new G4SubtractionSolid("poolI2",poolI,TC_pool_cutout,0,TC_in_poolPlacement);
  G4SubtractionSolid *pool = new G4SubtractionSolid("pool",poolI2,graphite_block,0,GB_in_poolPlacement);
  G4ThreeVector poolPlacement(0,0,floor_thick_pool/2.);
  G4LogicalVolume *pool_log = new G4LogicalVolume(pool,nistManager->FindOrBuildMaterial("G4_WATER"),"pool_log");
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
  G4Box *TC_BS_cutout1 = new G4Box("TC_BS_cutout1",(thermalC_length_small)/2.,(thermalC_height_small)/2. + SSPlate_thick,(thermalC_height_small)/2. + SSPlate_thick);
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
  G4LogicalVolume *bioShield_log = new G4LogicalVolume(bioShield,HDConcrete,"bioShield_log");


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
  G4LogicalVolume *poolSS_log = new G4LogicalVolume(poolSS,StainlessSteel,"poolSS_log");

  // SS TC liner 
  // SS only in the small box part
  // there's some SS that goes into the HDC that I'm not including yet
  G4Box *TC_SS_box = new G4Box("TC_SS_box", thermalC_length_small/2.,(thermalC_height_small)/2. + SSPlate_thick,(thermalC_height_small)/2. + SSPlate_thick);
  G4Box *TC_SS_cutout = new G4Box("TC_SS_cutout", thermalC_length_small/2.,thermalC_height_small/2.,thermalC_height_small/2.);
  G4ThreeVector TC_SS_cutout_Trans(SSPlate_thick,0.,0.);
  G4SubtractionSolid *TC_SS = new G4SubtractionSolid("TC_SS",TC_SS_box,TC_SS_cutout,0,TC_SS_cutout_Trans);
  G4ThreeVector TC_SS_placement(water_rad+bioShield_thick-thermalC_length_big-(thermalC_length_small/2.),0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  G4LogicalVolume *TC_SS_log =  new G4LogicalVolume(TC_SS,StainlessSteel,"TC_SS_log");

  // Al TC liner 
  // Al only in the big box part
  G4Box *TC_Al_box = new G4Box("TC_Al_box", thermalC_length_big/2.,(thermalC_height_big)/2. + aluminumPlate_thick,(thermalC_height_big)/2. + aluminumPlate_thick);
  G4Box *TC_Al_cutout = new G4Box("TC_Al_cutout", thermalC_length_big/2.,thermalC_height_big/2.,thermalC_height_big/2.);
  G4ThreeVector TC_Al_cutout_Trans(aluminumPlate_thick,0.,0.);
  G4SubtractionSolid *TC_AlI = new G4SubtractionSolid("TC_AlI",TC_Al_box,TC_Al_cutout,0,TC_Al_cutout_Trans);
  G4ThreeVector TC_Al_cutout_Trans2(-aluminumPlate_thick,0.,0.);
  G4SubtractionSolid *TC_Al = new G4SubtractionSolid("TC_Al",TC_AlI,TC_SS_cutout,0,TC_Al_cutout_Trans2);
  G4ThreeVector TC_Al_placement(water_rad+bioShield_thick-thermalC_length_big/2.,0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  G4LogicalVolume *TC_Al_log =  new G4LogicalVolume(TC_Al,nistManager->FindOrBuildMaterial("G4_Al"),"TC_Al_log");

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
  G4LogicalVolume *door_Liner_log = new G4LogicalVolume(doorLiner,Boraflex ,"door_Liner_log");
  G4LogicalVolume *door_Lead_log = new G4LogicalVolume(doorLead,shieldLead ,"door_Lead_log");
  G4LogicalVolume *door_Concrete_log = new G4LogicalVolume(doorConcrete,HDConcrete ,"door_Concrete_log");
  G4LogicalVolume *door_Steel_log = new G4LogicalVolume(doorSteel,StainlessSteel ,"door_Steel_log");


  G4AssemblyVolume *door = new G4AssemblyVolume();
  G4RotationMatrix *zeroRot = new G4RotationMatrix;
  G4ThreeVector zeroPos;
  door->AddPlacedVolume(door_Liner_log,doorLinerPos,zeroRot);
  door->AddPlacedVolume(door_Lead_log,doorLeadPos,zeroRot);
  door->AddPlacedVolume(door_Concrete_log,doorConcretePos,zeroRot);
  door->AddPlacedVolume(door_Steel_log,doorSteelPos,zeroRot);
//  G4ThreeVector posDoor(water_rad+bioShield_thick,0.,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  // move door back for now
  G4ThreeVector posDoor(water_rad+bioShield_thick + 3.*12.*in,0.,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);


///*

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
  G4Box *poly_TC = new G4Box("poly_TC", TCPoly_thick/2.,thermalC_height_small/2.,thermalC_height_small/2.);
  G4LogicalVolume *polyTC_log = new G4LogicalVolume(poly_TC,HDPE,"polyTC_log");

  // reactor Front lead shielding
  G4Box *lead_TC = new G4Box("lead_TC",TCLead_thick/2.,thermalC_height_big/2.,thermalC_height_big/2.);
  G4LogicalVolume *leadTC_log = new G4LogicalVolume(lead_TC,shieldLead,"leadTC_log");


  // icebox oversimplified for now, three concentric cans with a can thickness of space between each
  G4Tubs *icebox1_out = new G4Tubs("icebox1_out", 0., iceboxRadius, iceboxHeight/2.,0,360*deg);
  G4Tubs *icebox1_in = new G4Tubs("icebox1_in", 0., iceboxRadius - iceboxThick, iceboxHeight/2. - iceboxThick,0,360*deg);
  G4SubtractionSolid *icebox1 = new G4SubtractionSolid("icebox1",icebox1_out,icebox1_in);
  G4LogicalVolume *icebox1_log = new G4LogicalVolume(icebox1,copper,"icebox1_log");

  G4Tubs *icebox2_out = new G4Tubs("icebox2_out", 0., iceboxRadius - 2*iceboxThick, iceboxHeight/2. - 2*iceboxThick,0,360*deg);
  G4Tubs *icebox2_in = new G4Tubs("icebox2_in", 0., iceboxRadius - 3*iceboxThick, iceboxHeight/2. - 3*iceboxThick,0,360*deg);
  G4SubtractionSolid *icebox2 = new G4SubtractionSolid("icebox2",icebox2_out,icebox2_in);
  G4LogicalVolume *icebox2_log = new G4LogicalVolume(icebox2,copper,"icebox2_log");

  G4Tubs *icebox3_out = new G4Tubs("icebox3_out", 0., iceboxRadius - 4*iceboxThick, iceboxHeight/2. - 4*iceboxThick,0,360*deg);
  G4Tubs *icebox3_in = new G4Tubs("icebox3_in", 0., iceboxRadius - 5*iceboxThick, iceboxHeight/2. - 5*iceboxThick,0,360*deg);
  G4SubtractionSolid *icebox3 = new G4SubtractionSolid("icebox3",icebox3_out,icebox3_in);
  G4LogicalVolume *icebox3_log = new G4LogicalVolume(icebox3,copper,"icebox3_log");

  // Let's put some detectors in the icebox, they can just float for now
  G4Tubs *det = new G4Tubs("det",0.,2*in,detThick/2.,0,360*deg);
  G4LogicalVolume *det1_log = new G4LogicalVolume(det, det_Ge,"det1_log");
  G4LogicalVolume *det2_log = new G4LogicalVolume(det, det_Ge,"det2_log");
  G4LogicalVolume *det3_log = new G4LogicalVolume(det, det_Ge,"det3_log");
  G4LogicalVolume *det4_log = new G4LogicalVolume(det, det_Ge,"det4_log");
  G4LogicalVolume *det5_log = new G4LogicalVolume(det, det_Si,"det5_log");
  G4LogicalVolume *det6_log = new G4LogicalVolume(det, det_Si,"det6_log");
  G4LogicalVolume *det7_log = new G4LogicalVolume(det, det_Si,"det7_log");
  G4LogicalVolume *det8_log = new G4LogicalVolume(det, det_Si,"det8_log");
  // should have slightly more than 3" for spacing of detectors

  // poly around the icebox
  G4Box *poly_IBshield_out = new G4Box("poly_IBshield_out",iceboxRadius+iceboxPoly_thick,iceboxRadius+iceboxPoly_thick,iceboxHeight/2. + iceboxPoly_thick);
  G4Box *poly_IBshield_in = new G4Box("poly_IBshield_in",iceboxRadius,iceboxRadius,iceboxHeight/2.);
  G4SubtractionSolid *poly_IBshield = new G4SubtractionSolid("poly_IBshield",poly_IBshield_out,poly_IBshield_in);
  G4LogicalVolume *polyIBshield_log = new G4LogicalVolume(poly_IBshield,HDPE,"poly_IBshield_log");

  // lead around the icebox
  G4Box *lead_IBshield_out = new G4Box("lead_IBshield_out",iceboxRadius+iceboxPoly_thick+iceboxLead_thick,iceboxRadius+iceboxPoly_thick+iceboxLead_thick,iceboxHeight/2. +iceboxPoly_thick+ iceboxLead_thick);
  G4Box *lead_IBshield_in = new G4Box("lead_IBshield_in",iceboxRadius+iceboxPoly_thick,iceboxRadius+iceboxPoly_thick,iceboxHeight/2. + iceboxPoly_thick);
  G4SubtractionSolid *lead_IBshield = new G4SubtractionSolid("lead_IBshield",lead_IBshield_out,lead_IBshield_in);
  G4LogicalVolume *leadIBshield_log = new G4LogicalVolume(lead_IBshield,shieldLead,"lead_IBshield_log");

  // poly around the IB lead shielding
  G4Box *poly_outershield_out = new G4Box("poly_outershield_out",iceboxRadius+iceboxPoly_thick+iceboxLead_thick,iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick,iceboxHeight/2. + iceboxPoly_thick+iceboxLead_thick+outerPoly_thick);
  //current design has only 4" on doorside and no extra on backside
  G4Box *poly_outershield_in = new G4Box("poly_outershield_in",iceboxRadius+iceboxPoly_thick+iceboxLead_thick,iceboxRadius+iceboxPoly_thick+iceboxLead_thick,iceboxHeight/2.+iceboxPoly_thick+iceboxLead_thick);
  G4SubtractionSolid *poly_outershield_shell = new G4SubtractionSolid("poly_outershield_shell",poly_outershield_out,poly_outershield_in);
  G4ThreeVector  poly_outer_back_place(iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_back_thick/2.,0,0);  
  G4Box *poly_outershield_back = new G4Box("poly_outershield_back",outerPoly_back_thick/2.,iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick,iceboxHeight/2.+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick);
  G4UnionSolid *poly_outershield = new G4UnionSolid("poly_outershield",poly_outershield_shell,poly_outershield_back,0,poly_outer_back_place);
  G4LogicalVolume *polyoutershield_log = new G4LogicalVolume(poly_outershield,HDPE,"poly_outershield_log");

  //scintillator veto, just top and bottom, this means there will be some extra space on the sides
  G4Box *muVeto = new G4Box("muVeto",(thermalC_length_big-aluminumPlate_thick-TCLead_thick)/2.,iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick,muVeto_thick/2.);
  G4LogicalVolume *muVetoTop_log = new G4LogicalVolume(muVeto,NE213 ,"muVetoTop_log");
  G4LogicalVolume *muVetoBottom_log = new G4LogicalVolume(muVeto,NE213 ,"muVetoBottom_log");
  //back scintillator
  G4Box *backScint = new G4Box("backScint",backScint_thick/2.,iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick,iceboxHeight/2.+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick);
  G4LogicalVolume *backScint_log = new G4LogicalVolume(backScint,NE213 ,"backScint_log");


  // Lead next to the liner on outside of thermal column
  G4Box *lead_outershield_in = new G4Box("lead_outershield_in",2*m,thermalC_height_big/2. -outerLead_thick,thermalC_height_big/2. -outerLead_thick);
  G4Box *lead_outershield_out = new G4Box("lead_outershield_out",(thermalC_length_big-aluminumPlate_thick-TCLead_thick)/2.,thermalC_height_big/2.,thermalC_height_big/2.);
  G4SubtractionSolid *lead_outershield = new G4SubtractionSolid("lead_outershield",lead_outershield_out,lead_outershield_in);
  G4LogicalVolume *leadoutershield_log = new G4LogicalVolume(lead_outershield,shieldLead,"lead_outershield_log");

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

//*/

/*
  // ************************************************** //
  // ************************************************** //
  // ***********    Shielding setup #2  *************** //
  // ************************************************** //
  // ************************************************** //

  // TC all but 8" poly
  // 8" lead barrier after
  // 2" poly surounding icebox
  // 2" lead after that
  // 5" poly top/bottom, 6" back
  // 2" muon veto top and bottom, 8" nova-like detector back
  G4double iceboxThick = 3./8.*in;  // ???
  G4double iceboxRadius = 6*in;
  G4double iceboxHeight = 18*in;
  G4double iceboxLead_thick = 2*in;
  G4double iceboxPoly_thick = 2*in;
  G4double muVeto_thick = 2*in;
  G4double backScint_thick = 8*in;
  G4double outerLead_thick = 2*in;
  G4double outerPoly_thick = 5*in - SSPlate_thick;
  G4double outerPoly_back_thick = 5*in;
  G4double TCLead_thick = 8*in;
  G4double TCPoly_thick = thermalC_length_small-largeOverburden_length - SSPlate_thick - 0.5*TCLead_thick;
  G4double detThick = 1.33*in;
  G4double detSpace = 0.4*in;

  // reactor poly shielding, completely filling the small part for now
  G4Box *poly_TC = new G4Box("poly_TC", TCPoly_thick/2.,thermalC_height_small/2.,thermalC_height_small/2.);
  G4LogicalVolume *polyTC_log = new G4LogicalVolume(poly_TC,HDPE,"polyTC_log");

  // reactor Front lead shielding
  G4Box *lead_TC = new G4Box("lead_TC",TCLead_thick/2.,thermalC_height_small/2.,thermalC_height_small/2.);
  G4LogicalVolume *leadTC_log = new G4LogicalVolume(lead_TC,shieldLead,"leadTC_log");


  // icebox oversimplified for now, three concentric cans with a can thickness of space between each
  G4Tubs *icebox1_out = new G4Tubs("icebox1_out", 0., iceboxRadius, iceboxHeight/2.,0,360*deg);
  G4Tubs *icebox1_in = new G4Tubs("icebox1_in", 0., iceboxRadius - iceboxThick, iceboxHeight/2. - iceboxThick,0,360*deg);
  G4SubtractionSolid *icebox1 = new G4SubtractionSolid("icebox1",icebox1_out,icebox1_in);
  G4LogicalVolume *icebox1_log = new G4LogicalVolume(icebox1,copper,"icebox1_log");

  G4Tubs *icebox2_out = new G4Tubs("icebox2_out", 0., iceboxRadius - 2*iceboxThick, iceboxHeight/2. - 2*iceboxThick,0,360*deg);
  G4Tubs *icebox2_in = new G4Tubs("icebox2_in", 0., iceboxRadius - 3*iceboxThick, iceboxHeight/2. - 3*iceboxThick,0,360*deg);
  G4SubtractionSolid *icebox2 = new G4SubtractionSolid("icebox2",icebox2_out,icebox2_in);
  G4LogicalVolume *icebox2_log = new G4LogicalVolume(icebox2,copper,"icebox2_log");

  G4Tubs *icebox3_out = new G4Tubs("icebox3_out", 0., iceboxRadius - 4*iceboxThick, iceboxHeight/2. - 4*iceboxThick,0,360*deg);
  G4Tubs *icebox3_in = new G4Tubs("icebox3_in", 0., iceboxRadius - 5*iceboxThick, iceboxHeight/2. - 5*iceboxThick,0,360*deg);
  G4SubtractionSolid *icebox3 = new G4SubtractionSolid("icebox3",icebox3_out,icebox3_in);
  G4LogicalVolume *icebox3_log = new G4LogicalVolume(icebox3,copper,"icebox3_log");

  // Let's put some detectors in the icebox, they can just float for now
  G4Tubs *det = new G4Tubs("det",0.,2*in,detThick/2.,0,360*deg);
  G4LogicalVolume *det1_log = new G4LogicalVolume(det, det_Ge,"det1_log");
  G4LogicalVolume *det2_log = new G4LogicalVolume(det, det_Ge,"det2_log");
  G4LogicalVolume *det3_log = new G4LogicalVolume(det, det_Ge,"det3_log");
  G4LogicalVolume *det4_log = new G4LogicalVolume(det, det_Ge,"det4_log");
  G4LogicalVolume *det5_log = new G4LogicalVolume(det, det_Si,"det5_log");
  G4LogicalVolume *det6_log = new G4LogicalVolume(det, det_Si,"det6_log");
  G4LogicalVolume *det7_log = new G4LogicalVolume(det, det_Si,"det7_log");
  G4LogicalVolume *det8_log = new G4LogicalVolume(det, det_Si,"det8_log");
  // should have slightly more than 3" for spacing of detectors

  // poly around the icebox
  G4Box *poly_IBshield_out = new G4Box("poly_IBshield_out",iceboxRadius+iceboxPoly_thick,iceboxRadius+iceboxPoly_thick,iceboxHeight/2. + iceboxPoly_thick);
  G4Box *poly_IBshield_in = new G4Box("poly_IBshield_in",iceboxRadius,iceboxRadius,iceboxHeight/2.);
  G4SubtractionSolid *poly_IBshield = new G4SubtractionSolid("poly_IBshield",poly_IBshield_out,poly_IBshield_in);
  G4LogicalVolume *polyIBshield_log = new G4LogicalVolume(poly_IBshield,HDPE,"poly_IBshield_log");

  // lead around the icebox
  G4Box *lead_IBshield_out = new G4Box("lead_IBshield_out",iceboxRadius+iceboxPoly_thick+iceboxLead_thick,iceboxRadius+iceboxPoly_thick+iceboxLead_thick,iceboxHeight/2. +iceboxPoly_thick+ iceboxLead_thick);
  G4Box *lead_IBshield_in = new G4Box("lead_IBshield_in",iceboxRadius+iceboxPoly_thick,iceboxRadius+iceboxPoly_thick,iceboxHeight/2. + iceboxPoly_thick);
  G4SubtractionSolid *lead_IBshield = new G4SubtractionSolid("lead_IBshield",lead_IBshield_out,lead_IBshield_in);
  G4LogicalVolume *leadIBshield_log = new G4LogicalVolume(lead_IBshield,shieldLead,"lead_IBshield_log");

  // poly around the IB lead shielding
  G4Box *poly_outershield_out = new G4Box("poly_outershield_out",iceboxRadius+iceboxPoly_thick+iceboxLead_thick,iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick,iceboxHeight/2. + iceboxPoly_thick+iceboxLead_thick+outerPoly_thick);
  //current design has only 4" on doorside and no extra on backside
  G4Box *poly_outershield_in = new G4Box("poly_outershield_in",iceboxRadius+iceboxPoly_thick+iceboxLead_thick,iceboxRadius+iceboxPoly_thick+iceboxLead_thick,iceboxHeight/2.+iceboxPoly_thick+iceboxLead_thick);
  G4SubtractionSolid *poly_outershield_shell = new G4SubtractionSolid("poly_outershield_shell",poly_outershield_out,poly_outershield_in);
  G4ThreeVector  poly_outer_back_place(iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_back_thick/2.,0,0);
  G4Box *poly_outershield_back = new G4Box("poly_outershield_back",outerPoly_back_thick/2.,iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick,iceboxHeight/2.+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick);
  G4UnionSolid *poly_outershield = new G4UnionSolid("poly_outershield",poly_outershield_shell,poly_outershield_back,0,poly_outer_back_place);
  G4LogicalVolume *polyoutershield_log = new G4LogicalVolume(poly_outershield,HDPE,"poly_outershield_log");

  //scintillator veto, just top and bottom, this means there will be some extra space on the sides
  G4Box *muVeto = new G4Box("muVeto",iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_back_thick/2.,iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick,muVeto_thick/2.);
  G4LogicalVolume *muVetoTop_log = new G4LogicalVolume(muVeto,NE213 ,"muVetoTop_log");
  G4LogicalVolume *muVetoBottom_log = new G4LogicalVolume(muVeto,NE213 ,"muVetoBottom_log");
  //back scintillator
  G4Box *backScint = new G4Box("backScint",backScint_thick/2.,iceboxRadius+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick,iceboxHeight/2.+iceboxPoly_thick+iceboxLead_thick+outerPoly_thick);
  G4LogicalVolume *backScint_log = new G4LogicalVolume(backScint,NE213 ,"backScint_log");


  // Lead next to the liner on outside of thermal column
  G4Box *lead_outershield_in = new G4Box("lead_outershield_in",2*m,thermalC_height_big/2. -outerLead_thick,thermalC_height_big/2. -outerLead_thick);
  G4Box *lead_outershield_out = new G4Box("lead_outershield_out",(thermalC_length_big-aluminumPlate_thick)/2.,thermalC_height_big/2.,thermalC_height_big/2.);
  G4SubtractionSolid *lead_outershield = new G4SubtractionSolid("lead_outershield",lead_outershield_out,lead_outershield_in);
  G4LogicalVolume *leadoutershield_log = new G4LogicalVolume(lead_outershield,shieldLead,"lead_outershield_log");

  // ok, let's connect everything into an assembly
  G4AssemblyVolume *shielding = new G4AssemblyVolume();

  G4ThreeVector   posShield(water_rad+bioShield_thick +SSPlate_thick  - thermalC_length_big - (thermalC_length_big - TCPoly_thick - .5*TCLead_thick) + iceboxRadius+iceboxPoly_thick+iceboxLead_thick, 0, -(worldZ/2.)+floor_thick_out+floor_to_TC_out);

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
  G4ThreeVector   leadOut_place(largeOverburden_length +(thermalC_length_big-aluminumPlate_thick-TCLead_thick)/2. - (iceboxRadius+iceboxPoly_thick+iceboxLead_thick ),0.,0.);
  G4ThreeVector   muVetoTop_place(outerPoly_back_thick/2.,0,(iceboxHeight/2. +iceboxPoly_thick+ iceboxLead_thick + outerPoly_thick+muVeto_thick/2.));
  G4ThreeVector   muVetoBot_place(outerPoly_back_thick/2.,0,-1*(iceboxHeight/2. +iceboxPoly_thick+ iceboxLead_thick + outerPoly_thick+muVeto_thick/2.));
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
  fLogicDet = new G4LogicalVolume(target, det_Ge,"Craig_log");

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


  G4LogicalVolume *fLogicDetTubeA = new G4LogicalVolume(detTubeA, nistManager->FindOrBuildMaterial("G4_Al"),"detTubeA_log");
  G4LogicalVolume *fLogicDetTubeB = new G4LogicalVolume(detTubeB, nistManager->FindOrBuildMaterial("G4_Al"),"detTubeB_log");
  G4LogicalVolume *fLogicDetTubeC = new G4LogicalVolume(detTubeC, StainlessSteel,"detTubeC_log");
  G4LogicalVolume *fLogicDetTubeD = new G4LogicalVolume(detTubeD, StainlessSteel,"detTubeD_log");
  G4LogicalVolume *fLogicDetLiN2A = new G4LogicalVolume(liqN2_TubeA,liquidNitrogen,"liqN2_TubeA_log");
  G4LogicalVolume *fLogicDetLiN2B = new G4LogicalVolume(liqN2_TubeB,liquidNitrogen,"liqN2_TubeB_log");
  G4LogicalVolume *fLogicDetLiN2C = new G4LogicalVolume(liqN2_TubeC,liquidNitrogen,"liqN2_TubeC_log");
  G4LogicalVolume *fLogicDetLiN2D = new G4LogicalVolume(liqN2_TubeD,liquidNitrogen,"liqN2_TubeD_log");



  fullDet->AddPlacedVolume(fLogicDetTubeA,tubeAPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetTubeB,tubeBPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetTubeC,tubeCPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetTubeD,tubeDPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetLiN2A,tubeA_N2_Pos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetLiN2B,tubeB_N2_Pos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetLiN2C,tubeCPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDetLiN2D,tubeDPos,zeroRot);
  fullDet->AddPlacedVolume(fLogicDet,zeroPos,zeroRot);


  G4ThreeVector posDet = G4ThreeVector(water_rad+bioShield_thick + 12.*3.*in,0.,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);  // 3 ft outside thermal column


  // collecting placements here to make it easy to turn stuff off
  fUniverse_phys = new G4PVPlacement(0,G4ThreeVector(),"universe_P",universe_log,0,false,0);
  new G4PVPlacement(0,floorPlacement,floor_log,"floor_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,graphite_Placement,graphite_log,"graphite_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,poolPlacement,pool_log,"pool_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,bioShieldPlacement,bioShield_log,"bioShield_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,poolSSPlacement,poolSS_log,"poolSS_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TC_SS_placement,TC_SS_log,"TC_SS_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TC_Al_placement,TC_Al_log,"TC_Al_phys",universe_log,false,0,fCheckOverlaps);
  door->MakeImprint(universe_log,posDoor,zeroRot,0,fCheckOverlaps);
//  fullDet->MakeImprint(universe_log,posDet,RotDet,0,fCheckOverlaps);
  shielding->MakeImprint(universe_log,posShield,zeroRot,0,fCheckOverlaps);





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


  icebox1_log->SetVisAttributes(aVisAttCu);
  icebox2_log->SetVisAttributes(aVisAttCu);
  icebox3_log->SetVisAttributes(aVisAttCu);
  polyTC_log->SetVisAttributes(aVisAttPoly);
  leadTC_log->SetVisAttributes(aVisAttLead);
  det1_log->SetVisAttributes(aVisAttDet);
  det2_log->SetVisAttributes(aVisAttDet);
  det3_log->SetVisAttributes(aVisAttDet);
  det4_log->SetVisAttributes(aVisAttDet);
  det5_log->SetVisAttributes(aVisAttDet);
  det6_log->SetVisAttributes(aVisAttDet);
  det7_log->SetVisAttributes(aVisAttDet);
  det8_log->SetVisAttributes(aVisAttDet);
  polyIBshield_log->SetVisAttributes(aVisAttPoly);
  leadIBshield_log->SetVisAttributes(aVisAttLead);
  polyoutershield_log->SetVisAttributes(aVisAttPoly);
  leadoutershield_log->SetVisAttributes(aVisAttLead);
  muVetoTop_log->SetVisAttributes(aVisAttScint);
  muVetoBottom_log->SetVisAttributes(aVisAttScint);
  backScint_log->SetVisAttributes(aVisAttScint);

  return fUniverse_phys;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
