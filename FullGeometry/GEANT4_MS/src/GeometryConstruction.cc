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

  G4String vacName = "/AntiCompton/Vac";
  G4String vacHits = "ACVacHitsCollection";

  MinerSD* vac1SD = new MinerSD(vacName+"1", vacHits+"1");
  SDman->AddNewDetector(vac1SD);
  SetSensitiveDetector("Vac1_log", vac1SD, true);

  MinerSD* vac2SD = new MinerSD(vacName+"2", vacHits+"2");
  SDman->AddNewDetector(vac2SD);
  SetSensitiveDetector("Vac2_log", vac2SD, true);

  MinerSD* vac3SD = new MinerSD(vacName+"3", vacHits+"3");
  SDman->AddNewDetector(vac3SD);
  SetSensitiveDetector("Vac3_log", vac3SD, true);

  MinerSD* vac4SD = new MinerSD(vacName+"4", vacHits+"4");
  SDman->AddNewDetector(vac4SD);
  SetSensitiveDetector("Vac4_log", vac4SD, true);

  G4String detSDname = "/AntiCompton/Zip";
  G4String hitsName = "ACZipHitsCollection";
  MinerSD* aZipSD = new MinerSD(detSDname, hitsName);
  SDman->AddNewDetector(aZipSD);
  SetSensitiveDetector("Zip_log", aZipSD, true);

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

  // this is temporary until we have better understanding of HD concrete composition
  //G4double concreteDensity = 3.53*g/cm3;
  //G4Material* HDConcrete = nistManager->BuildMaterialWithNewDensity("HDConcrete","G4_CONCRETE",concreteDensity);

  //G4_POLYETHYLENE  (C_2H_4)_N-Polyethylene  0.94 g/cm3,  
  //             1     0.143711
  //             6     0.856289
  G4Material *boron   = nistManager->FindOrBuildMaterial("G4_B");       
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
  G4Material* fScintMat2 = new G4Material("NE213", densityS,2);

  G4double aH = 1.01*g/mole;
  G4Element* elH  = new G4Element("Hydrogen","H" ,1., aH);
  //might be able to replace with xylene?
  G4double aC = 12.01*g/mole;
  G4Element* elC  = new G4Element("Carbon"  ,"C" ,6., aC);
  fScintMat2->AddElement(elC,0.5479);
  fScintMat2->AddElement(elH,0.4521);

  



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
  G4LogicalVolume * universe_log = new G4LogicalVolume(universe_s,vac,"universe_L",0,0,0);


  G4double bioShield_thick = 1.7*m; // this includes SS plating which is why you'll see it subtracted out below
  G4double water_rad = 12.*4.5*in;  // from Sean's talk at TAMU CNS workshop
  G4double floor_to_TC_out = 48*in; // distance from ground floor to center of thermal column, ESTIMATE
  G4double floor_to_TC_pool = 60*in; // distance from pool floor to center of thermal column 
  G4double floor_thick_out = 0.5*m; // reasonable?
  G4double floor_thick_pool = (floor_thick_out-(floor_to_TC_pool-floor_to_TC_out)); // this includes SS plating which is why you'll see it subtracted out below
  G4double thermalC_height_big = 1.24*m; // larger part of thermal column, square face side
  G4double thermalC_height_small = 1.0*m; // smaller part of thermal column, square face side
  G4double thermalC_length_big = 1.3*m;
  G4double thermalC_length_small = 1.3*m;
  G4double bioShield_flat_width = thermalC_height_big + 1*m; // this is the width of the "flat" part of the bio shield where the thermal column is ESTIMATE
  G4double poolC_to_TC = water_rad - (25.8125+8.875)*in; // this is the distance from the center of the circular part of the pool to the back of the thermal column
  G4double aluminumPlate_thick = 0.5*in; // THIS IS A GUESS
  G4double SSPlate_thick = 0.5*in; // THIS IS A GUESS
  G4double gb_depth = 26.5625*in; // depth of graphite block
  G4double gb_width = 29.5*in; // width of graphite block ESTIMATE
  G4double gb_height = 29.5*in; // height of graphite block


  G4cout << water_rad+bioShield_thick-thermalC_length_big-thermalC_length_small-SSPlate_thick << 0. << -(worldZ/2.)+floor_thick_out+floor_to_TC_out << G4endl;

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
  G4Box *TC_pool_cutout = new G4Box("TC_pool_cutout",(thermalC_length_small+SSPlate_thick)/2.,(thermalC_height_small+SSPlate_thick)/2.,(thermalC_height_small+SSPlate_thick)/2.);
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
  G4Box *TC_BS_cutout1 = new G4Box("TC_BS_cutout1",(thermalC_length_small+SSPlate_thick)/2.,(thermalC_height_small+SSPlate_thick)/2.,(thermalC_height_small+SSPlate_thick)/2.);
  G4Box *TC_BS_cutout2 = new G4Box("TC_BS_cutout2",(thermalC_length_big+aluminumPlate_thick)/2.,(thermalC_height_big+aluminumPlate_thick)/2.,(thermalC_height_big+aluminumPlate_thick)/2.);
  //there's still some HDC left because of the curvature of the wall so cut that out as well
  G4Box *TC_BS_cutout3 = new G4Box("TC_BS_cutout3",(thermalC_length_big+SSPlate_thick)/2.,(thermalC_height_big+SSPlate_thick)/2.,(thermalC_height_big+SSPlate_thick)/2.);
  G4ThreeVector TCsmall_in_BSPlacement(water_rad+bioShield_thick-thermalC_length_big-(thermalC_length_small/2.),0,-((worldZ/2.)-(floor_thick_out/2.))+floor_to_TC_out);
  G4ThreeVector TCbig_in_BSPlacement(water_rad+bioShield_thick-(thermalC_length_big/2.),0,-((worldZ/2.)-(floor_thick_out/2.))+floor_to_TC_out);
  G4ThreeVector TCleftover_in_BSPlacement(water_rad+bioShield_thick,0,-((worldZ/2.)-(floor_thick_out/2.))+floor_to_TC_out);
  G4SubtractionSolid * TC_in_BS1 = new G4SubtractionSolid("TC_in_BS1",bioShieldI2,TC_BS_cutout1,0,TCsmall_in_BSPlacement);
  G4SubtractionSolid * TC_in_BS2 = new G4SubtractionSolid("TC_in_BS2",TC_in_BS1,TC_BS_cutout2,0,TCbig_in_BSPlacement);
  G4SubtractionSolid * bioShieldI3 = new G4SubtractionSolid("bioShieldI3",TC_in_BS2,TC_BS_cutout3,0,TCleftover_in_BSPlacement);
  G4double sagitta = (water_rad + bioShield_thick) - sqrt(pow((water_rad + bioShield_thick),2) - pow((bioShield_flat_width/2.),2));
  G4cout << sagitta << G4endl;
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
  G4ThreeVector TC_in_SSPlacement(water_rad+bioShield_thick-thermalC_length_big-(thermalC_length_small/2.), 0, -((worldZ/2.) - floor_thick_pool/2.) + floor_to_TC_pool); // move to bottom of floor then up to proper height
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
  G4Box *TC_SS_box = new G4Box("TC_SS_box", thermalC_length_small/2.,(thermalC_height_small+SSPlate_thick)/2.,(thermalC_height_small+SSPlate_thick)/2.);
  G4Box *TC_SS_cutout = new G4Box("TC_SS_cutout", thermalC_length_small/2.,thermalC_height_small/2.,thermalC_height_small/2.);
  G4ThreeVector TC_SS_cutout_Trans(SSPlate_thick,0.,0.);
  G4SubtractionSolid *TC_SS = new G4SubtractionSolid("TC_SS",TC_SS_box,TC_SS_cutout,0,TC_SS_cutout_Trans);
  G4ThreeVector TC_SS_placement(water_rad+bioShield_thick-thermalC_length_big-(thermalC_length_small/2.),0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  G4LogicalVolume *TC_SS_log =  new G4LogicalVolume(TC_SS,StainlessSteel,"TC_SS_log");

  // Al TC liner 
  // Al only in the big box part
  G4Box *TC_Al_box = new G4Box("TC_Al_box", thermalC_length_big/2.,(thermalC_height_big+aluminumPlate_thick)/2.,(thermalC_height_big+aluminumPlate_thick)/2.);
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
  G4LogicalVolume *door_Lead_log = new G4LogicalVolume(doorLead,nistManager->FindOrBuildMaterial("G4_Pb") ,"door_Lead_log");
  G4LogicalVolume *door_Concrete_log = new G4LogicalVolume(doorConcrete,HDConcrete ,"door_Concrete_log");
  G4LogicalVolume *door_Steel_log = new G4LogicalVolume(doorSteel,StainlessSteel ,"door_Steel_log");


  G4AssemblyVolume *door = new G4AssemblyVolume();
  G4RotationMatrix *zeroRot = new G4RotationMatrix;
  G4ThreeVector zeroPos;
  door->AddPlacedVolume(door_Liner_log,doorLinerPos,zeroRot);
  door->AddPlacedVolume(door_Lead_log,doorLeadPos,zeroRot);
  door->AddPlacedVolume(door_Concrete_log,doorConcretePos,zeroRot);
  door->AddPlacedVolume(door_Steel_log,doorSteelPos,zeroRot);
  G4ThreeVector posDoor(water_rad+bioShield_thick,0.,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);



  //these are dummy scoring meshes essentially until I learn how to do it for reals
  // these are overlapping something FIX IT FIX IT FIX IT
  G4Box * vac1 = new G4Box("vac1",0.0001*m,thermalC_height_small/2.,thermalC_height_small/2.);
  G4Box * vac2 = new G4Box("vac2",0.0001*m,thermalC_height_small/2.,thermalC_height_small/2.);
  G4Box * vac3 = new G4Box("vac3",0.0001*m,thermalC_height_big/2.,thermalC_height_big/2.);
  G4Box * vac4 = new G4Box("vac4",0.0001*m,thermalC_height_big/2.,thermalC_height_big/2.);

  G4ThreeVector vac1_placement(water_rad+bioShield_thick+SSPlate_thick-thermalC_length_big-thermalC_length_small + 0.0001*m,0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  G4ThreeVector vac2_placement(water_rad+bioShield_thick+SSPlate_thick-thermalC_length_big + 0.0001*m,0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  G4ThreeVector vac3_placement(water_rad+bioShield_thick-0.0002*m,0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);
  G4ThreeVector vac4_placement(water_rad+bioShield_thick+((3./8.)*in) + (4.*in) + (20.*in) + ((1./16.)*in) + 0.0002*m ,0,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);  

  G4LogicalVolume *Vac1_log = new G4LogicalVolume(vac1,vac,"Vac1_log");
  G4LogicalVolume *Vac2_log = new G4LogicalVolume(vac2,vac,"Vac2_log");
  G4LogicalVolume *Vac3_log = new G4LogicalVolume(vac3,vac,"Vac3_log");
  G4LogicalVolume *Vac4_log = new G4LogicalVolume(vac4,vac,"Vac4_log");



  G4ThreeVector posDet = G4ThreeVector(water_rad+bioShield_thick + 12.*3.*in,0.,-(worldZ/2.)+floor_thick_out+floor_to_TC_out);  // 3 ft outside thermal column

  // Craig's Ge detector
  // Need to put N2 into the detector still

  G4double detRad = 50.8*mm; // 4" diameter
  G4double detHalfZ = 16.891*mm; // 1.33" thick 

  G4VSolid* target = new G4Tubs("Target",0.,2.5*cm,(4.805/2.)*cm,0,360*deg);
  G4RotationMatrix* RotDet = new G4RotationMatrix;
  RotDet->rotateZ(-pi/4.);
  RotDet->rotateX(0.);
  RotDet->rotateY(-pi/2.);
  fLogicDet = new G4LogicalVolume(target, nistManager->FindOrBuildMaterial("G4_Ge"),"Zip_log");

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




  // collecting placements here to make it easy to turn stuff off
  fUniverse_phys = new G4PVPlacement(0,G4ThreeVector(),"universe_P",universe_log,0,false,0);
  new G4PVPlacement(0,floorPlacement,floor_log,"floor_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,graphite_Placement,graphite_log,"graphite_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,poolPlacement,pool_log,"pool_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,bioShieldPlacement,bioShield_log,"bioShield_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,poolSSPlacement,poolSS_log,"poolSS_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TC_SS_placement,TC_SS_log,"TC_SS_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TC_Al_placement,TC_Al_log,"TC_Al_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,vac1_placement,Vac1_log,"vac1_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,vac2_placement,Vac2_log,"vac2_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,vac3_placement,Vac3_log,"vac3_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,vac4_placement,Vac4_log,"vac4_phys",universe_log,false,0,fCheckOverlaps);
  door->MakeImprint(universe_log,posDoor,zeroRot,0,fCheckOverlaps);
  fullDet->MakeImprint(universe_log,posDet,RotDet,0,fCheckOverlaps);





//--------- Visualization attributes -------------------------------
  universe_log->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* aVisAttWater = new G4VisAttributes(G4Colour(0.0,0,1.0));
  G4VisAttributes* aVisAttGraphite = new G4VisAttributes(G4Colour(.2,0.2,.2));
  G4VisAttributes* aVisAttLead = new G4VisAttributes(G4Colour(.1,0.1,.1));
  G4VisAttributes* aVisAttPoly = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  G4VisAttributes* aVisAttAlTubes = new G4VisAttributes(G4Colour(.5,.5,.5));
  G4VisAttributes* aVisAttStainless = new G4VisAttributes(G4Colour(.5,.5,.5));
  G4VisAttributes* aVisAttHDConcrete = new G4VisAttributes(G4Colour(0.3,0.3,.3));
  G4VisAttributes* aVisAttBioShield = new G4VisAttributes(G4Colour(0.55,0.55,.55));
  G4VisAttributes* aVisAttSSLiner = new G4VisAttributes(G4Colour(0.4,0.4,.4));
  G4VisAttributes* aVisAttLiqN2 = new G4VisAttributes(G4Colour(0.0,0.0,.8));


  aVisAttAlTubes->SetForceWireframe(true);
  aVisAttWater->SetForceWireframe(true);
  aVisAttBioShield->SetForceWireframe(true);
//  aVisAttHDConcrete->SetForceWireframe(true);
//  aVisAttStainless->SetForceWireframe(true);
  aVisAttSSLiner->SetForceWireframe(true);
  aVisAttLiqN2->SetForceWireframe(true);

  pool_log->SetVisAttributes(aVisAttWater);
  floor_log->SetVisAttributes(aVisAttHDConcrete);
  bioShield_log->SetVisAttributes(aVisAttBioShield);
  poolSS_log->SetVisAttributes(aVisAttSSLiner);

  graphite_log->SetVisAttributes(aVisAttGraphite);

  TC_SS_log->SetVisAttributes(aVisAttSSLiner);
  TC_Al_log->SetVisAttributes(aVisAttSSLiner);


//  Vac1_log->SetVisAttributes(aVisAttSSLiner);
//  Vac2_log->SetVisAttributes(aVisAttSSLiner);
//  Vac3_log->SetVisAttributes(aVisAttSSLiner);
//  Vac4_log->SetVisAttributes(aVisAttSSLiner);

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


  return fUniverse_phys;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
