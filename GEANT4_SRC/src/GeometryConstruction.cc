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
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


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

  MinerSD* waterSD = new MinerSD("/AntiCompton/Water", "ACWaterHitsCollection");
  SDman->AddNewDetector(waterSD);
  SetSensitiveDetector("Water_log", waterSD, true);

  MinerSD* graphiteSD = new MinerSD("/AntiCompton/Graphite", "ACGraphiteHitsCollection");
  SDman->AddNewDetector(graphiteSD);
  SetSensitiveDetector("Graphite_log", graphiteSD, true);

  MinerSD* concreteSD = new MinerSD("/AntiCompton/Concrete", "ACConcreteHitsCollection");
  SDman->AddNewDetector(concreteSD);
  SetSensitiveDetector("Concrete_log", concreteSD, true);

  MinerSD* leadSD = new MinerSD("/AntiCompton/Lead", "ACLeadHitsCollection");
  SDman->AddNewDetector(leadSD);
  SetSensitiveDetector("Lead_log", leadSD, true);

  MinerSD* polySD = new MinerSD("/AntiCompton/Poly", "ACPolyHitsCollection");
  SDman->AddNewDetector(polySD);
  SetSensitiveDetector("Poly_log", polySD, true);

  G4String detSDname = "/AntiCompton/Zip";
  G4String hitsName = "ACZipHitsCollection";
  MinerSD* aZipSD = new MinerSD(detSDname, hitsName);
  SDman->AddNewDetector(aZipSD);
  SetSensitiveDetector("Zip_log", aZipSD, true);

  MinerSD* scintSD = new MinerSD("/AntiCompton/ScintBox", "ACScintHitsCollection");
  SDman->AddNewDetector(scintSD);
  SetSensitiveDetector("ScintBox_log", scintSD, true);

}


G4VPhysicalVolume* GeometryConstruction::Construct()
{

  G4NistManager* nistManager = G4NistManager::Instance();
  G4double in = 2.54*cm;
  
  G4bool fCheckOverlaps = true;
  G4Material* vac = new G4Material("Vacuum", 1., 1.01*g/mole, universe_mean_density,
                                 kStateGas,0.000017*kelvin,1.e-19*pascal);


  G4double world_r = 20*m;
  G4Sphere * universe_s  = new G4Sphere("universe_s", 0, world_r, 0, twopi, 0, pi);
  G4LogicalVolume * universe_log = new G4LogicalVolume(universe_s,vac,"universe_L",0,0,0);
  fUniverse_phys = new G4PVPlacement(0,G4ThreeVector(),"universe_P",universe_log,0,false,0);


  G4double detRad = 50.8*mm; // 4" diameter
  G4double detHalfZ = 16.891*mm; // 1.33" thick 

  // Andrew's Ge detector
  //G4double detRadA = (1./sqrt(pi))*cm; // from chamber center to center!
  //G4double detHalfZA = 0.5*cm;

  // Craig's Ge detector
  G4double detRadA = 2.5*cm; // from chamber center to center!
  G4double detHalfZA = (4.805/2.)*cm;


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

  G4double leadThick = 0.000000001*m;
  G4double polyThick = 0.000000001*m;
  G4double graphiteThick = 0.000000001*m;
  G4double waterThick = 0.000000001*m;
  G4double concreteThick = 0.000000001*m;


  //G4double leadThick = 0.2032*m;
  //G4double polyThick = 0.2032*m;
  //G4double graphiteThick = 0.6604*m;
  //G4double waterThick = 0.5*m;

  //G4double leadThick = 24*in;
  //G4double polyThick = 0.4064*m;
  //G4double graphiteThick = 0.6604*m;
  //G4double waterThick = 0.5*m;
  //G4double concreteThick = 10*in;

  //G4double leadThick = 8.0*in;


  // for Andrew's pool wall setup
  //G4double waterThick = 24*in;
  //G4double concreteThick = 6*12*in;


  G4double cavernSide = 1.2*m;

  G4Box * vac1 = new G4Box("vac1",0.01*m,cavernSide/2.,cavernSide/2.);
  G4Box * vac2 = new G4Box("vac2",0.01*m,cavernSide/2.,cavernSide/2.);
  G4Box * water = new G4Box("water",waterThick/2.,cavernSide/2.,cavernSide/2.);
  G4Box * graphite = new G4Box("graphite",graphiteThick/2.,cavernSide/2.,cavernSide/2.);
  G4Box * lead = new G4Box("lead",leadThick/2.,cavernSide/2.,cavernSide/2.);
  G4Box * poly = new G4Box("poly",polyThick/2.,cavernSide/2.,cavernSide/2.);
  G4Box * concrete = new G4Box("concrete",concreteThick/2.,cavernSide/2.,cavernSide/2.);

  fLogicVac1 = new  G4LogicalVolume(vac1, vac , "Vac1_log");
  fLogicVac2 = new  G4LogicalVolume(vac2, vac , "Vac2_log");
  fLogicWater = new G4LogicalVolume(water, nistManager->FindOrBuildMaterial("G4_WATER") , "Water_log");
  fLogicGraphite = new G4LogicalVolume(graphite, nistManager->FindOrBuildMaterial("G4_GRAPHITE") , "Graphite_log");
  fLogicLead = new G4LogicalVolume(lead,nistManager->FindOrBuildMaterial("G4_Pb") , "Lead_log");

  //fLogicPoly = new G4LogicalVolume(poly,nistManager->FindOrBuildMaterial("G4_POLYETHYLENE") , "Poly_log");
  fLogicPoly = new G4LogicalVolume(poly,boratedpoly05 , "Poly_log");

  fLogicConcrete = new G4LogicalVolume(concrete,HDConcrete , "Concrete_log");

  G4ThreeVector zTransV1(-5*m - waterThick/2. - 0.01*m,0,0);
  G4ThreeVector zTransW(-5.*m,0,0);
  G4ThreeVector zTransG(-5.*m + waterThick/2. + graphiteThick/2.,0,0);
  G4ThreeVector zTransC(-5.*m + waterThick/2. + graphiteThick + concreteThick/2.,0,0);
  G4ThreeVector zTransV2(-5*m + waterThick/2. + graphiteThick + concreteThick + leadThick + polyThick + 0.01*m,0,0);

  //lead then poly
  G4ThreeVector zTransL(-5.*m + waterThick/2. + graphiteThick + concreteThick + leadThick/2.,0,0);
  G4ThreeVector zTransP(-5.*m + waterThick/2. + graphiteThick + concreteThick + leadThick + polyThick/2.,0,0);

  //poly then lead
  //G4ThreeVector zTransP(-5.*m + waterThick/2. + graphiteThick + concreteThick + polyThick/2.,0,0);
  //G4ThreeVector zTransL(-5.*m + waterThick/2. + graphiteThick + concreteThick + polyThick + leadThick/2.,0,0);

  /*
  new G4PVPlacement(0,zTransV1,fLogicVac1,"vac1_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,zTransV2,fLogicVac2,"vac2_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,zTransW,fLogicWater,"water_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,zTransG,fLogicGraphite,"graphite_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,zTransC,fLogicConcrete,"concrete_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,zTransL,fLogicLead,"lead_phys",universe_log,false,0,fCheckOverlaps);
  new G4PVPlacement(0,zTransP,fLogicPoly,"poly_phys",universe_log,false,0,fCheckOverlaps);
  */

//  G4ThreeVector positionTarget = G4ThreeVector(-5.*m + waterThick/2. + graphiteThick + concreteThick + leadThick + polyThick + detHalfZA + 5*cm,0,0);

//  G4ThreeVector positionTarget = G4ThreeVector(-5.*m + waterThick/2. + graphiteThick + concreteThick + leadThick + polyThick + detRad + 24*in,0,0);

  G4ThreeVector positionTarget = G4ThreeVector(0,0,0);


  //G4VSolid* zip1 = new G4Tubs("Zip1",0.,detRad,detHalfZ,0,360*deg);

  // This is Craig's
  ///*
  G4VSolid* target = new G4Tubs("Target",0.,detRadA,detHalfZA,0,360*deg);
  G4RotationMatrix* RotDet = new G4RotationMatrix;
  RotDet->rotateZ(pi/2.);
  RotDet->rotateX(pi/2.);
  RotDet->rotateY(0);
  fLogicDet = new G4LogicalVolume(target, nistManager->FindOrBuildMaterial("G4_Ge"),"Zip_log");
  new G4PVPlacement(RotDet,               // no rotation
                    positionTarget,  // at (x,y,z)
                    fLogicDet,    // its logical volume
                    "Zip_phys",        // its name
                    universe_log,         // its mother volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps 
  //*/

  /*
  G4VSolid* zip1 = new G4Tubs("Zip1",0.,detRad,detHalfZ,0,360*deg);
  //fLogicDet = new G4LogicalVolume(zip1, nistManager->FindOrBuildMaterial("G4_Si"),"Zip_log");
  fLogicDet = new G4LogicalVolume(zip1, nistManager->FindOrBuildMaterial("G4_Ge"),"Zip_log");
//  new G4PVPlacement(0,               // no rotation
//                    positionTarget,  // at (x,y,z)
//                    fLogicDet,    // its logical volume
//                    "Zip_phys",        // its name
//                    universe_log,         // its mother volume
//                    false,           // no boolean operations
//                    0,               // copy number
//                    fCheckOverlaps); // checking overlaps 
  */


  ///*
  // 8 " Lead housing for Andrew's detector
  G4Box *bigBox = new G4Box("bigBox", 8*in, 8*in, 6*in);
  G4Tubs *subCyl = new G4Tubs("subCyl", 0., 0.5*in, 2.*in,0,360*deg);
  G4ThreeVector cavPos(0,0,-4.*in);
  G4SubtractionSolid *leadHouse = new G4SubtractionSolid("leadHouse",bigBox,subCyl,0,cavPos);

  fLogicLeadHouse = new G4LogicalVolume(leadHouse, nistManager->FindOrBuildMaterial("G4_Pb"),"LeadHouse_log");
  G4ThreeVector zTransLH(0,0,2*in);
  //new G4PVPlacement(0,zTransLH,fLogicLeadHouse,"leadHouse_phys",universe_log,false,0,fCheckOverlaps);
  //*/

  // 4 " Lead housing for Andrew's detector
  /*
  G4Box *bigBox = new G4Box("bigBox", 4*in, 4*in, 4*in);
  G4Tubs *subCyl = new G4Tubs("subCyl", 0., 0.5*in, 2.*in,0,360*deg);
  G4ThreeVector cavPos(0,0,-2.*in);
  G4SubtractionSolid *leadHouse = new G4SubtractionSolid("leadHouse",bigBox,subCyl,0,cavPos);

  fLogicLeadHouse = new G4LogicalVolume(leadHouse, nistManager->FindOrBuildMaterial("G4_Pb"),"LeadHouse_log");
  G4ThreeVector zTransLH(0,0,0);
  new G4PVPlacement(0,zTransLH,fLogicLeadHouse,"leadHouse_phys",universe_log,false,0,fCheckOverlaps);
  */

  // Test a liquid scintillator box around detector
  G4double densityS = 0.874*g/cm3;
  G4Material* fScintMat2 = new G4Material("NE213", densityS,2);

  G4double aH = 1.01*g/mole;
  G4Element* elH  = new G4Element("Hydrogen","H" ,1., aH);

  G4double aC = 12.01*g/mole;
  G4Element* elC  = new G4Element("Carbon"  ,"C" ,6., aC);
  fScintMat2->AddElement(elC,0.5479);
  fScintMat2->AddElement(elH,0.4521);

  G4Box *scintBoxO = new G4Box("scintBoxO", 7.*in, 7.*in, 7.*in);
  G4Box *scintBoxI = new G4Box("scintBoxI", 5.*in, 5.*in, 5.*in);
  G4SubtractionSolid *scintBox = new G4SubtractionSolid("scintBox",scintBoxO,scintBoxI);
  fLogicScintBox = new G4LogicalVolume(scintBox, fScintMat2, "ScintBox_log");
  //new G4PVPlacement(0,positionTarget,fLogicScintBox,"ScintBox_phys",universe_log,false,0,fCheckOverlaps);


  G4Box *leadBoxO = new G4Box("scintBoxO", 14.*in, 14.*in, 14.*in);
  G4Box *leadBoxI = new G4Box("scintBoxI", 6.*in, 6.*in, 6.*in);
  G4SubtractionSolid *leadBox = new G4SubtractionSolid("leadBox",leadBoxO,leadBoxI);
  fLogicLeadBox = new G4LogicalVolume(leadBox, nistManager->FindOrBuildMaterial("G4_Pb"), "LeadBox_log");
  //new G4PVPlacement(0,positionTarget,fLogicLeadBox,"LeadBox_phys",universe_log,false,0,fCheckOverlaps);

  G4Box *polyBoxO = new G4Box("polyBoxO", 20.*in, 20.*in, 20.*in);
  G4Box *polyBoxI = new G4Box("polyBoxI", 14.*in, 14.*in, 14.*in);
  G4SubtractionSolid *polyBox = new G4SubtractionSolid("polyBox",polyBoxO,polyBoxI);
  fLogicPolyBox = new G4LogicalVolume(polyBox,nistManager->FindOrBuildMaterial("G4_POLYETHYLENE") , "PolyBox_log");
  //new G4PVPlacement(0,positionTarget,fLogicPolyBox,"PolyBox_phys",universe_log,false,0,fCheckOverlaps);


//--------- Visualization attributes -------------------------------
  universe_log->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* aVisAttWater = new G4VisAttributes(G4Colour(0.0,0,1.0));
  G4VisAttributes* aVisAttGraphite = new G4VisAttributes(G4Colour(.5,0.5,.5));
  G4VisAttributes* aVisAttLead = new G4VisAttributes(G4Colour(1.0,0,1.0));
  G4VisAttributes* aVisAttPoly = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));


  fLogicWater->SetVisAttributes(aVisAttWater);
  fLogicGraphite->SetVisAttributes(aVisAttGraphite);
  fLogicLead->SetVisAttributes(aVisAttLead);
  fLogicPoly->SetVisAttributes(aVisAttPoly);

  return fUniverse_phys;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
