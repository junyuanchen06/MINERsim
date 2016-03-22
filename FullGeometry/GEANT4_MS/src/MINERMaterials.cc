#include "MINERMaterials.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


MINERMaterials* MINERMaterials::theInstance = 0;

const MINERMaterials* MINERMaterials::Instance() {
  if (!theInstance) theInstance = new MINERMaterials;
  return theInstance;
}

MINERMaterials::MINERMaterials()
{
  nistManager = G4NistManager::Instance();

  // material definitions go here:

  vac = new G4Material("Vacuum", 1., 1.01*g/mole, universe_mean_density,kStateGas,0.000017*kelvin,1.e-19*pascal);


  //Define Stainless Steel
  G4Element* C  = nistManager->FindOrBuildElement("C");
  G4Element* Si = nistManager->FindOrBuildElement("Si");
  G4Element* Cr = nistManager->FindOrBuildElement("Cr");
  G4Element* Mn = nistManager->FindOrBuildElement("Mn");
  G4Element* Fe = nistManager->FindOrBuildElement("Fe");
  G4Element* Ni = nistManager->FindOrBuildElement("Ni");
  StainlessSteel = new G4Material("StainlessSteel", 8.06*g/cm3, 6);
  StainlessSteel->AddElement(C, 0.001);
  StainlessSteel->AddElement(Si, 0.007);
  StainlessSteel->AddElement(Cr, 0.18);
  StainlessSteel->AddElement(Mn, 0.01);
  StainlessSteel->AddElement(Fe, 0.712);
  StainlessSteel->AddElement(Ni, 0.09);

  liquidNitrogen = new G4Material("LiquidNitrogen",7.,14.01*g/mole,0.808*g/cm3, kStateLiquid,77.*kelvin);
  air  = nistManager->FindOrBuildMaterial("G4_AIR");


  //G4_POLYETHYLENE  (C_2H_4)_N-Polyethylene  0.94 g/cm3,  
  //             1     0.143711
  //             6     0.856289
  boron   = nistManager->FindOrBuildMaterial("G4_B");
  copper   = nistManager->FindOrBuildMaterial("G4_Cu");
  aluminum = nistManager->FindOrBuildMaterial("G4_Al");

  graphiteMat = nistManager->FindOrBuildMaterial("G4_GRAPHITE");

  water = nistManager->FindOrBuildMaterial("G4_WATER");

  det_Ge = nistManager->FindOrBuildMaterial("G4_Ge");
  det_Si = nistManager->FindOrBuildMaterial("G4_Si");

  shieldLead = nistManager->FindOrBuildMaterial("G4_Pb");
  HDPE = nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
  //http://www.shieldwerx.com
  // 5% borated polyethylene = SWX203
  boratedpoly05 = new G4Material("BoratedPoly05",1.06*g/cm3, 2);
  boratedpoly05->AddMaterial(boron,0.05);
  boratedpoly05->AddMaterial(HDPE,0.95);

  // 30% borated polyethylene = SWX210
  boratedpoly30 = new G4Material("BoratedPoly30",1.19*g/cm3, 2);
  boratedpoly30->AddMaterial(boron,0.3);
  boratedpoly30->AddMaterial(HDPE,0.70);


  // stand in for boraflex
  Boraflex = new G4Material("Boraflex",1.64*g/cm3, 2);
  Boraflex->AddMaterial(boron,0.276);
  Boraflex->AddMaterial(HDPE,0.724);

  // HD Concrete in NSC
  G4Element* H1 = new G4Element("H1","H1", 1, 1.0078*g/mole);
  G4Element* O16 = new G4Element("O16","O16", 8, 15.9949*g/mole);
  G4Element* Mg24 = new G4Element("Mg24","MG24", 12, 23.985*g/mole);
  G4Element* Mg25 = new G4Element("Mg25","Mg25", 12, 24.9858*g/mole);
  G4Element* Mg26 = new G4Element("Mg26","Mg26", 12, 25.9826*g/mole);
  G4Element* Al27 = new G4Element("Al27","Al27", 13, 26.9815*g/mole);
  G4Element* Si28 = new G4Element("Si28","Si28", 14, 27.9769*g/mole);
  G4Element* Si29 = new G4Element("Si29","Si29", 14, 28.9765*g/mole);
  G4Element* Si30 = new G4Element("Si30","Si30", 14, 29.9738*g/mole);
  G4Element* S32 = new G4Element("S32","S32", 16, 31.9721*g/mole);
  G4Element* S33 = new G4Element("S33","S33", 16, 32.9745*g/mole);
  G4Element* S34 = new G4Element("S34","S34", 16, 33.9679*g/mole);
  G4Element* Ca40 = new G4Element("Ca40","Ca40", 20, 39.9626*g/mole);
  G4Element* Ca42 = new G4Element("Ca42","Ca42", 20, 41.9586*g/mole);
  G4Element* Ca43 = new G4Element("Ca43","Ca43", 20, 42.9588*g/mole);
  G4Element* Ca44 = new G4Element("Ca44","Ca44", 20, 43.9555*g/mole);
  G4Element* Ca46 = new G4Element("Ca46","Ca46", 20, 45.9537*g/mole);
  G4Element* Ca48 = new G4Element("Ca48","Ca48", 20, 47.9525*g/mole);
  G4Element* Ti46 = new G4Element("Ti46","Ti46", 22, 45.9526*g/mole);
  G4Element* Ti47 = new G4Element("Ti47","Ti47", 22, 46.9518*g/mole);
  G4Element* Ti48 = new G4Element("Ti48","Ti48", 22, 47.9479*g/mole);
  G4Element* Ti49 = new G4Element("Ti49","Ti49", 22, 48.9479*g/mole);
  G4Element* Ti50 = new G4Element("Ti50","Ti50", 22, 49.9448*g/mole);
  G4Element* V51 = new G4Element("V51","V51", 23, 50.9440*g/mole);
  G4Element* Cr50 = new G4Element("Cr50","Cr50", 24, 49.9460*g/mole);
  G4Element* Cr52 = new G4Element("Cr52","Cr52", 24, 51.9405*g/mole);
  G4Element* Cr53 = new G4Element("Cr53","Cr53", 24, 52.9407*g/mole);
  G4Element* Cr54 = new G4Element("Cr54","Cr54", 24, 53.9389*g/mole);
  G4Element* Mn55 = new G4Element("Mn55","Mn55", 25, 54.9380*g/mole);
  G4Element* Fe54 = new G4Element("Fe54","Fe54", 26, 53.9396*g/mole);
  G4Element* Fe56 = new G4Element("Fe56","Fe56", 26, 55.9308*g/mole);
  G4Element* Fe57 = new G4Element("Fe57","Fe57", 26, 56.9354*g/mole);
  G4Element* Fe58 = new G4Element("Fe58","Fe58", 26, 57.9333*g/mole);


  HDConcrete = new G4Material("HDConcrete",3.53*g/cm3, 33);
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
  NE213 = new G4Material("NE213",0.874*g/cm3 ,2);
  G4Element* elH  = new G4Element("Hydrogen","H" ,1., 1.01*g/mole);
  //might be able to replace with xylene?
  G4Element* elC  = new G4Element("Carbon"  ,"C" ,6., 12.01*g/mole);
  NE213->AddElement(elC,0.5479);
  NE213->AddElement(elH,0.4521);

  // Mu-metal
  mumetal = new G4Material("MuMetal", 8.75*g/cm3, 2);
  mumetal->AddElement(nistManager->FindOrBuildElement("Fe"), 0.19);
  mumetal->AddElement(nistManager->FindOrBuildElement("Ni"), 0.81);

  // pseudocumene
  psc = new G4Material("Pseudocumene", 0.876*g/cm3,2,kStateLiquid);
  psc->AddElement(nistManager->FindOrBuildElement("H"),12);
  psc->AddElement(nistManager->FindOrBuildElement("C"),9);

  // mineraloil
  minOil = new G4Material("Dodecane", 0.75*g/cm3,2,kStateLiquid);
  minOil->AddElement(nistManager->FindOrBuildElement("H"), 26);
  minOil->AddElement(nistManager->FindOrBuildElement("C"), 12);

  liqHe = new G4Material("lHe",0.125*g/cm3,1,kStateLiquid);
  liqHe->AddElement(nistManager->FindOrBuildElement("He"), 1.0);

  // Acrylic
  Acry = new G4Material("Acrylic", 1.18*g/cm3, 3);
  Acry->AddElement(nistManager->FindOrBuildElement("C"),5); 
  Acry->AddElement(nistManager->FindOrBuildElement("O"),2); 
  Acry->AddElement(nistManager->FindOrBuildElement("H"),8); 


}
