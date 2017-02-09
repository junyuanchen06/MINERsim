#ifndef MINERMaterials_hh
#define MINERMaterials_hh 1

class G4Material;
class G4NistManager;

class MINERMaterials {

public:


  static const MINERMaterials* Instance();
  static const MINERMaterials* GetInstance() { return Instance(); }

  G4Material* GetVacuum() const { return vac; }
  G4Material* GetStainlessSteel() const { return StainlessSteel; }
  G4Material* GetLiN2() const { return liquidNitrogen; }
  G4Material* GetAir() const { return air; }
  G4Material* GetBoron() const { return boron; }
  G4Material* GetCopper() const { return copper; }
  G4Material* GetAluminum() const { return aluminum; }
  G4Material* GetGraphite() const { return graphiteMat; }
  G4Material* GetWater() const { return water; }
  G4Material* GetDetGe() const { return det_Ge; }
  G4Material* GetDetSi() const { return det_Si; }
  G4Material* GetShieldLead() const { return shieldLead; }
  G4Material* GetHDPoly() const { return HDPE; }
  G4Material* GetBoratedPoly05() const { return boratedpoly05; }
  G4Material* GetBoratedPoly30() const { return boratedpoly30; }
  G4Material* GetBoraflex() const { return Boraflex; }
  G4Material* GetHDConcrete() const { return HDConcrete; }
  G4Material* GetNE213() const { return NE213; }
  G4Material* GetMuMetal() const { return mumetal; }
  G4Material* GetPseudocumene() const { return psc; }
  G4Material* GetMineralOil() const { return minOil; }
  G4Material* GetLiHe() const { return liqHe; }
  G4Material* GetAcrylic() const { return Acry; }
  G4Material* GetNeopreneBlend() const { return NeopreneBlend;}
  G4Material* GetTin() const { return tin; }

private:
  MINERMaterials();
  ~MINERMaterials();

  static MINERMaterials *theInstance;

  G4NistManager *nistManager;

  G4Material *vac;
  G4Material *StainlessSteel;
  G4Material *liquidNitrogen;
  G4Material *air;
  G4Material *boron;
  G4Material *copper;
  G4Material *aluminum;
  G4Material *graphiteMat;
  G4Material *water;
  G4Material *det_Ge;
  G4Material *det_Si;
  G4Material *shieldLead;
  G4Material *HDPE;
  G4Material *boratedpoly05;
  G4Material *boratedpoly30;
  G4Material *Boraflex;
  G4Material *HDConcrete;
  G4Material *NE213;
  G4Material *mumetal;
  G4Material *psc;
  G4Material *minOil;
  G4Material *liqHe;
  G4Material *Acry;
  G4Material *NeopreneBlend;
  G4Material *tin;

};

#endif /* MINERMaterals_hh */
