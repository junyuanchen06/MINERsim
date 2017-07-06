#include <sstream>
#include "RootIO.hh"
#include "BaseHit.hh"
#include "BaseTrack.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

#if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
#include "Cintex/Cintex.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static RootIO* instance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootIO::RootIO():fileName("hits.root")
{
  fMessenger = new RootIOMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootIO::~RootIO()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootIO::Setup()
{
  theFile = new TFile(fileName, "RECREATE");
  theFile->cd();

  hElapsedTime = new TH1D("hElapsedTime", "hElapsedTime", 1, 0, 1);

  theTree = new TTree("theTree", "theTree");
  sHits = new TClonesArray("BaseHit");
  theTree->Branch("sHits", &sHits, 6400, 0);
  hitC = 0;
  sTracks = new TClonesArray("BaseTrack");
  theTree->Branch("sTracks", &sTracks, 6400, 0);
  trackC = 0;

  theTree->Branch("hitC", &hitC, "hitC/I");
  theTree->Branch("trackC", &trackC, "trackC/I");
  theTree->Branch("event", &event, "event/I");
  theTree->Branch("eInc", &eInc, "eInc/F");

  hitC = 0;
  trackC = 0;
  event = 1;
  eInc = 0.;

  hists = new HistManager(theFile);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void RootIO::SetFileName(G4String file)
{
  fileName = file;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String RootIO::GetFileName()
{
  return fileName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootIO* RootIO::GetInstance()
{
  if (instance == 0 )
    instance = new RootIO();
  return instance;
}

void RootIO::SetIncomingE(G4double en)
{
  eInc = en;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootIO::AddHits(MinerHitsCollection * zipHits, G4int detID)
{
  for (G4int i = 0; i < zipHits->entries(); i++){
    BaseHit* hit = new ((*sHits)[hitC]) BaseHit;
    G4ThreeVector vec = (*zipHits)[i]->GetPos();
    G4ThreeVector mom = (*zipHits)[i]->GetMom();
    hit->SetPos(vec.x(), vec.y(), vec.z());
    hit->Setp3(mom.x(), mom.y(), mom.z());
    hit->SetEkin((*zipHits)[i]->GetParticleEnergy());
    hit->SetEdep((*zipHits)[i]->GetEdep());
    hit->SetPid((*zipHits)[i]->GetPDGID());
    hit->SetDetID(detID);
    hit->SetTime((*zipHits)[i]->GetTime());
    hit->SetWeight((*zipHits)[i]->GetWeight());
    hit->SetPreProcess((*zipHits)[i]->GetPreProcess());
    hit->SetPostProcess((*zipHits)[i]->GetPostProcess());
    hit->SetInc(0);
    hitC++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootIO::AddTrack(const G4Track*  trk)
{
  G4ThreeVector vertex = trk->GetVertexPosition ();
  G4double x = vertex.x(), y = vertex.y(), z = vertex.z();

  G4ThreeVector mom = trk->GetMomentum();
  G4double px = mom.x(), py = mom.y(), pz = mom.z();

  BaseTrack* tr = new ((*sTracks)[trackC]) BaseTrack;
  tr->Setp4(px, py, pz, trk->GetTotalEnergy());
  tr->SetPos(x, y, z);
  tr->SetPid(trk->GetDynamicParticle()->GetPDGcode());
  tr->SetTime(trk->GetGlobalTime());
  tr->SetInc(trk->GetParentID());
  trackC++;
}

void RootIO::FillMonitoring(MinerHitsCollection * zipHits, G4int detID)
{
  for (G4int i = 0; i < zipHits->entries(); i++) {
    G4ThreeVector vec = (*zipHits)[i]->GetPos();
    G4ThreeVector mom = (*zipHits)[i]->GetMom();

    std::ostringstream detNumber_ss;
    std::ostringstream pidNumber_ss;

    detNumber_ss << detID;
    pidNumber_ss << (*zipHits)[i]->GetPDGID();

    G4String detNumber = detNumber_ss.str();
    G4String pidNumber = pidNumber_ss.str();

    hists->fill1DHist((*zipHits)[i]->GetParticleEnergy(),
                      "Det" + detNumber + "_Ekin_PID" + pidNumber, "",
                      500, 0, 10,
                      (*zipHits)[i]->GetWeight(),
                      "Det"+detNumber + "_Monitoring");
    hists->fill2DHist(vec.y(), vec.z(),
                      "Det" + detNumber + "_pos_PID" + pidNumber, "",
                      120, -600, 600, 120, -1900, -700,
                      (*zipHits)[i]->GetWeight(),
                      "Det" + detNumber + "_Monitoring");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootIO::updateElapsedTime(G4double elapsedTime)
{
  hElapsedTime->SetBinContent(1, elapsedTime);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootIO::Write()
{
  event++;
  theTree->Fill();
  sHits->Clear();
  sTracks->Clear();
  hitC = 0;
  trackC = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootIO::Close()
{
  theFile->cd();
  theFile->Write();
  theFile->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
