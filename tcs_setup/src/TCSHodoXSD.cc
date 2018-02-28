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

#include "TCSHodoXSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "TCSTrackInformation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSHodoXSD::TCSHodoXSD(const G4String& name,
				   const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(0), lastID(-1)
{
  collectionName.insert(hitsCollectionName);

  G4cout << "TCSHodoXSD::TCSHodoXSD: constructed" << G4endl;
  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSHodoXSD::~TCSHodoXSD() 
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHodoXSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
  = new TCSHodoHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

  //  G4cout << "TCSHodoXSD::Initialize: initialized" << G4endl;
  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TCSHodoXSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  //  G4int pid = step->GetTrack()->GetDefinition()->GetPDGEncoding();
  TCSTrackInformation* info =
    (TCSTrackInformation*)(step->GetTrack()->GetUserInformation());
  G4int pid = info->GetOriginalParticle()->GetPDGEncoding();

  G4ThreeVector pos = step->GetTrack()->GetPosition();

  // Particle entering hodoscope, save kinetic energy for flux calc-s.
  if (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==
      "hodoXWorld_PV" &&
      step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==
      "HXBar_PV" &&
      step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
    G4int chan =step->GetPostStepPoint()->GetTouchableHandle()->GetCopyNumber();
    G4double ekin = step->GetTrack()->GetKineticEnergy();
    TCSHodoHit* hit = new TCSHodoHit(chan, pid, ekin, pos, 1);
    fHitsCollection->insert( hit );
    ////    return true;
  }

  // Particle in the bars, save energy deposit for dose calc-s.
  //  if (step->GetTrack()->GetVolume()->GetName() == "HXBar_PV") {

  G4TouchableHandle touch = step->GetPreStepPoint()->GetTouchableHandle();
  if (touch->GetVolume(1)->GetName() == "HodoXAssembly_PV" &&
      touch->GetVolume()->GetName() == "HXBar_PV") {
    G4int chan = touch->GetCopyNumber();
    G4double edep = step->GetTotalEnergyDeposit();
    TCSHodoHit* hit = new TCSHodoHit(chan, pid, edep, pos, 0);
    fHitsCollection->insert( hit );

    //G4cout << "TCSHodoXSD::ProcessHits:" << G4endl;
    //G4cout << "  stepvol: " << step->GetTrack()->GetVolume()->GetName()
    //<< "  prestepvol: "
    //<< step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()
    //<< "  prestepvol1: "
    //<< step->GetPreStepPoint()->GetTouchableHandle()->GetVolume(1)->GetName()
    //	   << G4endl;
    //    getchar();

    //G4cout << "TCSHodoXSD::ProcessHits:" << G4endl;
    //G4cout << "  chan = "
    //<< step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber()
    //	 << "  chan1 = "
    //<< step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1)
    //	 << "  chan2 = "
    //<< step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2)
    //	 << G4endl;
    //    hit->Print();
    //    getchar();
    ////    return true;
  }

  //  G4cout << "TCSHodoXSD::ProcessHits:" << G4endl;
  //  G4cout << " PreStepPoint volume:" <<
  //    step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()
  //  	 << G4endl;
  //  G4cout << " PostStepPoint volume:" <<
  //    step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName()
  //  	 << G4endl;

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHodoXSD::EndOfEvent(G4HCofThisEvent*)
{
  //  if ( verboseLevel>1 ) { 
  //  G4int nofHits = fHitsCollection->entries();
  //   G4cout << G4endl
  //          << "-------->Hits Collection: in this event they are " << nofHits 
  //          << " hits in the tracker chambers: " << G4endl;
  //  for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->delete();
  //}

  //  G4cout << "TCSHodoXSD::EndOfEvent: end of event" << G4endl;
  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
