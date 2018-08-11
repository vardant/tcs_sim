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

#include "TCSTrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "TCSTrackInformation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSTrackerSD::TCSTrackerSD(const G4String& name,
			   const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(0), lastID(-1)
{
  collectionName.insert(hitsCollectionName);

  //  G4cout << "TCSTrackerXSD::TCSTrackerXSD: constructed" << G4endl;
  //  G4cout << "   name = " << name << G4endl;
  //  G4cout << "   hitsCollectionName = " << hitsCollectionName << G4endl;
  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSTrackerSD::~TCSTrackerSD() 
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSTrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
  = new TCSTrackerHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

  //  G4cout << "TCSTrackerSD::Initialize: initialized" << G4endl;
  // G4cout << "   SensitiveDetectorName = " << SensitiveDetectorName << G4endl;
  //  G4cout << "   collectionName = " << collectionName[0] << G4endl;
  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TCSTrackerSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  TCSTrackInformation* info =
    (TCSTrackInformation*)(step->GetTrack()->GetUserInformation());

  // Particle entering tracker.
  if (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==
      "trackerWorld_PV" &&
      step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==
      "Insulator_PV" &&
      step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {

    //G4int chan =step->GetPostStepPoint()->GetTouchableHandle()->GetCopyNumber();
    G4int pid = step->GetTrack()->GetDefinition()->GetPDGEncoding();
    G4ThreeVector mom = step->GetTrack()->GetMomentum();
    G4double P = sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
    G4ThreeVector pos = step->GetTrack()->GetPosition();
    TCSTrackerHit* hit = new TCSTrackerHit(pid, P, pos);
    fHitsCollection->insert( hit );
    ////    return true;
  }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSTrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  //  if ( verboseLevel>1 ) { 
  //  G4int nofHits = fHitsCollection->entries();
  //   G4cout << G4endl
  //          << "-------->Hits Collection: in this event they are " << nofHits 
  //          << " hits in the tracker chambers: " << G4endl;
  //  for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->delete();
  //}

  //  G4cout << "TCSTrackerXSD::EndOfEvent: end of event" << G4endl;
  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
