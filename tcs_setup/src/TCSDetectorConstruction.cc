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

#include "TCSDetectorConstruction.hh"
#include "TCSTargetSD.hh"
#include "TCSCalorimeterSD.hh"
#include "TCSHodoXSD.hh"
#include "TCSHodoYSD.hh"
////#include "TCSTrackerXSD.hh"
////#include "TCSTrackerYSD.hh"
#include "G4SDManager.hh"

// **** Magnetic field ******
// New include files - used for magnetic field
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "SimpleField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4NystromRK4.hh"
#include "G4SimpleHeum.hh"
#include "G4SimpleRunge.hh"
//#include "G4ChordFinder.hh"
//#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4UserLimits.hh"

// Shapes...
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"
//#include "UExtrudedSolid.hh"
#include "G4SystemOfUnits.hh"

// Others...
#include "G4AutoDelete.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSDetectorConstruction::TCSDetectorConstruction()
  : G4VUserDetectorConstruction(), fScoringVolume(0), fField(0), fEquation(0),
    fStepper(0), fChordFinder(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSDetectorConstruction::~TCSDetectorConstruction()
{
  delete fField;
  delete fEquation;
  delete fStepper;
  delete fChordFinder;
  delete fScoringVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* TCSDetectorConstruction::Construct()
{

  // Setup Magnetic Field here!!!
  ConstructField();

  // Setup sensitive detector!
  ConstructSDandField();

  //always return the physical World
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSDetectorConstruction::ConstructField() 
{
  static G4TransportationManager* trMgr= 
    G4TransportationManager::GetTransportationManager(); 

  // A field object is held by a field manager
  // Find the global Field Manager
  G4FieldManager* globalFieldMgr= trMgr->GetFieldManager();

  static G4bool fieldIsInitialized = false;
  if(!fieldIsInitialized)    {

    fField = new SimpleField();

    //The ChordFinder is an helper class to track particles 
    //in magnetic fields, it sets the accuracy to be used.

    fEquation = new G4Mag_UsualEqRhs (fField);

    fStepper = new G4ClassicalRK4 (fEquation);    //Default choice.
    //    fStepper = new G4SimpleHeum (fEquation);    //300 ev/min
    //    fStepper = new G4ImplicitEuler (fEquation);    //300 ev/min
    //    fStepper = new G4ExplicitEuler (fEquation);    //<300 ev/min
    //    fStepper = new G4SimpleRunge (fEquation);    //300 ev/min

    //Mag. field
    //    fStepper = new G4HelixImplicitEuler (fEquation);    //does not work
    //    fStepper = new G4HelixExplicitEuler (fEquation);    //slow
    //    fStepper = new G4HelixSimpleRunge (fEquation);    //does not work
    //    fStepper = new G4NystromRK4 (fEquation);    //slow

    fChordFinder = new G4ChordFinder(fField,1e-4*m,fStepper);
    globalFieldMgr->SetChordFinder(fChordFinder);
    globalFieldMgr->SetDetectorField(fField);
    globalFieldMgr->GetChordFinder()->SetDeltaChord(1e-4*m);
    globalFieldMgr->SetDeltaIntersection(1e-4*m);
    globalFieldMgr->SetDeltaOneStep(1e-4*m);

    G4cout << "Magnetic field has been constructed " << 
      "in TCSDetectorConstruction::ConstructField()" << G4endl;
    fieldIsInitialized = true; 
  }
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void TCSDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  //Avoid double initialization.
  static G4ThreadLocal G4bool initialized = false;
  if ( ! initialized ) {

    // Calorimeter SD

    TCSCalorimeterSD* caloSD = new TCSCalorimeterSD("CalorimeterSD",
						   "CalorimeterHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(caloSD);

    // Register the messenger for deleting
    //  G4AutoDelete::Register(caloSD);
  
    //  SetSensitiveDetector("CalorimeterAssembly", tcsSD, true);
    SetSensitiveDetector("Block", caloSD, true);
    SetSensitiveDetector("caloWorld", caloSD, true);

    // Hodoscope SD

    TCSHodoXSD* hodoxSD = new TCSHodoXSD("HodoscopeXSD", "HodoXHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(hodoxSD);
    SetSensitiveDetector("HXBar", hodoxSD, true);
    SetSensitiveDetector("hodoXWorld", hodoxSD, true);

    TCSHodoYSD* hodoySD = new TCSHodoYSD("HodoscopeYSD", "HodoYHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(hodoySD);
    SetSensitiveDetector("HYBar", hodoySD, true);
    SetSensitiveDetector("hodoYWorld", hodoySD, true);

    // Tracker SD

    ////    TCSTrackerXSD* trackerxSD = new TCSTrackerXSD("TrackerXSD",
    ////					  "TrackerXHitsCollection");
    ////    G4SDManager::GetSDMpointer()->AddNewDetector(trackerxSD);
    ////    SetSensitiveDetector("TXBar", trackerxSD, true);
    ////    SetSensitiveDetector("trackerXWorld", trackerxSD, true);

    ////    TCSTrackerYSD* trackerySD = new TCSTrackerYSD("TrackerYSD",
    ////					  "TrackerYHitsCollection");
    ////    G4SDManager::GetSDMpointer()->AddNewDetector(trackerySD);
    ////    SetSensitiveDetector("TYBar", trackerySD, true);
    ////    SetSensitiveDetector("trackerYWorld", trackerySD, true);

    TCSTargetSD* targetSD = new TCSTargetSD("TargetSD",
					    "TargetHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(targetSD);
    SetSensitiveDetector("TargetAssembly", targetSD, true);

    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    //  G4ThreeVector fieldValue = G4ThreeVector();
    //fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    //fMagFieldMessenger->SetVerboseLevel(1);
  
    // Register the field messenger for deleting
    //G4AutoDelete::Register(fMagFieldMessenger);

    initialized=true;
  }

}
