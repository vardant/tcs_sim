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
#include "TCSCalorimeterConstruction.hh"
#include "TCSTargetSD.hh"
#include "TCSCalorimeterSD.hh"
#include "TCSHodoXSD.hh"
#include "TCSHodoYSD.hh"
#include "TCSTrackerSD.hh"
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

#include "G4Material.hh"
#include "G4MaterialTable.hh"

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

  G4NistManager* man = G4NistManager::Instance();
  //  man->SetVerbose(1);
  G4Material* Air = man->FindOrBuildMaterial("G4_AIR");

  // Create the hall
  //
  G4Box * WorldBox = new G4Box("WorldBox",fXWorld/2, fYWorld/2, fZWorld/2);
  G4LogicalVolume * WorldLV = new G4LogicalVolume(WorldBox, Air,"WorldLV");
  physWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), WorldLV, "World",
				0, false, 0);

  // Construct calorimeter.

  TCSCalorimeterConstruction CalorimeterConstruction;
  CalorimeterConstruction.Construct();
  G4LogicalVolume* Calorimeter_log = CalorimeterConstruction.GetCalorimeter();

  //  for (int quarter=0; quarter<4; quarter++)
  for (int quarter=0; quarter<1; quarter++)
    PositionCalorimeter(Calorimeter_log, quarter);


  // Read trackers from the gdml files.

  G4LogicalVolume* Tracker1_log = GetGDMLVolume(
		   "tcs_gdmls/pointer_referenced/tracker1_ref.gdml",
		   "TrackerAssembly0xe91030");
  //		   "tracker1World0xe915c0");
  for (int quarter=0; quarter<4; quarter++)
    PositionTracker(Tracker1_log, quarter, 1);

  G4LogicalVolume* Tracker2_log = GetGDMLVolume(
		   "tcs_gdmls/pointer_referenced/tracker2_ref.gdml",
		   "TrackerAssembly0x1a10030");
		   //		   "tracker2World0x1a105c0");
  for (int quarter=0; quarter<4; quarter++)
    PositionTracker(Tracker2_log, quarter, 2);

  G4LogicalVolume* Tracker3_log = GetGDMLVolume(
		   "tcs_gdmls/pointer_referenced/tracker3_ref.gdml",
		   "TrackerAssembly0x1760030");
		   //		   "tracker3World0x17605c0");
  for (int quarter=0; quarter<4; quarter++)
    PositionTracker(Tracker3_log, quarter, 3);

  //Scattering chamber.

  fParser.ReadModule("tcs_gdmls/scattering_chamber.gdml");
  G4LogicalVolume* chamber_log = fParser.GetVolume("ChamberAssembly");
  G4RotationMatrix*  cham_rot = new G4RotationMatrix();
  cham_rot->rotateY(90.*degree);
  new G4PVPlacement(cham_rot,
		    G4ThreeVector(0.,0.,0.),
		    chamber_log,
		    "ChamberAssembly_PV",
		    physWorld->GetLogicalVolume(), //its mother  volume
		    false,
		    0,
                    false);

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
    ////    SetSensitiveDetector("Block", caloSD, true);
    ////    SetSensitiveDetector("caloWorld", caloSD, true);

    SetSensitiveDetector("Block_log", caloSD, true);
    SetSensitiveDetector("Counter_log", caloSD, true);
    SetSensitiveDetector("Calorimeter_LV", caloSD, true);

    // Hodoscope SD
    /*
    TCSHodoXSD* hodoxSD = new TCSHodoXSD("HodoscopeXSD", "HodoXHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(hodoxSD);
    SetSensitiveDetector("HXBar", hodoxSD, true);
    SetSensitiveDetector("hodoXWorld", hodoxSD, true);

    TCSHodoYSD* hodoySD = new TCSHodoYSD("HodoscopeYSD", "HodoYHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(hodoySD);
    SetSensitiveDetector("HYBar", hodoySD, true);
    SetSensitiveDetector("hodoYWorld", hodoySD, true);
    */
    // Tracker SD

    TCSTrackerSD* trackerSD = new TCSTrackerSD("TrackerSD",
					       "TrackerHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(trackerSD);
    //    SetSensitiveDetector("Drift", trackerSD, true);
    SetSensitiveDetector("Drift0xe15800", trackerSD, true);   //tracker 1
    SetSensitiveDetector("Drift0x1994830", trackerSD, true);  //tracker 2
    SetSensitiveDetector("Drift0x16e4830", trackerSD, true);  //tracker 1
    //    SetSensitiveDetector("tracker1World0xe915c0", trackerSD, true);
    //    SetSensitiveDetector("tracker2World0x1a105c0", trackerSD, true);
    //    SetSensitiveDetector("tracker3World0x17605c0", trackerSD, true);

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

//==============================================================================

void TCSDetectorConstruction::PositionCalorimeter(
			    G4LogicalVolume* Calorimeter_log, int quarter) {

  // Positioning of the calorimeter quarters.

  double phi, theta_tilt, theta_pos;

  switch (quarter) {
  case 0:
    phi        = Calo.RotationAngle;
    theta_tilt = -Calo.TiltAngle;
    theta_pos  =  Calo.PositionAngle;
    break;
  case 1:
    phi        = -Calo.RotationAngle;
    theta_tilt = -Calo.TiltAngle;
    theta_pos  =  Calo.PositionAngle;
    break;
  case 2:
    phi        = -Calo.RotationAngle;
    theta_tilt =  Calo.TiltAngle;
    theta_pos  = -Calo.PositionAngle;
    break;
  case 3:
    phi        =  Calo.RotationAngle;
    theta_tilt =  Calo.TiltAngle;
    theta_pos  = -Calo.PositionAngle;
    break;
  default:
    phi        = 0.;
    theta_tilt = 0.;
    theta_pos  = 0.;
  }


  /*
  //This is how it should be.
 // u, v, w are the daughter axes, projected on the mother frame
  G4ThreeVector u = G4ThreeVector(cos(phi), 0.,-sin(phi));
  G4ThreeVector v = G4ThreeVector(-sin(theta_tilt)*sin(phi), cos(theta_tilt),
				  -sin(theta_tilt)*cos(phi));
  G4ThreeVector w = G4ThreeVector(sin(phi), 0., cos(phi));
  G4RotationMatrix rotm = G4RotationMatrix(u, v, w);
  G4cout << "Direct rotation matrix : ";
  rotm.print(G4cout);     
  */

  //This is consistent with gdml coding.
  G4RotationMatrix rotm;
  rotm.rotateY(phi);
  rotm.rotateX(theta_tilt);
  rotm.rotateZ(0.);

  G4ThreeVector position=G4ThreeVector(sin(phi)*cos(theta_pos), sin(theta_pos),
				       cos(phi)*cos(theta_pos));
  position *= Calo.Distance;

  G4Transform3D transform = G4Transform3D(rotm, position);

  new G4PVPlacement(transform,                     //position, rotation        
                    Calorimeter_log,               //logical volume
                    "Calorimeter_PV",              //name
		    physWorld->GetLogicalVolume(), //its mother  volume
                    false,                         //no boolean operation
                    quarter);                      //copy number
}
//==============================================================================

void TCSDetectorConstruction::PositionTracker(G4LogicalVolume* Tracker_log,
					      int quarter, int layer) {

  // Positioning of the tracker quarters.

  double phi, theta_tilt, theta_pos;

  switch (quarter) {
  case 0:
    phi        =  Tracker.RotationAngle;
    theta_tilt = -Tracker.TiltAngle;
    theta_pos  =  Tracker.PositionAngle;
    break;
  case 1:
    phi        = -Tracker.RotationAngle;
    theta_tilt = -Tracker.TiltAngle;
    theta_pos  =  Tracker.PositionAngle;
    break;
  case 2:
    phi        = -Tracker.RotationAngle;
    theta_tilt =  Tracker.TiltAngle;
    theta_pos  = -Tracker.PositionAngle;
    break;
  case 3:
    phi        =  Tracker.RotationAngle;
    theta_tilt =  Tracker.TiltAngle;
    theta_pos  = -Tracker.PositionAngle;
    break;
  default:
    phi        = 0.;
    theta_tilt = 0.;
    theta_pos  = 0.;
  }

  /*
  //This is how it should be.
 // u, v, w are the daughter axes, projected on the mother frame
  G4ThreeVector u = G4ThreeVector(cos(phi), 0.,-sin(phi));
  G4ThreeVector v = G4ThreeVector(-sin(theta_tilt)*sin(phi), cos(theta_tilt),
				  -sin(theta_tilt)*cos(phi));
  G4ThreeVector w = G4ThreeVector(sin(phi), 0., cos(phi));
  G4RotationMatrix rotm = G4RotationMatrix(u, v, w);
  G4cout << "Direct rotation matrix : ";
  rotm.print(G4cout);     
  */

  //This is consistent with gdml coding.
  G4RotationMatrix rotm;
  rotm.rotateY(phi);
  rotm.rotateX(theta_tilt);
  rotm.rotateZ(0.);

  G4ThreeVector position=G4ThreeVector(sin(phi)*cos(theta_pos), sin(theta_pos),
				       cos(phi)*cos(theta_pos));
  position *= Tracker.Distance[layer-1];

  G4Transform3D transform = G4Transform3D(rotm, position);

  new G4PVPlacement(transform,                     //position, rotation        
                    Tracker_log,                   //logical volume
                    "Tracker"+to_string(layer),    //name
		    physWorld->GetLogicalVolume(), //its mother  volume
                    false,                         //no boolean operation
                    quarter);                      //copy number
}

//==============================================================================

G4LogicalVolume* TCSDetectorConstruction::GetGDMLVolume(const string file_name,
							const string vol_name) {
  fParser.ReadModule(file_name);
  return fParser.GetVolume(vol_name);
}
