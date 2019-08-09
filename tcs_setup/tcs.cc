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
//
// Neutral Particle Spectrometer (based on Gabi's Amagnet.cc).
// Roadmap:
// define geometry using gdml
// define magnetic field
// define interpolated mag. field
// define scoring meshes
// define ROOT histograms and trees

#include "TCSDetectorConstruction.hh"
//#include "TCSActionInitialization.hh"
#include "TCSPrimaryGeneratorAction.hh"
#include "TCSRunAction.hh"
#include "TCSEventAction.hh"
#include "TCSTrackingAction.hh"
#include "TCSStackingAction.hh"
#include "TCSSteppingAction.hh"
#include "TCSHistoManager.hh"

// attempt to use the GDML parser
#include <vector>
////#include "G4GDMLParser.hh"

////#ifdef G4MULTITHREADED
////#include "G4MTRunManager.hh"
////#else
#include "G4RunManager.hh"
////#endif

#include "G4LogicalVolumeStore.hh"
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
////#include "QGSP_BERT.hh"
#include "NPSAddOptics.hh"

#include "G4VModularPhysicsList.hh"
////#include "G4EmStandardPhysics.hh"
//#include "G4EmStandardPhysics_option4.hh"
//#include "G4EmLivermorePhysics.hh"

#include "G4StepLimiterPhysics.hh"
#include "G4ScoringManager.hh"
#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

/*
void print_aux(const G4GDMLAuxListType* auxInfoList, G4String prepend="|")
{  
  for(std::vector<G4GDMLAuxStructType>::const_iterator
      iaux = auxInfoList->begin(); iaux != auxInfoList->end(); iaux++ )
    {
      G4String str=iaux->type;
      G4String val=iaux->value;
      G4String unit=iaux->unit;
      
      G4cout << prepend << str << " : " << val  << " " << unit << G4endl;

      if (iaux->auxList) print_aux(iaux->auxList, prepend + "|");
    }
  return;
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{

  cout << "==> TCS setup: argc = " << argc << endl;
  for (int i=0; i<argc; i++) {
    cout << "==>  argv " << i << " = " << argv[i] << endl;
  }

#ifdef G4UI_USE

  // Detect interactive mode (if no arguments) and define UI session

  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
    cout << "==>  UI session defined." << endl;
  }

#endif

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4int myseed = 1234;
  G4Random::setTheSeed(myseed);
 
  // Construct the default run manager

  ////#ifdef G4MULTITHREADED
  ////  G4MTRunManager* runManager = new G4MTRunManager;
  ////#else
  G4RunManager* runManager = new G4RunManager;
  ////#endif

  // Setup scoring stuff
  //  G4ScoringManager* scoringManager = G4ScoringManager::GetScoringManager();

  // GDML parser.

  ////  G4GDMLParser parser;

  // Uncomment the following if wish to avoid names stripping
  // parser.SetStripFlag(false);

  ///  parser.Read("tcs_gdmls/tcs_setup.gdml");
  ///  runManager->SetUserInitialization(new TCSDetectorConstruction(
  ///					parser.GetWorldVolume()));
  ////  runManager->SetUserInitialization(new TCSDetectorConstruction(
  ////  				    parser.GetWorldVolume(),
  ////  				    parser.GetVolume("Block")));

  TCSDetectorConstruction* setup = new TCSDetectorConstruction;

  runManager->SetUserInitialization(setup);

  // Example how to retrieve Auxiliary Information

  //  G4Colour  white   ()              ;  // white
  G4Colour  white   (1.0, 1.0, 1.0) ;  // white
  G4Colour  grey    (0.5, 0.5, 0.5) ;  // grey
  G4Colour  black   (0.0, 0.0, 0.0) ;  // black
  G4Colour  red     (1.0, 0.0, 0.0) ;  // red
  G4Colour  green   (0.0, 1.0, 0.0) ;  // green
  G4Colour  blue    (0.0, 0.0, 1.0) ;  // blue
  G4Colour  cyan    (0.0, 1.0, 1.0) ;  // cyan
  G4Colour  magenta (1.0, 0.0, 1.0) ;  // magenta 
  G4Colour  yellow  (1.0, 1.0, 0.0) ;  // yellow
  G4Colour  copper  (0.5, 0.5, 0.0) ;
  G4Colour  orange  (1.0, 0.5, 0.0) ;
  /*
  const G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
  std::vector<G4LogicalVolume*>::const_iterator lvciter;

  for( lvciter = lvs->begin(); lvciter != lvs->end(); lvciter++ ) {

    G4GDMLAuxListType auxInfo = parser.GetVolumeAuxiliaryInformation(*lvciter);
    std::vector<G4GDMLAuxPairType>::const_iterator ipair = auxInfo.begin();

    for( ipair = auxInfo.begin(); ipair != auxInfo.end(); ipair++ ) {

      G4String str=ipair->type;
      G4String val=ipair->value;
      //      G4cout << " Auxiliary Information is found for Logical Volume :  "
      //	     << (*lvciter)->GetName() << G4endl;
      //      G4cout << " Name of Auxiliary type is     :  "
      //	     << str << G4endl;
      //      G4cout << " Associated Auxiliary value is :  "
      //	     << val << G4endl;

      if (str.compare("Color")==0) {

	G4VisAttributes *ttt =
	  new G4VisAttributes((*lvciter)->GetVisAttributes());

	// can't use switch on strings so we need this ugly thing!
	if (val.contains("trans")) {
	  ttt->SetVisibility(false);
	  (*lvciter)->SetVisAttributes(ttt);
	}
	if (val.contains("red")) {
	  (*lvciter)->SetVisAttributes(new G4VisAttributes(red));
	}
	if (val.contains("yellow")) {
	  (*lvciter)->SetVisAttributes(new G4VisAttributes(yellow));
	}
	if (val.contains("green")) {
	  (*lvciter)->SetVisAttributes(new G4VisAttributes(green));
	}
	if (val.contains("blue")) {
	  (*lvciter)->SetVisAttributes(new G4VisAttributes(blue));
	}
	if (val.contains("copper")) {
	  (*lvciter)->SetVisAttributes(new G4VisAttributes(copper));
	}
	if (val.contains("orange")) {
	  (*lvciter)->SetVisAttributes(new G4VisAttributes(orange));
	}
	if (val.contains("cyan")) {
	  (*lvciter)->SetVisAttributes(new G4VisAttributes(G4Colour::Cyan));
	}

	if (val.contains("wire")) {
	  ttt->SetVisibility(true);
	  ttt->SetForceWireframe(true);
	  (*lvciter)->SetVisAttributes(ttt);
	}

      }   //if Color

    }     //ipair
  }       //lvciter
  */
  /*
  const G4GDMLAuxMapType* auxmap = parser.GetAuxMap();
  G4cout << "Found " << auxmap->size()
	 << " volume(s) with auxiliary information."
	 << G4endl << G4endl;
  for(G4GDMLAuxMapType::const_iterator iter=auxmap->begin();
      iter!=auxmap->end(); iter++) 
    {
      G4cout << "Volume " << ((*iter).first)->GetName()
	     << " has the following list of auxiliary information: "
	     << G4endl << G4endl;
      for (G4GDMLAuxListType::const_iterator vit=(*iter).second.begin();
	   vit!=(*iter).second.end(); vit++)
	{
	  G4cout << "--> Type: " << (*vit).type
		 << " Value: " << (*vit).value << G4endl;
	  if ((*vit).type=="Color")
	    {
	      G4Colour tint(G4Colour::White());
	      if(G4Colour::GetColour((*vit).value, tint)) {
		((*iter).first)->SetVisAttributes(new G4VisAttributes(tint));
		G4cout << "    Assigned color " << (*vit).value
		       << " to volume " << ((*iter).first)->GetName()
		       <<  G4endl << G4endl;
	      }
	      else
		G4cout << "     Color " << (*vit).value << " not found."
		       << G4endl;
	      ;
	    }
	}
    }
  G4cout << G4endl;
   
  // now the 'global' auxiliary info
  G4cout << std::endl;
  G4cout << "Global auxiliary info:" << std::endl;
  G4cout << std::endl;

  print_aux(parser.GetAuxList());   

  G4cout << std::endl;
   
  // End of Auxiliary Information block
  */

  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  ////  G4VModularPhysicsList* physicsList = new QGSP_BERT;

  //  G4VModularPhysicsList* physicsList = new G4VModularPhysicsList;
  ////  physicsList->RegisterPhysics(new G4EmStandardPhysics());
  ////  physicsList->RegisterPhysics(new G4EmStandardPhysics_option4());
  //  physicsList->RegisterPhysics(new G4EmLivermorePhysics());

  physicsList->RegisterPhysics(new NPSAddOptics("Cerenkov & Scintillation") );
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);

  // User action initialization
  //  runManager->SetUserInitialization(new TCSActionInitialization());

  // Set a HistoManager.
  // Look for root filename in the input line, otherwise do the default.

  char const *kinFileName="tcs_gen.kin_data";
  char const *rootFileName="tcs_setup.root";
  if ( argc >= 3 )    kinFileName=argv[2];
  if ( argc >= 4 )    rootFileName=argv[3];

  TCSHistoManager*  histo = new TCSHistoManager(kinFileName,rootFileName);

  // Set user action classes.

  TCSPrimaryGeneratorAction* gen_action = new TCSPrimaryGeneratorAction();
  runManager->SetUserAction(gen_action);

  TCSRunAction* run_action = new TCSRunAction(histo);  
  runManager->SetUserAction(run_action);

  TCSEventAction* event_action = new TCSEventAction(histo);
  runManager->SetUserAction(event_action);

  TCSTrackingAction* tracking_action = new TCSTrackingAction();
  runManager->SetUserAction(tracking_action);

  TCSStackingAction* stacking_action = new TCSStackingAction();
  runManager->SetUserAction(stacking_action);

  TCSSteppingAction* stepping_action = new TCSSteppingAction(event_action);
  runManager->SetUserAction(stepping_action);

  // User Limits
  /*
  G4double maxStep1=0.5*cm, maxStep2 = 2.5*cm, maxTime = 5.*ms,
    minEkin = 1.*MeV, minEkinWorld=0.1*MeV;

  parser.GetVolume("BeamPipeVol")->SetUserLimits(new G4UserLimits(DBL_MAX,
								  DBL_MAX,
								  maxTime,
								  minEkin,
								  maxStep2));
  parser.GetVolume("World")->SetUserLimits(new G4UserLimits(DBL_MAX,
								  DBL_MAX,
								  maxTime,
								  minEkinWorld,
								  maxStep2));
  parser.GetVolume("CP1Vol")->SetUserLimits(new G4UserLimits(DBL_MAX,
								  DBL_MAX,
								  maxTime,
								  minEkin,
								  maxStep1));
  parser.GetVolume("CP2Vol")->SetUserLimits(new G4UserLimits(DBL_MAX,
								  DBL_MAX,
								  maxTime,
								  minEkin,
								  maxStep1));
  parser.GetVolume("CoilVol")->SetUserLimits(new G4UserLimits(DBL_MAX,
								  DBL_MAX,
								  maxTime,
								  minEkin,
								  maxStep2));
  parser.GetVolume("CoreVol")->SetUserLimits(new G4UserLimits(DBL_MAX,
								  DBL_MAX,
								  maxTime,
								  minEkin,
								  maxStep2));
  parser.GetVolume("InsertVol")->SetUserLimits(new G4UserLimits(DBL_MAX,
								  DBL_MAX,
								  maxTime,
								  minEkin,
								  maxStep2));
  parser.GetVolume("TopShieldVol")->SetUserLimits(new G4UserLimits(DBL_MAX,
								  DBL_MAX,
								  maxTime,
								  minEkin,
								  maxStep2));
  parser.GetVolume("SideShieldVol")->SetUserLimits(new G4UserLimits(DBL_MAX,
								  DBL_MAX,
								  maxTime,
								  minEkin,
								  maxStep2));
  */

#ifdef G4VIS_USE
  // Visualization manager construction
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
  cout << "==>  VisManager defined & initialized." << endl;
#endif
    
  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( argc != 1 ) {

    // execute an argument macro file if exist
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);

  }
  else {

#ifdef G4UI_USE

#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute init_vis.mac"); 
    cout << "==>  init_vis.mac executed." << endl;
#else
    UImanager->ApplyCommand("/control/execute init.mac"); 
    cout << "==>  init.mac executed." << endl;
#endif

    if (ui->IsGUI()) {
      ///         UImanager->ApplyCommand("/control/execute gui.mac");
    }     

    // start interactive session
    ui->SessionStart();
    cout << "==>  UI session started." << endl;
    delete ui;

#endif

  }   //argc!=1

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
