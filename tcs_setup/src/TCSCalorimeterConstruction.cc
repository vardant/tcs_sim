#include "TCSCalorimeterConstruction.hh"
#include "NPSModuleConstruction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"

#include <iostream>

using namespace std;

TCSCalorimeterConstruction::TCSCalorimeterConstruction() {
  Construct();
}

TCSCalorimeterConstruction::~TCSCalorimeterConstruction() {;}

G4LogicalVolume* TCSCalorimeterConstruction::GetCalorimeter() {
  return fCalorimeter;
}

void TCSCalorimeterConstruction::Construct() {

  G4NistManager* man = G4NistManager::Instance();
  //  man->SetVerbose(1);

  NPSModuleConstruction* ModuleConstruction = new NPSModuleConstruction(man);

  double xmod = ModuleConstruction->GetSizeX();
  double ymod = ModuleConstruction->GetSizeY();
  double zmod = ModuleConstruction->GetSizeZ();

  cout << "TCSCalorimeterConstruction::GetCalorimeter: module sizes:" << endl;
  cout << "  xmod = " << xmod/cm << " cm" << endl;
  cout << "  ymod = " << ymod/cm << " cm" << endl;
  cout << "  zmod = " << zmod/cm << " cm" << endl;

  double xstep = xmod + fFrameThick;
  double ystep = ymod + fFrameThick;

  double xcol = xstep*fNCOL + fFrameThick;
  double ycol = ystep*fNCOL + fFrameThick;
  double zcol = zmod;
  cout << "TCSCalorimeterConstruction::GetCalorimeter: calo sizes:" << endl;
  cout << "  xcol = " << xcol/cm << " cm" << endl;
  cout << "  ycol = " << ycol/cm << " cm" << endl;
  cout << "  zcol = " << zcol/cm << " cm" << endl;
  getchar();

  G4Material* Air = man->FindOrBuildMaterial("G4_Air");

  G4Box* caloBox = new G4Box("caloBox", xcol/2, ycol/2, zcol/2);
  fCalorimeter = new G4LogicalVolume(caloBox, Air, "Calo", 0, 0, 0);

  int copy_number = 0;
  const double zpos = 0.;

  double ypos = -ycol/2 + ystep/2;
  for (int irow = 0; irow<fNROW; irow++) {

    double xpos = -xcol/2 + xstep/2;
    for (int icol = 0; icol<fNCOL; icol++) {
      G4LogicalVolume* Module_log = ModuleConstruction->GetModule();
      new G4PVPlacement(0,                                //no rotation
			G4ThreeVector(xpos,ypos,zpos),    //at position
			Module_log,                       //its logical volume
			"Module",                         //its name
			fCalorimeter,                     //its mother  volume
			false,                            //no boolean operation
			copy_number,                      //copy number
			0);                               //overlaps checking
      copy_number++;
      xpos += xstep;
    }

    ypos += ystep;
  }

  fCalorimeter->SetVisAttributes (G4VisAttributes::Invisible);
}
