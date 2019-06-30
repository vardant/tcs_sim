#include "TCSCalorimeterConstruction.hh"
#include "NPSModuleConstruction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
//#include "G4MultiUnion.hh"
#include "G4SubtractionSolid.hh"

#include <iostream>

using namespace std;

TCSCalorimeterConstruction::TCSCalorimeterConstruction() {
  //  Construct();  call separately
}

TCSCalorimeterConstruction::~TCSCalorimeterConstruction() {;}

G4LogicalVolume* TCSCalorimeterConstruction::GetCalorimeter() {
  return fCalorimeter;
}

void TCSCalorimeterConstruction::Construct() {

  G4NistManager* man = G4NistManager::Instance();
  //  man->SetVerbose(1);

  // Module construction and its sizes.

  NPSModuleConstruction* ModuleConstruction = new NPSModuleConstruction(man);

  double xmod = ModuleConstruction->GetSizeX();
  double ymod = ModuleConstruction->GetSizeY();
  double zmod = ModuleConstruction->GetSizeZ();

  cout << "TCSCalorimeterConstruction::GetCalorimeter: module sizes:" << endl;
  cout << "  xmod = " << xmod/cm << " cm" << endl;
  cout << "  ymod = " << ymod/cm << " cm" << endl;
  cout << "  zmod = " << zmod/cm << " cm" << endl;

  // Calorimeter construction.

  double xstep = xmod + fFrameThick;
  double ystep = ymod + fFrameThick;

  double xcol = xstep*fNCOL + fFrameThick;
  double ycol = ystep*fNCOL + fFrameThick;
  double zcol = zmod;
  cout << "TCSCalorimeterConstruction::GetCalorimeter: calo sizes:" << endl;
  cout << "  xcol = " << xcol/cm << " cm" << endl;
  cout << "  ycol = " << ycol/cm << " cm" << endl;
  cout << "  zcol = " << zcol/cm << " cm" << endl;
  //  getchar();

  G4Material* Air = man->FindOrBuildMaterial("G4_AIR");

  ///  G4Box* caloBox = new G4Box("caloBox", xcol/2, ycol/2, zcol/2);
  //Calo box little bit longer to allow for flux hit count in the CalorimeterSD.
  G4Box* caloBox = new G4Box("caloBox", xcol/2, ycol/2, zcol/2+0.1*mm);
  fCalorimeter = new G4LogicalVolume(caloBox, Air, "Calorimeter_LV", 0, 0, 0);

  // Positioning of modules in the calorimeter.

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
			"Module_phys",                    //its name
			fCalorimeter,                     //its mother  volume
			false,                            //no boolean operation
			copy_number,                      //copy number
			0);                               //overlaps checking

      copy_number++;
      xpos += xstep;
    }

    ypos += ystep;
  }

  // Support frames, made from carbon composite (carbon fiber + epoxy).

  G4Element* C  = man->FindOrBuildElement("C");
  G4Element* O  = man->FindOrBuildElement("O");
  G4Element* H  = man->FindOrBuildElement("H");

  //Carbon fiber of average density.
  G4Material* carbon_fiber = new G4Material("Carbon fiber", 1.84*g/cm3, 1);
  carbon_fiber->AddElement(C, 1);

  //Bisphenol A epoxy resin, with typical epoxy density.
  G4Material* epoxy = new G4Material("Bisphenol A", 1.25*g/cm3, 3);
  epoxy->AddElement(C, 17);
  epoxy->AddElement(O,  3);
  epoxy->AddElement(H, 19);

  G4Material* carbon_composite = new G4Material("Carbon composite",1.4*g/cm3,2);
  carbon_composite->AddMaterial(carbon_fiber, 50.*perCent);
  carbon_composite->AddMaterial(epoxy,        50.*perCent);

  /*
  G4Box* cellBox = new G4Box("cellBox", xmod/2, ymod/2, fFrameWidth/2+0.5*mm);
  G4MultiUnion* cellsSolid = new G4MultiUnion("cellsSolid");

  ypos = -ycol/2 + ystep/2;
  for (int irow = 0; irow<fNROW; irow++) {

    double xpos = -xcol/2 + xstep/2;
    for (int icol = 0; icol<fNCOL; icol++) {

      G4RotationMatrix rotm  = G4RotationMatrix();
      G4ThreeVector pos = G4ThreeVector(xpos, ypos, zpos);
      G4Transform3D tr  = G4Transform3D(rotm, pos);

      cellsSolid->AddNode(*cellBox, tr);

      xpos += xstep;
    }

    ypos += ystep;
  }

  cellsSolid->Voxelize();

  G4Box* frameBox = new G4Box("frameBox", xcol/2, ycol/2, fFrameWidth/2);

  G4SubtractionSolid* frameSolid =
    new G4SubtractionSolid("frameSolid", frameBox, cellsSolid);

  G4LogicalVolume* frameLog = new G4LogicalVolume(frameSolid, carbon_composite,
						  "frame", 0, 0, 0);

  new G4PVPlacement(0,                                          //no rotation
		    G4ThreeVector(0.,0.,-zcol/2+fFrameWidth/2), //at position
		    frameLog,                                   //logical volume
		    "front frame",                              //its name
		    fCalorimeter,                               //mother  volume
		    false,                                      //no bool. op.
		    0,                                          //copy number
		    0);                                         //overlaps check

  double zpos_frame = -zcol/2+ModuleConstruction->GetBlockSizeZ()-fFrameWidth/2;
  new G4PVPlacement(0,                                          //no rotation
		    G4ThreeVector(0.,0.,zpos_frame),            //at position
		    frameLog,                                   //logical volume
		    "back frame",                               //its name
		    fCalorimeter,                               //mother  volume
		    false,                                      //no bool. op.
		    0,                                          //copy number
		    0);                                         //overlaps check
  */

  G4Box* voidBox = new G4Box("voidBox", xmod/2, ymod/2, fFrameWidth/2+0.5*mm);
  G4Box* frameBox = new G4Box("frameBox", (xmod+fFrameThick)/2,
			                  (ymod+fFrameThick)/2,
			                  fFrameWidth/2);
  G4SubtractionSolid* frameSolid =
    new G4SubtractionSolid("frameSolid", frameBox, voidBox);

  G4LogicalVolume* frameLog = new G4LogicalVolume(frameSolid, carbon_composite,
						  "frame", 0, 0, 0);

  copy_number = 0;

  ypos = -ycol/2 + ystep/2;
  for (int irow = 0; irow<fNROW; irow++) {

    double xpos = -xcol/2 + xstep/2;
    for (int icol = 0; icol<fNCOL; icol++) {
  
      new G4PVPlacement(0,                                //no rotation
			G4ThreeVector(xpos,ypos,-zcol/2+fFrameWidth/2), //at pos
			frameLog,                         //its logical volume
			"front frame",                    //its name
			fCalorimeter,                     //its mother  volume
			false,                            //no boolean operation
			copy_number,                      //copy number
			0);                               //overlaps checking

      double zpos_frame =
	     -zcol/2+ModuleConstruction->GetBlockSizeZ()-fFrameWidth/2;
      new G4PVPlacement(0,                                //no rotation
			G4ThreeVector(xpos,ypos,zpos_frame), //at pos
			frameLog,                         //its logical volume
			"back frame",                     //its name
			fCalorimeter,                     //its mother  volume
			false,                            //no boolean operation
			copy_number,                      //copy number
			0);                               //overlaps checking

      copy_number++;
      xpos += xstep;
    }

    ypos += ystep;
  }

  //  fCalorimeter->SetVisAttributes (G4VisAttributes::Invisible);
}
