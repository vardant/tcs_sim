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

  double xcal = xstep*fNCOL + fFrameThick;
  double ycal = ystep*fNROW + fFrameThick;
  double zcal = zmod;
  cout << "TCSCalorimeterConstruction::GetCalorimeter: calo matrix sizes:"
       << endl;
  cout << "  xcal = " << xcal/cm << " cm" << endl;
  cout << "  ycal = " << ycal/cm << " cm" << endl;
  cout << "  zcal = " << zcal/cm << " cm" << endl;

  //Calorimeter's case sizes.
  double xcase = xcal + 1.*cm;
  double ycase = ycal + 1.*cm;
  double zcase = zcal + 1.*cm;
  cout << "                              Calorimeter's case sizes:" << endl;
  cout << "  xcase = " << xcase/cm << " cm" << endl;
  cout << "  ycase = " << ycase/cm << " cm" << endl;
  cout << "  zcase = " << zcase/cm << " cm" << endl;

  //  getchar();

  G4Material* Air = man->FindOrBuildMaterial("G4_AIR");

  G4Box* caloBox = new G4Box("caloBox", xcase/2, ycase/2, zcase/2);
  fCalorimeter = new G4LogicalVolume(caloBox, Air, "Calorimeter_LV", 0, 0, 0);

  // Positioning of modules in the calorimeter. Numbering of rows from bottom
  // top, columns from left to right.

  int copy_number = 0;
  const double zpos = 0.;

  double ypos = -ycal/2 + ystep/2;
  for (int irow = 0; irow<fNROW; irow++) {

    ///    double xpos = -xcal/2 + xstep/2;
    double xpos = xcal/2 - xstep/2;
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
      ///      xpos += xstep;
      xpos -= xstep;
    }

    ypos += ystep;
  }

  // Support frames, made from carbon composite (carbon fiber + epoxy,
  // front frame) and mu-metal (rear frame).

  G4Element* C  = man->FindOrBuildElement("C");
  G4Element* O  = man->FindOrBuildElement("O");
  G4Element* H  = man->FindOrBuildElement("H");

  //Carbon fiber of average density.
  G4Material* carbon_fiber = new G4Material("Carbon fiber", 1.84*g/cm3, 1);
  carbon_fiber->AddElement(C, 1);

  G4Element* Ni  = man->FindOrBuildElement("Ni");
  G4Element* Mo  = man->FindOrBuildElement("Mo");
  G4Element* Fe  = man->FindOrBuildElement("Fe");

  //mu-metal (ASTM A753 Alloy 4) according to mu-metal.com.
  G4Material* mu_metal = new G4Material("mu-metal", 8.7*g/cm3, 3);
  mu_metal->AddElement(Ni, 80.*perCent);
  mu_metal->AddElement(Mo,  5.*perCent);
  mu_metal->AddElement(Fe, 15.*perCent);

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

  ypos = -ycal/2 + ystep/2;
  for (int irow = 0; irow<fNROW; irow++) {

    double xpos = -xcal/2 + xstep/2;
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

  G4Box* frameBox = new G4Box("frameBox", xcal/2, ycal/2, fFrameWidth/2);

  G4SubtractionSolid* frameSolid =
    new G4SubtractionSolid("frameSolid", frameBox, cellsSolid);

  G4LogicalVolume* frameLog = new G4LogicalVolume(frameSolid, carbon_composite,
						  "frame", 0, 0, 0);

  new G4PVPlacement(0,                                          //no rotation
		    G4ThreeVector(0.,0.,-zcal/2+fFrameWidth/2), //at position
		    frameLog,                                   //logical volume
		    "front frame",                              //its name
		    fCalorimeter,                               //mother  volume
		    false,                                      //no bool. op.
		    0,                                          //copy number
		    0);                                         //overlaps check

  double zpos_frame = -zcal/2+ModuleConstruction->GetBlockSizeZ()-fFrameWidth/2;
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

  G4LogicalVolume* FrontFrameLog = new G4LogicalVolume(frameSolid,
				   carbon_composite, "front frame", 0, 0, 0);

  G4LogicalVolume* RearFrameLog = new G4LogicalVolume(frameSolid,
				  mu_metal, "rear frame", 0, 0, 0);

  copy_number = 0;

  ypos = -ycal/2 + ystep/2;
  for (int irow = 0; irow<fNROW; irow++) {

    ///    double xpos = -xcal/2 + xstep/2;
    double xpos = xcal/2 - xstep/2;
    for (int icol = 0; icol<fNCOL; icol++) {
  
      new G4PVPlacement(0,                                //no rotation
			G4ThreeVector(xpos,ypos,-zcal/2+fFrameWidth/2), //at pos
			FrontFrameLog,                    //its logical volume
			"front frame",                    //its name
			fCalorimeter,                     //its mother  volume
			false,                            //no boolean operation
			copy_number,                      //copy number
			0);                               //overlaps checking

      double zpos_frame =
	     -zcal/2+ModuleConstruction->GetBlockSizeZ()-fFrameWidth/2;
      new G4PVPlacement(0,                                //no rotation
			G4ThreeVector(xpos,ypos,zpos_frame), //at pos
			RearFrameLog,                     //its logical volume
			"rear frame",                     //its name
			fCalorimeter,                     //its mother  volume
			false,                            //no boolean operation
			copy_number,                      //copy number
			0);                               //overlaps checking

      copy_number++;
      ///      xpos += xstep;
      xpos -= xstep;
    }

    ypos += ystep;
  }

  // Case.

  G4Material* Al = man->FindOrBuildMaterial("G4_Al");

  G4Box* caseBox = new G4Box("caloBox", xcase/2, ycase/2, zcase/2);

  const double caseThick = 1.*mm;

  G4Box* caseVoid = new G4Box("caloBox", xcase/2-caseThick, ycase/2-caseThick,
			      zcase/2-caseThick);
  G4SubtractionSolid* caseSolid =
    new G4SubtractionSolid("caseSolid", caseBox, caseVoid);

  G4LogicalVolume* caseLog = new G4LogicalVolume(caseSolid, Al, "case_log",
						 0, 0, 0);

  new G4PVPlacement(0,                                //no rotation
		    G4ThreeVector(),                  //at pos
		    caseLog,                          //its logical volume
		    "CaloCase phys",                  //its name
		    fCalorimeter,                     //its mother  volume
		    false,                            //no boolean operation
		    0,                                //copy number
		    0);                               //overlaps checking

  //  fCalorimeter->SetVisAttributes (G4VisAttributes::Invisible);
}
