#ifndef NPSModuleConstruction_h
#define NPSModuleConstruction_h 1

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;

class NPSModuleConstruction {

public:
  NPSModuleConstruction(G4NistManager* man);
  ~NPSModuleConstruction();

public:

  //  const G4VPhysicalVolume* Get_block() {return block_phys;};
  G4LogicalVolume* GetModule();
  G4double GetSizeX();
  G4double GetSizeY();
  G4double GetSizeZ();
  G4double GetBlockSizeZ();

private:

  G4double tedlar_x;
  G4double tedlar_y;
  G4double tedlar_z;

  G4double mylar_x;
  G4double mylar_y;
  G4double mylar_z;

  G4double tedlar_thick;
  G4double mylar_thick;
  G4double glue_thick;
  G4double air_gap;

  G4int refFlag;
  G4String refName;
  G4int refNumData;
  G4double* refWL;
  G4double* refReIndex;
  G4double* refImIndex;
  G4double subRefrIndex;
  G4double* refRefl;
  bool fFrontCoverFlag;
  
  G4double block_x;
  G4double block_y;
  G4double block_z;

  G4double PMT_diameter;
  G4double PMTWin_thick;
  G4double Cathode_diam;
  G4double Cathode_thick;

  G4double counter_x;
  G4double counter_y;
  G4double counter_z;

  ///  G4double expHall_x;
  ///  G4double expHall_y;
  ///  G4double expHall_z;

  G4LogicalVolume* block_log;
  G4LogicalVolume* mylar_log;
  G4LogicalVolume* PMTWin_right_log;
  G4LogicalVolume* PMTHouse_log;
  G4LogicalVolume* Cathode_log;
  G4LogicalVolume* glue_log;
  G4LogicalVolume* tedlar_log;
  G4LogicalVolume* counter_log;
  ///  G4VPhysicalVolume* counter_phys;
  ///  G4LogicalVolume* expHall_log;

  void Construct(G4NistManager* man);
};

#endif /*NPSModuleConstruction_h*/
