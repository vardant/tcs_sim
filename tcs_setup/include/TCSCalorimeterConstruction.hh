#ifndef TCSCalorimeterConstruction_h
#define TCSCalorimeterConstruction_h 1

#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

class G4LogicalVolume;
class G4Material;

class TCSCalorimeterConstruction {

public:
  TCSCalorimeterConstruction();
  ~TCSCalorimeterConstruction();
  G4LogicalVolume* GetCalorimeter();

private:

  const int fNCOL = 23;
  const int fNROW = 23;

  const double fFrameThick = 0.5*mm;
  const double fFrameWidth = 2.*cm;

  void Construct();

  G4LogicalVolume* fCalorimeter;
};

#endif
