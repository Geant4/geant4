///////////////////////////////////////////////////////////////////////////////
// File: CCalHall.cc
// Description: CCalHall Geometry factory class for the experimental hall
///////////////////////////////////////////////////////////////////////////////
#include "CCalHall.hh"

#include <fstream.h>
#include "CCalutils.hh"

#include "CCalHcal.hh"
#include "CCalEcal.hh"

//#define debug

CCalHall::CCalHall(const G4String &name): CCalDetector(name) {}

CCalHall::~CCalHall() {}

int CCalHall::readFile() {
  ///////////////////////////////////////////////////////////////
  //Let's open the file
  cout << " ==> Opening file " << File() << " to read elements..."
       << endl;

  ifstream is;
  bool ok = openGeomFile(is, pathName, File());
  if (!ok)
    return 0;

  // Find *DO HcalTB96
  findDO(is, G4String("HcalTB96"));

  // Calorimeter boundaries
  readName(is,genMaterial);
  is >> dy_2Hall >> dx_2Hall >> jump;
#ifdef debug
  cout << tab << "General material: " << genMaterial  << " Size " << dy_2Hall
       << ", " << dx_2Hall << endl;
#endif
  
  ///////////////////////////////////////////////////////////////
  // Close the file
  cout << " ==> Closing file " << File() << endl;
  is.close();

  return 1;

}

void CCalHall::constructDaughters() {
  CCalHcal* hcal = new CCalHcal("HadronCalorimeter");
  addDetector(hcal);
  
  CCalEcal* ecal = new CCalEcal("CrystalMatrixModule");
  addDetector(ecal);
}
