///////////////////////////////////////////////////////////////////////////////
// File: HcalTB96.cc
// Date: 08/00 S.Banerjee
// Modifications:
///////////////////////////////////////////////////////////////////////////////
#include "HcalTB96.hh"

#include <fstream.h>
#include "utils.hh"

#include "HcalTB96HCal.hh"
#include "CrystalMatrix.hh"

//#define debug

HcalTB96::HcalTB96(const G4String &name): CMSDetector(name) {}

HcalTB96::~HcalTB96() {}

int HcalTB96::readFile() {
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

void HcalTB96::constructDaughters() {
  HcalTB96HCal* hcal = new HcalTB96HCal("HadronCalorimeter");
  addDetector(hcal);
  
  CrystalMatrix* xtalmod = new CrystalMatrix("CrystalMatrixModule");
  addDetector(xtalmod);
}
