//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
///////////////////////////////////////////////////////////////////////////////
// File: CCalHall.cc
// Description: CCalHall Geometry factory class for the experimental hall
///////////////////////////////////////////////////////////////////////////////
#include "CCalHall.hh"

#include "g4std/fstream"
#include "CCalutils.hh"

#include "CCalHcal.hh"
#include "CCalEcal.hh"

//#define debug

CCalHall::CCalHall(const G4String &name): CCalDetector(name) {}

CCalHall::~CCalHall() {}

int CCalHall::readFile() {
  ///////////////////////////////////////////////////////////////
  //Let's open the file
  G4cout << " ==> Opening file " << File() << " to read elements..."
       << G4endl;

  G4std::ifstream is;
  bool ok = openGeomFile(is, pathName, File());
  if (!ok)
    return 0;

  // Find *DO HcalTB96
  findDO(is, G4String("HcalTB96"));

  // Calorimeter boundaries
  readName(is,genMaterial);
  is >> dy_2Hall >> dx_2Hall >> jump;
#ifdef debug
  G4cout << tab << "General material: " << genMaterial  << " Size " << dy_2Hall
       << ", " << dx_2Hall << G4endl;
#endif
  
  ///////////////////////////////////////////////////////////////
  // Close the file
  G4cout << " ==> Closing file " << File() << G4endl;
  is.close();

  return 1;

}

void CCalHall::constructDaughters() {
  CCalHcal* hcal = new CCalHcal("HadronCalorimeter");
  addDetector(hcal);
  
  CCalEcal* ecal = new CCalEcal("CrystalMatrixModule");
  addDetector(ecal);
}
