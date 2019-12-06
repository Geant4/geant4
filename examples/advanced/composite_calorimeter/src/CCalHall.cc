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
///////////////////////////////////////////////////////////////////////////////
// File: CCalHall.cc
// Description: CCalHall Geometry factory class for the experimental hall
///////////////////////////////////////////////////////////////////////////////
#include "CCalHall.hh"

#include <fstream>
#include "CCalutils.hh"

#include "CCalHcal.hh"
#include "CCalEcal.hh"

//#define debug

CCalHall::CCalHall(const G4String &name): CCalDetector(name) {}

CCalHall::~CCalHall() {}

G4int CCalHall::readFile() {
  ///////////////////////////////////////////////////////////////
  //Let's open the file
  G4cout << " ==> Opening file " << File() << " to read elements..."
       << G4endl;

  std::ifstream is;
  G4bool ok = openGeomFile(is, pathName, File());
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
