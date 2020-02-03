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
#include "G4HadronicException.hh"

#include <iostream>

G4HadronicException::G4HadronicException(G4String in, G4int at, G4String mess)
    : theMessage{mess}, theName{in}, theLine{at} {

  std::ostringstream os;
  Report(os);
  whatString = os.str();

  G4cout << whatString;

  if (std::getenv("DumpCoreOnHadronicException")) {
    G4Exception("G4HadronicException", "007", FatalException,
                "Fatal problem in above location");
  }
}

G4HadronicException::~G4HadronicException() throw()
{}

const char* G4HadronicException::what() const noexcept {
  return whatString.c_str();
}

void G4HadronicException::Report(std::ostream &aS) const {
  aS << "In " << theName << ", line " << theLine << ": " << std::endl;
  aS << "===> " << theMessage << std::endl;
}
