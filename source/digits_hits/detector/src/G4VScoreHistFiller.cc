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
// Class G4VScoreHistFiller
//
// Class description:
//  Abstract base class that defines functions to fill histogram and
//  profile plot. This class avoids the direct dependency to G4Analysis
//  category. Template class G4TScoreHistFiller should be used by the
//  user to specify the type of the file format.
//
//  Created: M. Asai (Sept. 2020)
//
//

#include "G4VScoreHistFiller.hh"

G4VScoreHistFiller* G4VScoreHistFiller::fgMasterInstance = nullptr;
G4ThreadLocal G4VScoreHistFiller* G4VScoreHistFiller::fgInstance = nullptr;

G4VScoreHistFiller* G4VScoreHistFiller::Instance()
{
  // This function invokes creating the objects on workes,
  // The master instance should be created by the user
  // via the concrete class constructor

  G4bool isMaster = ! G4Threading::IsWorkerThread();

  if ((! isMaster) && (fgInstance == nullptr)) {
    if (fgMasterInstance != nullptr) {
      fgInstance = fgMasterInstance->CreateInstance();
    }
  }

  return fgInstance;
}

G4VScoreHistFiller::G4VScoreHistFiller()
{
  G4bool isMaster = ! G4Threading::IsWorkerThread();

  if (isMaster && (fgMasterInstance != nullptr)) {
    G4ExceptionDescription description;
    description << "      "
                << "G4VScoreHistFiller on master already exists."
                << "Cannot create another instance.";
    G4Exception(
      "G4VScoreHistFiller::G4VScoreHistFiller()", "Analysis_F001", FatalException, description);
  }
  if (fgInstance != nullptr) {
    G4ExceptionDescription description;
    description << "      "
                << "G4VScoreHistFiller on worker already exists."
                << "Cannot create another instance.";
    G4Exception(
      "G4VScoreHistFiller::G4VScoreHistFiller()", "Analysis_F001", FatalException, description);
  }
  if (isMaster) fgMasterInstance = this;
  fgInstance = this;
}

G4VScoreHistFiller::~G4VScoreHistFiller() { fgInstance = nullptr; }
