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
//
// Author: Ivana Hrivnacova, 11/09/2018  (ivana@ipno.in2p3.fr)
// ---------------------------------------------------------------------

#include "G4VScoreNtupleWriter.hh"

//_____________________________________________________________________________
G4VScoreNtupleWriter* G4VScoreNtupleWriter::fgMasterInstance         = nullptr;
G4ThreadLocal G4VScoreNtupleWriter* G4VScoreNtupleWriter::fgInstance = nullptr;

//_____________________________________________________________________________
G4VScoreNtupleWriter* G4VScoreNtupleWriter::Instance()
{
  // This function invokes creating the objects on workes,
  // The master instance should be created by the user
  // via the concrete class constructor

  G4bool isMaster = !G4Threading::IsWorkerThread();

  if((!isMaster) && (fgInstance == nullptr))
  {
    if(fgMasterInstance != nullptr)
    {
      fgInstance = fgMasterInstance->CreateInstance();
    }
  }

  return fgInstance;
}

//_____________________________________________________________________________
G4VScoreNtupleWriter::G4VScoreNtupleWriter()
{
  G4bool isMaster = !G4Threading::IsWorkerThread();

  if(isMaster && (fgMasterInstance != nullptr))
  {
    G4ExceptionDescription description;
    description << "      "
                << "G4VScoreNtupleWriter on master already exists."
                << "Cannot create another instance.";
    G4Exception("G4VScoreNtupleWriter::G4VScoreNtupleWriter()", "Analysis_F001",
                FatalException, description);
  }
  if(fgInstance != nullptr)
  {
    G4ExceptionDescription description;
    description << "      "
                << "G4VScoreNtupleWriter on worker already exists."
                << "Cannot create another instance.";
    G4Exception("G4VScoreNtupleWriter::G4VScoreNtupleWriter()", "Analysis_F001",
                FatalException, description);
  }
  if(isMaster)
    fgMasterInstance = this;
  fgInstance = this;
}

//_____________________________________________________________________________
G4VScoreNtupleWriter::~G4VScoreNtupleWriter() { fgInstance = nullptr; }
