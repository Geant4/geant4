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
// 18 April 2024   V.Ivanchenko
//
// This is a singleton class to store shared EM physics data.
// The access to data is possible only by data name.
// Once created data pointers does not changed or deleted until the
// end of the job, this class is responsible for destruction.
//

#ifndef G4EmDataRegistry_h
#define G4EmDataRegistry_h 1

#include <vector>
#include "G4EmDataHandler.hh"

class G4EmDataRegistry
{
 public:

  static G4EmDataRegistry* Instance();

  ~G4EmDataRegistry();

  // create or access data handler
  // nTable - number of tables in a new handler
  G4EmDataHandler* GetHandlerByName(const G4String&, std::size_t nTable);

  // register a new handler
  void Register(G4EmDataHandler*);

  // register a new handler
  void DeRegister(G4EmDataHandler*);
  
  // hide assignment operator 
  G4EmDataRegistry& operator=(const G4EmDataRegistry &right) = delete;
  G4EmDataRegistry(const G4EmDataRegistry&) = delete;

 private:

  G4EmDataRegistry();

  G4EmDataHandler* EmDataHandler(const G4String&);

  static G4EmDataRegistry* instance;

  std::vector<G4EmDataHandler*> fDataHandlers;
};

#endif
