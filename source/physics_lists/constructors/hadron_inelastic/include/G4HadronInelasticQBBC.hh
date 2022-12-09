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
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronInelasticQBBC
//
// Author: 31 August 2009 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4HadronInelasticQBBC_h
#define G4HadronInelasticQBBC_h 1

#include "globals.hh"
#include "G4VHadronPhysics.hh"

class G4HadronInelasticQBBC : public G4VHadronPhysics
{
public: 

  // constructor
  G4HadronInelasticQBBC(G4int ver = 1);

  G4HadronInelasticQBBC(const G4String& name,
			G4int ver = 1, G4bool ftf = false, G4bool bert = false,
			G4bool chips = false, G4bool hp = false, 
			G4bool glauber = false);

  virtual ~G4HadronInelasticQBBC() = default;

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type 
  void ConstructProcess() override;

  // copy constructor and hide assignment operator
  G4HadronInelasticQBBC(G4HadronInelasticQBBC &) = delete;
  G4HadronInelasticQBBC & operator=(const G4HadronInelasticQBBC &right) = delete;
};

#endif
