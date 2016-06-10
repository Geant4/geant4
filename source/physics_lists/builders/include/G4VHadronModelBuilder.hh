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
// $Id: G4VHadronModelBuilder.hh 66892 2013-01-17 10:57:59Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:  G4VHadronModelBuilder
//
// Author: 28 June 2009 V.Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef G4VHadronModelBuilder_h
#define G4VHadronModelBuilder_h 1

#include "G4HadronicInteraction.hh"
#include "globals.hh"

class G4VHadronModelBuilder 
{
public: 

  G4VHadronModelBuilder(const G4String& name ="");

  virtual ~G4VHadronModelBuilder();

  G4HadronicInteraction* GetModel();

  inline const G4String& GetName() const;

protected:

  virtual G4HadronicInteraction* BuildModel() = 0;

private:

  // copy constructor and hide assignment operator
  G4VHadronModelBuilder(G4VHadronModelBuilder &);
  G4VHadronModelBuilder & operator=(const G4VHadronModelBuilder &right);

  G4HadronicInteraction* model;
  G4String name; 
};

inline const G4String& G4VHadronModelBuilder::GetName() const
{
  return name;
}

#endif

