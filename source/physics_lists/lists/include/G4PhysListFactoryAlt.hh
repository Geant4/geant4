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
// $Id: G4PhysListFactoryAlt.hh 67015 2013-01-29 17:10:02Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:  G4PhysListFactoryAlt
//
// Author:  R. Hatcher  2014-10-15
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef G4PhysListFactoryAlt_h
#define G4PhysListFactoryAlt_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

namespace g4alt {

class G4PhysListFactory
{
public:

  G4PhysListFactory();

  ~G4PhysListFactory();

  G4VModularPhysicsList* GetReferencePhysList(const G4String&);
  // instantiate PhysList by name

  G4VModularPhysicsList* ReferencePhysList();
  // instantiate PhysList by environment variable "PHYSLIST"

  void SetDefaultReferencePhysList(const G4String& name="");
  // set a prefered list in case where $PHYSLIST isn't defined
  // if not set (or called with "") this falls back to system default

  G4bool IsReferencePhysList(const G4String&);
  // check if the name is in the list of PhysLists names

  const std::vector<G4String>& AvailablePhysLists() const;
  // list of avalable base Phys Lists

  const std::vector<G4String>& AvailablePhysListsEM() const;
  // list of avalable EM options

  void PrintAvailablePhysLists() const;
  // print what users can select

  void  SetVerbose(G4int val);
  G4int GetVerbose() const;

  void  SetUnknownFatal(G4int val);
  G4int GetUnknownFatal() const;
  // throw an exception if requested list is unsatisfiable?

private:

  // no data
};

}  // end-of-namespace 'g4alt'

#endif



