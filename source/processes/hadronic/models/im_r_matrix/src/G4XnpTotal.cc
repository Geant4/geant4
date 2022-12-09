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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//      For information related to this code contact:
//
//      File name:     G4XnpTotal
//
//      Author:        
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4XnpTotal.hh"
#include "G4XnpTotalLowE.hh"
#include "G4XPDGTotal.hh"
#include "G4VCrossSectionSource.hh"
#include "G4KineticTrack.hh"
#include "G4ParticleDefinition.hh"
#include "G4CrossSectionVector.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

G4XnpTotal::G4XnpTotal()
{ 
  components = new G4CrossSectionVector;

  G4VCrossSectionSource* xnpTotalLowE = new G4XnpTotalLowE;
  components->push_back(xnpTotalLowE);

  G4VCrossSectionSource* xnpTotalHighE = new G4XPDGTotal;
  components->push_back(xnpTotalHighE);
}


G4XnpTotal::~G4XnpTotal()
{ 
  if (components != nullptr) 
    {
      std::size_t nComponents = GetComponents()->size();
      for (std::size_t i=0; i<nComponents; ++i)
	{
          G4CrossSectionSourcePtr componentPtr = (*components)[i];
	  G4VCrossSectionSource* component = componentPtr();
	  delete component;
	  component = nullptr;
	  componentPtr = nullptr;
	}
    }
  delete components;
  components = nullptr;
}


G4bool G4XnpTotal::operator==(const G4XnpTotal &right) const
{
  return (this == (G4XnpTotal*) &right);
}


G4bool G4XnpTotal::operator!=(const G4XnpTotal &right) const
{
  return (this != (G4XnpTotal*) &right);
}


G4String G4XnpTotal::Name() const
{
  G4String name("npTotal");
  return name;
}
