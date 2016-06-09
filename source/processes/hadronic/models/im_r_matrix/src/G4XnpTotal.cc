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
  if (components != 0) 
    {
      G4int nComponents = this->GetComponents()->size();
      G4int i;
      for (i=0; i<nComponents; i++)
	{
          G4CrossSectionSourcePtr componentPtr = (*components)[i];
	  G4VCrossSectionSource* component = componentPtr();
	  delete component;
	  component = 0;
	  componentPtr = 0;
	}
    }
  delete components;
  components = 0;
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
