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

#include "globals.hh"
#include "G4XNNTotal.hh"
#include "G4XNNTotalLowE.hh"
#include "G4XPDGTotal.hh"
#include "G4VCrossSectionSource.hh"
#include "G4KineticTrack.hh"
#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"
#include "G4CrossSectionVector.hh"


G4XNNTotal::G4XNNTotal()
{ 
  components = new G4CrossSectionVector;

  G4VCrossSectionSource* xNNTotalLowE = new G4XNNTotalLowE;
  components->push_back(xNNTotalLowE);

  G4VCrossSectionSource* xNNTotalHighE = new G4XPDGTotal;
  components->push_back(xNNTotalHighE);
}


G4XNNTotal::~G4XNNTotal()
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


G4bool G4XNNTotal::operator==(const G4XNNTotal &right) const
{
  return (this == (G4XNNTotal*) &right);
}


G4bool G4XNNTotal::operator!=(const G4XNNTotal &right) const
{
  return (this != (G4XNNTotal*) &right);
}


G4String G4XNNTotal::Name() const
{
  G4String name("NNTotal");
  return name;
}
