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
//      File name:     G4XnpElastic
//
//      Author:        
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4XnpElastic.hh"
#include "G4XnpElasticLowE.hh"
#include "G4XPDGElastic.hh"
#include "G4VCrossSectionSource.hh"
#include "G4KineticTrack.hh"
#include "G4ParticleDefinition.hh"
#include "G4CrossSectionVector.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

G4XnpElastic::G4XnpElastic()
{ 
  components = new G4CrossSectionVector;

  G4VCrossSectionSource* xnpElasticLowE = new G4XnpElasticLowE;
  components->push_back(xnpElasticLowE);

//  G4ParticleDefinition* proton = G4Proton::ProtonDefinition();
//  G4ParticleDefinition* neutron = G4Neutron::NeutronDefinition();
  G4VCrossSectionSource* xnpElasticHighE = new G4XPDGElastic;
  components->push_back(xnpElasticHighE);
}


G4XnpElastic::~G4XnpElastic()
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


G4bool G4XnpElastic::operator==(const G4XnpElastic &right) const
{
  return (this == (G4XnpElastic*) &right);
}


G4bool G4XnpElastic::operator!=(const G4XnpElastic &right) const
{
  return (this != (G4XnpElastic*) &right);
}


G4String G4XnpElastic::Name() const
{
  G4String name("npElastic");
  return name;
}

