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
//      File name:     G4CrossSectionComposite
//
//      Author:        
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4CrossSectionComposite.hh"
#include "G4VCrossSectionSource.hh"
#include "G4KineticTrack.hh"
#include "G4CrossSectionVector.hh"

G4CrossSectionComposite::G4CrossSectionComposite()
{ }


G4CrossSectionComposite::~G4CrossSectionComposite()
{ }


G4bool G4CrossSectionComposite::operator==(const G4CrossSectionComposite &right) const
{
  return (this == (G4CrossSectionComposite*) &right);
}


G4bool G4CrossSectionComposite::operator!=(const G4CrossSectionComposite &right) const
{
  return (this != (G4CrossSectionComposite*) &right);
}


G4double G4CrossSectionComposite::CrossSection(const G4KineticTrack& trk1, 
					       const G4KineticTrack& trk2) const
{
  // Cross section of composite is the sum of components cross sections

  G4double crossSection = 0.;
  G4double ecm = (trk1.Get4Momentum() + trk2.Get4Momentum()).mag();

  const G4CrossSectionVector* components = GetComponents();
  if (components != 0)
    {
      G4int nComponents = this->GetComponents()->size();
//      G4int nValid = 0;
      
      G4int i;
      for (i=0; i<nComponents; i++)
	{
	  G4CrossSectionSourcePtr componentPtr = (*components)[i];
	  G4VCrossSectionSource* component = componentPtr();
	  if (component->IsValid(ecm))
	    {
	      crossSection += component->CrossSection(trk1,trk2);
	    }
	}
    }

  return crossSection;
}


G4bool G4CrossSectionComposite::IsValid(G4double e) const
{
  // The composite is valid if any of its components are valid
  G4bool answer = false;
  const G4CrossSectionVector* components = GetComponents();
  if (components != 0)
    {
      G4int n = components->size();
      G4int i;
      for (i=0; i<n; i++)
	{
	  G4CrossSectionSourcePtr componentPtr = (*components)[i];
	  G4VCrossSectionSource* component = componentPtr();
	  if (component->IsValid(e)) 
	    {
	      answer = true;
	      break;
	    }
	}
    }
  return answer;
}

