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
      std::size_t nComponents = GetComponents()->size();
      for (std::size_t i=0; i<nComponents; ++i)
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
      std::size_t n = components->size();
      for (std::size_t i=0; i<n; ++i)
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

