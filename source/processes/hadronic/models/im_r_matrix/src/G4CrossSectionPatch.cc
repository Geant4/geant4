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
//      File name:     G4CrossSectionPatch
//
//      Author:        
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4CrossSectionPatch.hh"
#include "G4VCrossSectionSource.hh"
#include "G4KineticTrack.hh"
#include "G4CrossSectionVector.hh"

G4CrossSectionPatch::G4CrossSectionPatch()
{ }


G4CrossSectionPatch::~G4CrossSectionPatch()
{ }


G4bool G4CrossSectionPatch::operator==(const G4CrossSectionPatch &right) const
{
  return (this == (G4CrossSectionPatch*) &right);
}


G4bool G4CrossSectionPatch::operator!=(const G4CrossSectionPatch &right) const
{
  return (this != (G4CrossSectionPatch*) &right);
}


G4double G4CrossSectionPatch::CrossSection(const G4KineticTrack& trk1, 
					   const G4KineticTrack& trk2) const
{
  // The cross section is provided by one of the components, according to their energy
  // validity range

  G4double crossSection = 0.;
  G4double ecm = (trk1.Get4Momentum() + trk2.Get4Momentum()).mag();

  const G4CrossSectionVector* components = GetComponents();
  if (components != 0)
    {
      std::size_t nComponents = this->GetComponents()->size();
      
      for (std::size_t i=0; i<nComponents; ++i)
	{
          G4CrossSectionSourcePtr componentPtr = (*components)[i];
	  G4VCrossSectionSource* component = componentPtr();
	  if (component->IsValid(ecm))
	    {
	      crossSection = component->CrossSection(trk1,trk2);
	    }
	  else if (i < (nComponents - 1) )
	    {
	      G4CrossSectionSourcePtr nextPtr = (*components)[i+1];
	      G4VCrossSectionSource* next = nextPtr();
	      if (ecm > component->HighLimit() && ecm < next->LowLimit())
		{
		  // Merge cross-sections in transition region between two validity ranges
		  crossSection = Transition(trk1,trk2,component,next);
		}
	    }
	}
    }

  return crossSection;
}


G4bool G4CrossSectionPatch::IsValid(G4double e) const
{
  // The Patch is valid if any of its components are valid
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


G4double G4CrossSectionPatch::Transition(const G4KineticTrack& trk1, const G4KineticTrack& trk2,
					      const G4VCrossSectionSource* comp1, 
					      const G4VCrossSectionSource* comp2) const
{
  //Merge two cross sections in the transition region between their validity ranges

  G4double crossSection = 0.;

  G4double ecm = (trk1.Get4Momentum() + trk2.Get4Momentum()).mag();
  G4double sigma1 = comp1->CrossSection(trk1,trk2);
  G4double sigma2 = comp2->CrossSection(trk1,trk2);
  G4double denom = comp2->LowLimit() - comp1->HighLimit();
  G4double diff = ecm - comp1->HighLimit();
  if (denom > 0. && diff > 0.)
    {
      G4double ratio = diff / denom;
      crossSection = (1.- ratio) * sigma1 + ratio * sigma2;
    }

  return crossSection;
}


G4double G4CrossSectionPatch::Transition(G4double ecm, 
					 G4double sigma1, G4double sigma2, 
					 G4double e1, G4double e2) const
{
  //Merge two cross sections in the transition region between their validity ranges

  G4double crossSection = 0.;
  
  G4double denom = e2 - e1;
  G4double diff = ecm - e1;
  if (denom > 0. && diff > 0.)
    {
      G4double ratio = diff / denom;
      crossSection = (1.- ratio) * sigma1 + ratio * sigma2;
    }

  return crossSection;
}

