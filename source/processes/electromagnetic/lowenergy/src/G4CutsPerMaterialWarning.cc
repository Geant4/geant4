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
// $Id: G4CutsPerMaterialWarning.cc,v 1.1 2001-11-07 22:39:02 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 05 Oct 2001   MGP        Created
//
// -------------------------------------------------------------------

#include "G4CutsPerMaterialWarning.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"

void G4CutsPerMaterialWarning::PrintWarning(const G4ParticleDefinition* particle) const
{
  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  size_t nMaterials = materialTable->size();
  if (nMaterials > 1)
    {
      G4Material* material = (*materialTable)[0]; 
      G4double cut0 = particle->GetRangeThreshold(material);  
      G4double cut = cut0;
      G4bool different = false;
      size_t mat = 0;
      while ((!different) && mat < (nMaterials-1))
	{
	  mat++;
	  G4Material* material = (*materialTable)[mat]; 
	  cut = particle->GetRangeThreshold(material);
	  if (cut != cut0) different = true;
	} 
    
      
      if (different)
	{
	  G4cout << "========================== W A R N I N G ============================ " << G4endl
		 << " "                                                                      << G4endl
		 << "You are using different range thresholds for different materials"       << G4endl
		 << "This is an UNSUPPORTED feature temporarily implemented in Geant4"       << G4endl
		 << "Geant4 Low Energy Electromagnetic Physics Processes are not supported," << G4endl
		 << "if this feature is activated and you may get inconsistent results"      << G4endl
		 << "Please define the same range threshold for all materials"               << G4endl
		 << " "                                                                      << G4endl
		 << "========================== W A R N I N G ============================ " << G4endl;
	}
    }
}
