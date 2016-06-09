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
// $Id: G4CrossSectionIonisationBorn.cc,v 1.2 2007/11/08 18:51:34 pia Exp $
// GEANT4 tag $Name: geant4-09-01 $
// 
// Contact Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// Reference: TNS Geant4-DNA paper
// Reference for implementation model: NIM. 155, pp. 145-156, 1978

// History:
// -----------
// Date         Name              Modification
// 28 Apr 2007  M.G. Pia          Created in compliance with design described in TNS paper
//
// -------------------------------------------------------------------

// Class description:
// Geant4-DNA Cross total cross section for electron elastic scattering in water
// Reference: TNS Geant4-DNA paper
// S. Chauvie et al., Geant4 physics processes for microdosimetry simulation:
// design foundation and implementation of the first set of models,
// IEEE Trans. Nucl. Sci., vol. 54, no. 6, Dec. 2007.
// Further documentation available from http://www.ge.infn.it/geant4/dna

// -------------------------------------------------------------------


#include "G4CrossSectionIonisationBorn.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4Track.hh"
#include "G4LogLogInterpolation.hh"
#include "G4SystemOfUnits.hh"


G4CrossSectionIonisationBorn::G4CrossSectionIonisationBorn()
{

  name = "IonisationBorn";
  
  // Default energy limits (defined for protection against anomalous behaviour only)
  lowEnergyLimitDefault = 25 * eV;
  highEnergyLimitDefault = 30 * keV;

  G4String fileElectron("dna/sigma_ionisation_e_born");
  G4String fileProton("dna/sigma_ionisation_p_born");

  G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();
  G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();

  G4String electron;
  G4String proton;
  
  // Factor to scale microscopic/macroscopic cross section data in water
  // ---- MGP ---- Hardcoded (taken from prototype code); to be replaced with proper calculation
  G4double scaleFactor = (1.e-22 / 3.343) * m*m;


  // Data members for electrons

  if (electronDef != 0)
    {
      electron = electronDef->GetParticleName();
      tableFile[electron] = fileElectron;

      // Energy limits
      lowEnergyLimit[electron] = 25. * eV;
      highEnergyLimit[electron] = 30. * keV;

      // Create data set with electron cross section data and load values stored in file
      G4DNACrossSectionDataSet* tableE = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
      tableE->LoadData(fileElectron);
      
      // Insert key-table pair in map
      tableData[electron] = tableE;
    }
  else
    {
      G4Exception("G4CrossSectionIonisationBorn Constructor: electron is not defined");
    }

  // Data members for protons

  if (protonDef != 0)
    {
      proton = protonDef->GetParticleName();
      tableFile[proton] = fileProton;

      // Energy limits
      lowEnergyLimit[proton] = 500. * keV;
      highEnergyLimit[proton] = 10. * MeV;

      // Create data set with proton cross section data and load values stored in file
      G4DNACrossSectionDataSet* tableP = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
      tableP->LoadData(fileProton);
      
      // Insert key-table pair in map
      tableData[proton] = tableP;
    }
  else
    {
      G4Exception("G4CrossSectionIonisationBorn Constructor: proton is not defined");
    }
}


G4CrossSectionIonisationBorn::~G4CrossSectionIonisationBorn()
{
   // Destroy the content of the map
  std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
  for (pos = tableData.begin(); pos != tableData.end(); ++pos)
    {
      G4DNACrossSectionDataSet* table = pos->second;
      delete table;
    }
}



G4double G4CrossSectionIonisationBorn::CrossSection(const G4Track& track )
{
  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double k = particle->GetKineticEnergy();
  
  // Cross section = 0 outside the energy validity limits set in the constructor
  G4double sigma = 0.;

  // ---- MGP ---- Better handling of these limits to be set in a following design iteration
  G4double lowLim = lowEnergyLimitDefault;
  G4double highLim = highEnergyLimitDefault;

  const G4String& particleName = particle->GetDefinition()->GetParticleName();

  // Retrieve energy limits for the current particle type

    std::map< G4String,G4double,std::less<G4String> >::iterator pos1;
    pos1 = lowEnergyLimit.find(particleName);

    // Lower limit
    if (pos1 != lowEnergyLimit.end())
      {
	lowLim = pos1->second;
      }

    // Upper limit
    std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
    pos2 = highEnergyLimit.find(particleName);

    if (pos2 != highEnergyLimit.end())
      {
	highLim = pos2->second;
      }

    // Verify that the current track is within the energy limits of validity of the cross section model
    if (k > lowLim && k < highLim)
      {
	std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
	pos = tableData.find(particleName);
	
	if (pos != tableData.end())
	  {
	    G4DNACrossSectionDataSet* table = pos->second;
	    if (table != 0)
	      {
		// Cross section
		sigma = table->FindValue(k);
	      }
	  }
	else
	  {
	    // The track does not corresponds to a particle pertinent to this model
	    G4Exception("G4CrossSectionIonisationBorn: attempting to calculate cross section for wrong particle");
	  }
      }

    return sigma;
}

