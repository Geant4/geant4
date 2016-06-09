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
// $Id: G4CrossSectionIonisationRudd.cc,v 1.2 2007/11/09 16:20:16 pia Exp $
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


#include "G4CrossSectionIonisationRudd.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4Track.hh"
#include "G4LogLogInterpolation.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAGenericIonsManager.hh"

G4CrossSectionIonisationRudd::G4CrossSectionIonisationRudd()
{

  name = "IonisationRudd";
  
  // Default energy limits (defined for protection against anomalous behaviour only)
  // ZERO LOW ENERGY LIMIT FOR ALLOWED KILLING
  lowEnergyLimitDefault = 0 * eV;
  highEnergyLimitDefault = 100 * MeV;

  G4String fileProton("dna/sigma_ionisation_p_rudd");
  G4String fileHydrogen("dna/sigma_ionisation_h_rudd");
  G4String fileAlphaPlusPlus("dna/sigma_ionisation_alphaplusplus_rudd");
  G4String fileAlphaPlus("dna/sigma_ionisation_alphaplus_rudd");
  G4String fileHelium("dna/sigma_ionisation_he_rudd");

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();
  G4ParticleDefinition* hydrogenDef = instance->GetIon("hydrogen");
  G4ParticleDefinition* alphaPlusPlusDef = instance->GetIon("alpha++");
  G4ParticleDefinition* alphaPlusDef = instance->GetIon("alpha+");
  G4ParticleDefinition* heliumDef = instance->GetIon("helium");

  G4String proton;
  G4String hydrogen;
  G4String alphaPlusPlus;
  G4String alphaPlus;
  G4String helium;

  // Factor to scale microscopic/macroscopic cross section data in water
  // ---- MGP ---- Hardcoded (taken from prototype code); to be replaced with proper calculation
  G4double scaleFactor = 1 * m*m;

  // Data members for protons

  if (protonDef != 0)
    {
      proton = protonDef->GetParticleName();
      tableFile[proton] = fileProton;

      // Energy limits
      lowEnergyLimit[proton] = 0. * eV;
      highEnergyLimit[proton] = 500. * keV;

      // Create data set with proton cross section data and load values stored in file
      G4DNACrossSectionDataSet* tableProton = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
									   eV,
									   scaleFactor );
      tableProton->LoadData(fileProton);
      
      // Insert key-table pair in map
      tableData[proton] = tableProton;
    }
  else
    {
      G4Exception("G4CrossSectionIonisationRudd Constructor: proton is not defined");
    }

  // Data members for hydrogen

  if (hydrogenDef != 0)
    {
      hydrogen = hydrogenDef->GetParticleName();
      tableFile[hydrogen] = fileHydrogen;

      // Energy limits
      lowEnergyLimit[hydrogen] = 0. * eV;
      highEnergyLimit[hydrogen] = 100. * MeV;

      // Create data set with hydrogen cross section data and load values stored in file
      G4DNACrossSectionDataSet* tableHydrogen = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
									     eV,
									     scaleFactor );
      tableHydrogen->LoadData(fileHydrogen);
      
      // Insert key-table pair in map
      tableData[hydrogen] = tableHydrogen;
    }
  else
    {
      G4Exception("G4CrossSectionIonisationRudd Constructor: hydrogen is not defined");
    }

  // Data members for alphaPlusPlus

  if (alphaPlusPlusDef != 0)
    {
      alphaPlusPlus = alphaPlusPlusDef->GetParticleName();
      tableFile[alphaPlusPlus] = fileAlphaPlusPlus;

      // Energy limits
      lowEnergyLimit[alphaPlusPlus] = 0. * keV;
      highEnergyLimit[alphaPlusPlus] = 10. * MeV;

      // Create data set with hydrogen cross section data and load values stored in file
      G4DNACrossSectionDataSet* tableAlphaPlusPlus = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
										  eV,
										  scaleFactor );
      tableAlphaPlusPlus->LoadData(fileAlphaPlusPlus);
      
      // Insert key-table pair in map
      tableData[alphaPlusPlus] = tableAlphaPlusPlus;
    }
  else
    {
      G4Exception("G4CrossSectionIonisationRudd Constructor: alphaPlusPlus is not defined");
    }

  // Data members for alphaPlus

  if (alphaPlusDef != 0)
    {
      alphaPlus = alphaPlusDef->GetParticleName();
      tableFile[alphaPlus] = fileAlphaPlus;

      // Energy limits
      lowEnergyLimit[alphaPlus] = 0. * keV;
      highEnergyLimit[alphaPlus] = 10. * MeV;

      // Create data set with hydrogen cross section data and load values stored in file
      G4DNACrossSectionDataSet* tableAlphaPlus = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
									      eV,
									      scaleFactor );
      tableAlphaPlus->LoadData(fileAlphaPlus);
      
      // Insert key-table pair in map
      tableData[alphaPlus] = tableAlphaPlus;
    }
  else
    {
      G4Exception("G4CrossSectionIonisationRudd Constructor: alphaPlus is not defined");
    }

  // Data members for helium

  if (heliumDef != 0)
    {
      helium = heliumDef->GetParticleName();
      tableFile[helium] = fileHelium;

      // Energy limits
      lowEnergyLimit[helium] = 0. * keV;
      highEnergyLimit[helium] = 10. * MeV;

      // Create data set with hydrogen cross section data and load values stored in file
      G4DNACrossSectionDataSet* tableHelium = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
									   eV,
									   scaleFactor );
      tableHelium->LoadData(fileHelium);
      
      // Insert key-table pair in map
      tableData[helium] = tableHelium;
    }
  else
    {
      G4Exception("G4CrossSectionIonisationRudd Constructor: helium is not defined");
    }
}


G4CrossSectionIonisationRudd::~G4CrossSectionIonisationRudd()
{
  // Destroy the content of the map
  std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
  for (pos = tableData.begin(); pos != tableData.end(); ++pos)
    {
      G4DNACrossSectionDataSet* table = pos->second;
      delete table;
    }
}



G4double G4CrossSectionIonisationRudd::CrossSection(const G4Track& track )
{
  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double k = particle->GetKineticEnergy();
  
  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  if (
      track.GetDefinition() != G4Proton::ProtonDefinition()
      &&
      track.GetDefinition() != instance->GetIon("hydrogen")
      &&
      track.GetDefinition() != instance->GetIon("alpha++")
      &&
      track.GetDefinition() != instance->GetIon("alpha+")
      &&
      track.GetDefinition() != instance->GetIon("helium")
      )
   	    
    G4Exception("G4CrossSectionIonisationRudd: attempting to calculate cross section for wrong particle");

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
  if (k >= lowLim && k <= highLim)
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


	      // BEGIN ELECTRON CORRECTION
	      // add ONE or TWO electron-water excitation for alpha+ and helium
   
	      if ( particle->GetDefinition() == instance->GetIon("alpha+") 
		   ||
		   particle->GetDefinition() == instance->GetIon("helium")
		   ) 
		{
      
		  G4DNACrossSectionDataSet* electronDataset = new G4DNACrossSectionDataSet 
		    (new G4LogLogInterpolation, eV, (1./3.343e22)*m*m);
       
		  electronDataset->LoadData("dna/sigma_ionisation_e_born");

		  G4double kElectron = k * 0.511/3728;

		  if ( particle->GetDefinition() == instance->GetIon("alpha+") ) 
		    {
		      G4double tmp1 = table->FindValue(k) + electronDataset->FindValue(kElectron);
		      delete electronDataset;
		      return tmp1;
		    }

		  if ( particle->GetDefinition() == instance->GetIon("helium") ) 
		    {
		      G4double tmp2 = table->FindValue(k) +  2. * electronDataset->FindValue(kElectron);
		      delete electronDataset;
		      return tmp2;
		    }
		}      

	      // END ELECTRON CORRECTION

	    }
	}
      else
	{
	  // The track does not corresponds to a particle pertinent to this model
	  G4Exception("G4CrossSectionIonisationRudd: attempting to calculate cross section for wrong particle");
	}
    }

  return sigma;
}

