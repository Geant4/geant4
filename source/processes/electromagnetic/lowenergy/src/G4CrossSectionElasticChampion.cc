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
// $Id: G4CrossSectionElasticChampion.cc,v 1.7 2009-06-11 15:47:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// -------------------------------------------------------------------

#include "G4CrossSectionElasticChampion.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CrossSectionElasticChampion::G4CrossSectionElasticChampion()
{
  lowEnergyLimit = 0.*eV; 
  highEnergyLimit = 10. * keV;

  G4double scaleFactor = 1e-16*cm*cm;

  G4String fileElectron("dna/sigma_elastic_e_champion");

  G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();
  G4String electron;
  
  if (electronDef != 0)
  {
    electron = electronDef->GetParticleName();
    tableFile[electron] = fileElectron;

    G4DNACrossSectionDataSet* tableE = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
    tableE->LoadData(fileElectron);
    tableData[electron] = tableE;
  }
  else
  {
    G4Exception("G4CrossSectionElasticChampion constructor: electron is not defined");
  }

   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class G4CrossSectionElasticChampion is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CrossSectionElasticChampion::~G4CrossSectionElasticChampion()
{
  std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
  for (pos = tableData.begin(); pos != tableData.end(); ++pos)
  {
    G4DNACrossSectionDataSet* table = pos->second;
    delete table;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CrossSectionElasticChampion::CrossSection(const G4Track& track )
{
  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double k = particle->GetKineticEnergy();
  
  G4double sigma = 0.;

  const G4String& particleName = particle->GetDefinition()->GetParticleName();
 
  if (k >= lowEnergyLimit && k < highEnergyLimit)
  {
	std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
	pos = tableData.find(particleName);
	
	if (pos != tableData.end())
	{
	  G4DNACrossSectionDataSet* table = pos->second;
	  if (table != 0)
	  {
	    sigma = table->FindValue(k);
	  }
	}
	else
	{
	    G4Exception("G4CrossSectionElasticChampion::CrossSection: attempting to calculate cross section for wrong particle");
	}
  }
  return sigma;
}
