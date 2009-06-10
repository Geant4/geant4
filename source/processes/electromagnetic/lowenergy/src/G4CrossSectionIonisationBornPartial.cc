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
// $Id: G4CrossSectionIonisationBornPartial.cc,v 1.6 2009-06-10 13:32:36 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4CrossSectionIonisationBornPartial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionIonisationBornPartial::G4CrossSectionIonisationBornPartial()
{
  lowEnergyLimitDefault = 12.61 * eV;
  highEnergyLimitDefault = 30 * keV;

  G4String fileElectron("dna/sigma_ionisation_e_born");
  G4String fileProton("dna/sigma_ionisation_p_born");

  G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();
  G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();

  G4String electron;
  G4String proton;
  
  G4double scaleFactor = (1.e-22 / 3.343) * m*m;

  if (electronDef != 0)
  {
    electron = electronDef->GetParticleName();
    tableFile[electron] = fileElectron;

    lowEnergyLimit[electron] = 12.61 * eV;
    highEnergyLimit[electron] = 30. * keV;

    G4DNACrossSectionDataSet* tableE = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
    tableE->LoadData(fileElectron);
      
    tableData[electron] = tableE;
  }
  else
  {
    G4Exception("G4CrossSectionIonisationBornPartial Constructor: electron is not defined");
  }

  if (protonDef != 0)
  {
    proton = protonDef->GetParticleName();
    tableFile[proton] = fileProton;

    lowEnergyLimit[proton] = 500. * keV;
    highEnergyLimit[proton] = 10. * MeV;

    G4DNACrossSectionDataSet* tableP = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
    tableP->LoadData(fileProton);
      
    tableData[proton] = tableP;
  }
  else
  {
    G4Exception("G4CrossSectionIonisationBornPartial Constructor: proton is not defined");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionIonisationBornPartial::~G4CrossSectionIonisationBornPartial()
{
  std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
  for (pos = tableData.begin(); pos != tableData.end(); ++pos)
  {
    G4DNACrossSectionDataSet* table = pos->second;
    delete table;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4CrossSectionIonisationBornPartial::RandomSelect(G4double k, const G4String& particle )
{   
  G4int level = 0;

  std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
  pos = tableData.find(particle);

  if (pos != tableData.end())
  {
    G4DNACrossSectionDataSet* table = pos->second;

    if (table != 0)
    {
      G4double* valuesBuffer = new G4double[table->NumberOfComponents()];
      const size_t n(table->NumberOfComponents());
      size_t i(n);
      G4double value = 0.;
	    
      while (i>0)
      { 
	i--;
	valuesBuffer[i] = table->GetComponent(i)->FindValue(k);
	value += valuesBuffer[i];
      }
	    
      value *= G4UniformRand();
    
      i = n;
	    
      while (i > 0)
      {
	i--;

	if (valuesBuffer[i] > value)
        {
          delete[] valuesBuffer;
          return i;
        }
	value -= valuesBuffer[i];
      }
	    
      if (valuesBuffer) delete[] valuesBuffer;
    
    }
  }
  else
  {
    G4Exception("G4CrossSectionIonisationBornPartial: attempting to calculate cross section for wrong particle");
  }
      
  return level;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CrossSectionIonisationBornPartial::CrossSection(const G4Track& track )
{
  G4double sigma = 0.;

  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double k = particle->GetKineticEnergy();
  
  G4double lowLim = lowEnergyLimitDefault;
  G4double highLim = highEnergyLimitDefault;

  const G4String& particleName = particle->GetDefinition()->GetParticleName();

  std::map< G4String,G4double,std::less<G4String> >::iterator pos1;
  pos1 = lowEnergyLimit.find(particleName);

  if (pos1 != lowEnergyLimit.end())
  {
    lowLim = pos1->second;
  }

  std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
  pos2 = highEnergyLimit.find(particleName);

  if (pos2 != highEnergyLimit.end())
  {
     highLim = pos2->second;
  }

  if (k > lowLim && k < highLim)
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
       G4Exception("G4CrossSectionIonisationBornPartial: attempting to calculate cross section for wrong particle");
     }
  }

  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CrossSectionIonisationBornPartial::Sum(G4double /* energy */, const G4String& /* particle */)
{
  return 0;
}
