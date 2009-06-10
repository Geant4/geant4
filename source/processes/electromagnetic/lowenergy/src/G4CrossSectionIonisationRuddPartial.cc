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
// $Id: G4CrossSectionIonisationRuddPartial.cc,v 1.5 2009-06-10 13:32:36 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4CrossSectionIonisationRuddPartial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionIonisationRuddPartial::G4CrossSectionIonisationRuddPartial()
{
  lowEnergyLimitDefault = 100 * eV;
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

  G4double scaleFactor = 1 * m*m;

  if (protonDef != 0)
  {
    proton = protonDef->GetParticleName();
    tableFile[proton] = fileProton;

    lowEnergyLimit[proton] = 100. * eV;
    highEnergyLimit[proton] = 500. * keV;

    G4DNACrossSectionDataSet* tableProton = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
    tableProton->LoadData(fileProton);
      
    tableData[proton] = tableProton;
  }
  else
  {
    G4Exception("G4CrossSectionIonisationRudd Constructor: proton is not defined");
  }

  if (hydrogenDef != 0)
  {
    hydrogen = hydrogenDef->GetParticleName();
    tableFile[hydrogen] = fileHydrogen;

    lowEnergyLimit[hydrogen] = 100. * eV;
    highEnergyLimit[hydrogen] = 100. * MeV;

    G4DNACrossSectionDataSet* tableHydrogen = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
    tableHydrogen->LoadData(fileHydrogen);
      
    tableData[hydrogen] = tableHydrogen;
  }
  else
  {
    G4Exception("G4CrossSectionIonisationRudd Constructor: hydrogen is not defined");
  }

  if (alphaPlusPlusDef != 0)
  {
    alphaPlusPlus = alphaPlusPlusDef->GetParticleName();
    tableFile[alphaPlusPlus] = fileAlphaPlusPlus;

    lowEnergyLimit[alphaPlusPlus] = 1. * keV;
    highEnergyLimit[alphaPlusPlus] = 10. * MeV;

    G4DNACrossSectionDataSet* tableAlphaPlusPlus = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
    tableAlphaPlusPlus->LoadData(fileAlphaPlusPlus);
      
    tableData[alphaPlusPlus] = tableAlphaPlusPlus;
  }
  else
  {
    G4Exception("G4CrossSectionIonisationRudd Constructor: alphaPlusPlus is not defined");
  }

  if (alphaPlusDef != 0)
  {
    alphaPlus = alphaPlusDef->GetParticleName();
    tableFile[alphaPlus] = fileAlphaPlus;

    lowEnergyLimit[alphaPlus] = 1. * keV;
    highEnergyLimit[alphaPlus] = 10. * MeV;

    G4DNACrossSectionDataSet* tableAlphaPlus = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
    tableAlphaPlus->LoadData(fileAlphaPlus);

    tableData[alphaPlus] = tableAlphaPlus;
  }
  else
  {
    G4Exception("G4CrossSectionIonisationRudd Constructor: alphaPlus is not defined");
  }

  if (heliumDef != 0)
  {
    helium = heliumDef->GetParticleName();
    tableFile[helium] = fileHelium;

    lowEnergyLimit[helium] = 1. * keV;
    highEnergyLimit[helium] = 10. * MeV;

    G4DNACrossSectionDataSet* tableHelium = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
    tableHelium->LoadData(fileHelium);
      
    tableData[helium] = tableHelium;
  }
  else
  {
    G4Exception("G4CrossSectionIonisationRudd Constructor: helium is not defined");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionIonisationRuddPartial::~G4CrossSectionIonisationRuddPartial()
{
  std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
  for (pos = tableData.begin(); pos != tableData.end(); ++pos)
  {
    G4DNACrossSectionDataSet* table = pos->second;
    delete table;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4CrossSectionIonisationRuddPartial::RandomSelect(G4double k, const G4String& particle )
{   
  
  // BEGIN PART 1/2 OF ELECTRON CORRECTION

  // add ONE or TWO electron-water excitation for alpha+ and helium
   
  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  G4double kElectron(0);
  G4double electronComponent(0);
  G4DNACrossSectionDataSet * electronDataset = new G4DNACrossSectionDataSet (new G4LogLogInterpolation, eV, (1./3.343e22)*m*m);
 
  if ( particle == instance->GetIon("alpha+")->GetParticleName()
       ||
       particle == instance->GetIon("helium")->GetParticleName()
       ) 
  {     
      electronDataset->LoadData("dna/sigma_ionisation_e_born");

      kElectron = k * 0.511/3728;
       
      electronComponent = electronDataset->FindValue(kElectron);
       
  }      
  
  delete electronDataset;
  
  // END PART 1/2 OF ELECTRON CORRECTION
  
  G4int level = 0;

  // Retrieve data table corresponding to the current particle type  

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

	      // BEGIN PART 2/2 OF ELECTRON CORRECTION

	      if (particle == instance->GetIon("alpha+")->GetParticleName()) 
		{valuesBuffer[i]=table->GetComponent(i)->FindValue(k) + electronComponent; }
     
	      if (particle == instance->GetIon("helium")->GetParticleName()) 
		{valuesBuffer[i]=table->GetComponent(i)->FindValue(k) + 2*electronComponent; }
      
	      // BEGIN PART 2/2 OF ELECTRON CORRECTION

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
    G4Exception("G4CrossSectionIonisationRuddPartial: attempting to calculate cross section for wrong particle");
  }
      
  return level;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CrossSectionIonisationRuddPartial::CrossSection(const G4Track& track )
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

  if (k >= lowLim && k <= highLim)
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
	  G4Exception("G4CrossSectionIonisationRuddPartial: attempting to calculate cross section for wrong particle");
      }
  }

  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CrossSectionIonisationRuddPartial::Sum(G4double /* energy */, const G4String& /* particle */)
{
  return 0;
}
