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
// $Id: G4CrossSectionIonisationRudd.cc,v 1.7 2009-06-11 15:47:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4CrossSectionIonisationRudd.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionIonisationRudd::G4CrossSectionIonisationRudd()
{
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

  G4double scaleFactor = 1 * m*m;

  if (protonDef != 0)
  {
    proton = protonDef->GetParticleName();
    tableFile[proton] = fileProton;

    lowEnergyLimit[proton] = 0. * eV;
    highEnergyLimit[proton] = 500. * keV;

    G4DNACrossSectionDataSet* tableProton = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
									   eV,
									   scaleFactor );
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

    lowEnergyLimit[hydrogen] = 0. * eV;
    highEnergyLimit[hydrogen] = 100. * MeV;

    G4DNACrossSectionDataSet* tableHydrogen = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
									     eV,
									     scaleFactor );
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

    lowEnergyLimit[alphaPlusPlus] = 0. * eV;
    highEnergyLimit[alphaPlusPlus] = 10. * MeV;

    G4DNACrossSectionDataSet* tableAlphaPlusPlus = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
										  eV,
										  scaleFactor );
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

    lowEnergyLimit[alphaPlus] = 0. * eV;
    highEnergyLimit[alphaPlus] = 10. * MeV;

    G4DNACrossSectionDataSet* tableAlphaPlus = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
									      eV,
									      scaleFactor );
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

    lowEnergyLimit[helium] = 0. * eV;
    highEnergyLimit[helium] = 10. * MeV;

    G4DNACrossSectionDataSet* tableHelium = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
									   eV,
									   scaleFactor );
    tableHelium->LoadData(fileHelium);
    tableData[helium] = tableHelium;
    }
  else
  {
    G4Exception("G4CrossSectionIonisationRudd Constructor: helium is not defined");
  }

   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class G4CrossSectionIonisationRudd is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionIonisationRudd::~G4CrossSectionIonisationRudd()
{
  std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
  for (pos = tableData.begin(); pos != tableData.end(); ++pos)
  {
    G4DNACrossSectionDataSet* table = pos->second;
    delete table;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

  G4double sigma = 0.;

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
	G4Exception("G4CrossSectionIonisationRudd: attempting to calculate cross section for wrong particle");
    }
  }

  return sigma;
}

