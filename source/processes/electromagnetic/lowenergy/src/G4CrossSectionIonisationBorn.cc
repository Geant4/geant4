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
// $Id: G4CrossSectionIonisationBorn.cc,v 1.7 2009-06-11 15:47:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4CrossSectionIonisationBorn.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionIonisationBorn::G4CrossSectionIonisationBorn()
{
  lowEnergyLimitDefault = 12.61 * eV; // SI: i/o 25 eV
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

    lowEnergyLimit[electron] = 12.61 * eV; // SI: i/o 25 eV
    highEnergyLimit[electron] = 30. * keV;

    G4DNACrossSectionDataSet* tableE = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
    tableE->LoadData(fileElectron);
      
    tableData[electron] = tableE;
  }
  else
  {
    G4Exception("G4CrossSectionIonisationBorn Constructor: electron is not defined");
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
    G4Exception("G4CrossSectionIonisationBorn Constructor: proton is not defined");
  }

   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class G4CrossSectionIonisationBorn is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionIonisationBorn::~G4CrossSectionIonisationBorn()
{
  std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
  for (pos = tableData.begin(); pos != tableData.end(); ++pos)
  {
    G4DNACrossSectionDataSet* table = pos->second;
    delete table;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CrossSectionIonisationBorn::CrossSection(const G4Track& track )
{
  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double k = particle->GetKineticEnergy();
  
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
      G4Exception("G4CrossSectionIonisationBorn: attempting to calculate cross section for wrong particle");
    }
  }

  return sigma;
}

