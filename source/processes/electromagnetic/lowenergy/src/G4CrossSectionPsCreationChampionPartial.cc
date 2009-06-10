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
// $Id: G4CrossSectionPsCreationChampionPartial.cc,v 1.2 2009-06-10 13:32:36 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// -------------------------------------------------------------------

#include "G4CrossSectionPsCreationChampionPartial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CrossSectionPsCreationChampionPartial::G4CrossSectionPsCreationChampionPartial()
{
  numberOfPartialCrossSections=2;

  G4double scaleFactor = 1e-16*cm*cm;
  
  G4DNACrossSectionDataSet* tablePositronium1s = 
    new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV, scaleFactor);
  
  G4DNACrossSectionDataSet* tablePositronium2s = 
    new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV, scaleFactor);
  
  G4String positronium1s("positronium1s");
  G4String filePositronium1s("dna/sigma_positronium1s_creation_champion");

  G4String positronium2s("positronium2s");
  G4String filePositronium2s("dna/sigma_positronium2s_creation_champion");

  tablePositronium1s->LoadData(filePositronium1s);
  tablePositronium2s->LoadData(filePositronium2s);

  tableData[positronium1s] = tablePositronium1s;
  tableData[positronium2s] = tablePositronium2s;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CrossSectionPsCreationChampionPartial::~G4CrossSectionPsCreationChampionPartial()
{
  std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
  for (pos = tableData.begin(); pos != tableData.end(); ++pos)
  {
      G4DNACrossSectionDataSet* table = pos->second;
      delete table;
  }
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CrossSectionPsCreationChampionPartial::CrossSection
  (G4double k, G4int index, const G4ParticleDefinition* particleDefinition)
{
  if (particleDefinition != G4Positron::PositronDefinition())
    G4Exception("G4CrossSectionPsCreationChampionPartial::CrossSection: attempting to calculate cross section for wrong particle");

  G4double sigma = 0.;

  std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;

  if (index==0)
  { 
    pos = tableData.find("positronium1s");
    if (pos != tableData.end())
    {
      G4DNACrossSectionDataSet* table = pos->second;
      if (table != 0) sigma = table->FindValue(k);
    }
  }
  
  if (index==1)
  { 
    pos = tableData.find("positronium2s");
    if (pos != tableData.end())
    {
      G4DNACrossSectionDataSet* table = pos->second;
      if (table != 0) sigma = table->FindValue(k); 
    }
  }
  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4CrossSectionPsCreationChampionPartial::RandomSelectState
  (G4double k, const G4ParticleDefinition* particleDefinition)
{
  const G4int n = numberOfPartialCrossSections;
  G4double* values(new G4double[n]);
  G4double value = 0;
  G4int i = n;
  
  while (i>0)
  {
      i--;
      values[i]=CrossSection(k, i, particleDefinition);
      value+=values[i];
  }
  
  value*=G4UniformRand();
  
  i=n;
  while (i>0)
  {
      i--;
   
      if (values[i]>value)
	break;
  
      value-=values[i];
  }
  
  delete[] values;

  return i;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4CrossSectionPsCreationChampionPartial::RandomSelectShell
  (G4double k, const G4ParticleDefinition*, G4int state)
{
  G4int level = 0;

  std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;

  if (state==0) pos = tableData.find("positronium1s");
  if (state==1) pos = tableData.find("positronium2s");

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
      G4Exception("G4CrossSectionPsCreationChampionPartial::RandomSelectShell: attempting to calculate cross section for wrong particle");
  }
    
  return level;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CrossSectionPsCreationChampionPartial::Sum
  (G4double k, const G4ParticleDefinition* particleDefinition)
{
  G4double totalCrossSection = 0.;

  for (G4int i=0; i<numberOfPartialCrossSections; i++)
  {
      totalCrossSection += CrossSection(k,i,particleDefinition);
  }
  return totalCrossSection;
}
