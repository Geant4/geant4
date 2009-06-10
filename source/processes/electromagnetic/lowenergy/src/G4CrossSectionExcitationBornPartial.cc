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
// $Id: G4CrossSectionExcitationBornPartial.cc,v 1.4 2009-06-10 13:32:36 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4CrossSectionExcitationBornPartial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionExcitationBornPartial::G4CrossSectionExcitationBornPartial()
{
  table = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionExcitationBornPartial::~G4CrossSectionExcitationBornPartial()
{
  delete table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4CrossSectionExcitationBornPartial::RandomSelect(G4double k)
{   
  G4int level = 0;

  if (table == 0)
  {
    table = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,(1e-22/3.343)*m*m );
    table->LoadData("dna/sigma_excitation_p_born");
  }

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

  return level;
}

