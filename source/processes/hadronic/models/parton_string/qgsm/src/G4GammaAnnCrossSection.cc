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
#include "G4GammaAnnCrossSection.hh"
#include "G4Gamma.hh"

G4GammaAnnCrossSection::
G4GammaAnnCrossSection()
{
  // pho0 Nucleon
  theGammaNucXSections.push_back(new G4ASCCrossSection(113,  2212, 13.7, 35.9, 0.45, 0.079));
  // omega0 Nucleon
  theGammaNucXSections.push_back(new G4ASCCrossSection(223,  2212, 13.7, 35.9, 0.45, 0.079));
  // phi0 Nucleon
  theGammaNucXSections.push_back(new G4ASCCrossSection(333,  2212, 12.2, 26.4, 0.50, 0.079));
}

G4bool G4GammaAnnCrossSection::
InCharge(G4int aCode, G4int bCode)
{
  G4bool result = false;
  if(aCode==G4Gamma::Gamma()->GetPDGEncoding())
  {
    result=true;
  }
  else if(bCode==G4Gamma::Gamma()->GetPDGEncoding())
  {
    result = true;
  }
  return result;
}    

G4double G4GammaAnnCrossSection::
GetXsec(G4double S)
{
  G4double result = 0;
  // ratios from Phys.Lett.B40:121-126,1972; 22% assigned to higher resonances

  typedef std::vector<G4ASCCrossSection*>::iterator iter;
  iter i;
  for(i=theGammaNucXSections.begin(); i!=theGammaNucXSections.end(); i++)
  {
    result += (*i)->GetXsec(S);
  }

  // Account for higher resonances.
  result /= 0.78;

  return result;
}

