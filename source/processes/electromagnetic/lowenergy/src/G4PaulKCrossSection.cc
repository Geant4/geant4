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
// History:
// -----------
//  21 Apr 2008   H. Abdelohauwed - 1st implementation
//  29 Apr 2009   ALF  Major Design Revision
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics, Cross section, p ionisation, K shell
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>
//#include "G4CompositeEMDataSet.hh"
//#include "G4ShellEMDataSet.hh"
#include "G4EMDataSet.hh"
//#include "G4VEMDataSet.hh"
//#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "G4PaulKCrossSection.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"


G4PaulKCrossSection::G4PaulKCrossSection()
{ 

  
  interpolation = new G4LogLogInterpolation();

  /*
    G4String path = getenv("G4LEDATA");
 
    if (!path)
    G4Exception("G4paulKCrossSection::G4paulKCrossSection: G4LEDATA environment variable not set");
    G4cout << path + "/kcsPaul/kcs-" << G4endl;
  */


    for (G4int i=4; i<93; i++) {
      protonDataSetMap[i] = new G4EMDataSet(i,interpolation);
      protonDataSetMap[i]->LoadData("pixe/kpcsPaul/kcs-");
    }
    for (G4int i=6; i<93; i++) {
      alphaDataSetMap[i] = new G4EMDataSet(i,interpolation);
      alphaDataSetMap[i]->LoadData("pixe/kacsPaul/kacs-");
    }




}

G4PaulKCrossSection::~G4PaulKCrossSection()
{ 

  protonDataSetMap.clear();
  alphaDataSetMap.clear();

}

G4double G4PaulKCrossSection::CalculateKCrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)
{
  
  G4Proton* aProtone = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();
  
  G4double sigma = 0;

  if (massIncident == aProtone->GetPDGMass() )
    {
      
      sigma = protonDataSetMap[zTarget]->FindValue(energyIncident/MeV); 
      
    }
  else
    {
      if (massIncident == aAlpha->GetPDGMass())
	{
	  
          sigma = alphaDataSetMap[zTarget]->FindValue(energyIncident/MeV); 
	  
	}
      else
	{ 
	  G4cout << "we can treat only Proton or Alpha incident particles " << G4endl;
	  sigma = 0.;
	}
    }
  
  
  // sigma is in internal units (mm^2)
  return sigma;
}












