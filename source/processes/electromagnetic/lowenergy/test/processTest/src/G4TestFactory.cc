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
// $Id: G4TestFactory.cc,v 1.4 2006-06-29 19:48:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 07 Oct 2001   MGP        Created
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4TestFactory.hh"
#include "G4ProcessTest.hh"
#include "G4ComptonTest.hh"
#include "G4GammaConversionTest.hh"
#include "G4PhotoelectricTest.hh"
#include "G4RayleighTest.hh"
#include "G4BremsstrahlungTest.hh"
#include "G4eIonisationTest.hh"

G4ProcessTest* G4TestFactory::createTestProcess(const G4String& type, 
						const G4String& category, 
						G4bool isPolarised)
{
  G4ProcessTest* test = 0;
  
  if (type == "compton")
    {
      test = new G4ComptonTest(category,isPolarised);
      G4cout << "Testing Compton scattering" << G4endl;
    }
  if (type == "conversion")
    {
      test = new G4GammaConversionTest(category,isPolarised);
      G4cout << "Testing gamma conversion" << G4endl;
    }
  if (type == "photoelectric")
    {
      test = new G4PhotoelectricTest(category);
      G4cout << "Testing photoelectric effect" << G4endl;
    }
  if (type == "rayleigh")
    {
      test = new G4RayleighTest(category,isPolarised);
      G4cout << "Testing Rayleigh scattering" << G4endl;
    }
  if (type == "bremsstrahlung")
    {
      test = new G4BremsstrahlungTest(category);
      G4cout << "Testing Bremsstrahlung" << G4endl;
    }
  if (type == "ionisation")
    {
      test = new G4eIonisationTest(category);
      G4cout << "Testing electron ionisation" << G4endl;
    }
  return test;
}



