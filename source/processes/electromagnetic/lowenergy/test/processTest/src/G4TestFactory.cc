//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4TestFactory.cc,v 1.1 2001-10-29 09:30:01 pia Exp $
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

const G4ProcessTest* G4TestFactory::createTestProcess(const G4String& type, 
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



