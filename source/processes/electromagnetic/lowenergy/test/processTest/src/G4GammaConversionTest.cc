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
// $Id: G4GammaConversionTest.cc,v 1.1 2001-10-29 09:30:00 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 07 Oct 2001   MGP        Created
//
// -------------------------------------------------------------------
// Class description:
// Test DoIt method of physics processes
// Further documentation available from http://www.ge.infn.it/geant4/lowE/index.html

// -------------------------------------------------------------------

#include "globals.hh"
#include "G4GammaConversionTest.hh"
#include "G4VProcess.hh"

#include "G4LowEnergyGammaConversion.hh"
//#include "G4LowEnergyPolarizedGammaConversion.hh"

#include "G4GammaConversion.hh"

#include "G4LowEnergyBremsstrahlung.hh"
#include "G4eBremsstrahlung.hh"

#include "G4LowEnergyIonisation.hh"
#include "G4eIonisation.hh"

G4GammaConversionTest::G4GammaConversionTest(const G4String& category, G4bool isPolarised)
  :type(category), polarised(isPolarised)
{ }

G4GammaConversionTest:: ~G4GammaConversionTest()
{ }

G4VProcess* G4GammaConversionTest::createElectronIonisation()
{
  G4VProcess* process = 0;
  if (type == "lowE")
    {
      process = new G4LowEnergyIonisation;
    }
  else if (type == "standard")
    {
      process = new G4eIonisation;
    }
  return process;
}

G4VProcess* G4GammaConversionTest::createBremsstrahlung()
{
  G4VProcess* process = 0;
  if (type == "lowE")
    {
      process = new G4LowEnergyBremsstrahlung;
    }
  else if (type == "standard")
    {
      process = new G4eBremsstrahlung;
    }
  return process;  
}

G4VProcess* G4GammaConversionTest::createProcess()
{
  G4VProcess* process = 0;
  if (type == "lowE" && (!polarised))
    {
      process = new G4LowEnergyGammaConversion;
    }
  else if (type == "lowE" && polarised)
    {
      G4Exception("Polarised LowE process not available yet");
      //      process = new G4LowEnergyPolarizedGammaConversion;
    }  
  else if (type == "standard" && (!polarised))
    {
      process = new G4GammaConversion;
    }
  else if (type == "standard" && polarised)
    {
      G4Exception("Standard process not available");
    }
  return process;  
}


