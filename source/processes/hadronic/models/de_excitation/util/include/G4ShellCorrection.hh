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
// $Id$
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4ShellCorrection_h
#define G4ShellCorrection_h 1

#include "globals.hh"
#include "G4CookShellCorrections.hh"
#include "G4CameronGilbertShellCorrections.hh"
//#include "G4CameronTruranHilfShellCorrections.hh"


class G4ShellCorrection
{
private:
  
  // Dummy constructor
  G4ShellCorrection();
  
  static G4ShellCorrection* theInstance;
  
public:
  static G4ShellCorrection* GetInstance();
  ~G4ShellCorrection();

  G4double GetShellCorrection(const G4int A, const G4int Z) const 
  {
    G4double SCorrection = 0.0;
    if (theCookShellCorrections->IsInTableThisN(A-Z) || 
	theCookShellCorrections->IsInTableThisZ(Z)) 
      SCorrection = theCookShellCorrections->GetShellCorrection(A,Z);
    else if (theCameronGilbertShellCorrections->IsInTableThisN(A-Z) || 
	     theCameronGilbertShellCorrections->IsInTableThisZ(Z))
      SCorrection = theCameronGilbertShellCorrections->GetShellCorrection(A,Z);
    
    return SCorrection;
  }

private:

  G4CookShellCorrections* theCookShellCorrections;
  G4CameronGilbertShellCorrections* theCameronGilbertShellCorrections;
//  G4CameronTruranHilfShellCorrections* theCameronTruranHilfShellCorrections;


};
#endif
