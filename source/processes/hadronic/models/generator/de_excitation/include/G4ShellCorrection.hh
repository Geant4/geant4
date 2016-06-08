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
// $Id: G4ShellCorrection.hh,v 1.2 2002/12/12 19:17:11 gunter Exp $
// GEANT4 tag $Name: geant4-05-00 $
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
  ~G4ShellCorrection() {};

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
