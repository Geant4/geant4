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
// $Id: G4ShellCorrection.cc,v 1.1 2003/08/26 18:50:45 lara Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#include "G4ShellCorrection.hh"


G4ShellCorrection* G4ShellCorrection::theInstance = 0;

G4ShellCorrection::G4ShellCorrection()
{
  theCookShellCorrections = G4CookShellCorrections::GetInstance();
  theCameronGilbertShellCorrections = G4CameronGilbertShellCorrections::GetInstance();
//  theCameronTruranHilfShellCorrections = G4CameronTruranHilfShellCorrections::GetInstance();
}


G4ShellCorrection* G4ShellCorrection::GetInstance()
{
  if (!theInstance) theInstance = new G4ShellCorrection();
  return theInstance;
}   
