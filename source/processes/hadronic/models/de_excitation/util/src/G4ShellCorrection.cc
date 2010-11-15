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
// $Id: G4ShellCorrection.cc,v 1.7 2010-11-15 11:47:18 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
// 

#include "G4ShellCorrection.hh"


G4ShellCorrection* G4ShellCorrection::theInstance = 0;

G4ShellCorrection::G4ShellCorrection()
{
  theCookShellCorrections = G4CookShellCorrections::GetInstance();
  theCameronGilbertShellCorrections = G4CameronGilbertShellCorrections::GetInstance();
//  theCameronTruranHilfShellCorrections = G4CameronTruranHilfShellCorrections::GetInstance();
}

G4ShellCorrection::~G4ShellCorrection()
{
}

G4ShellCorrection* G4ShellCorrection::GetInstance()
{
  if (!theInstance)  { 
    static G4ShellCorrection theCorrections;
    theInstance = &theCorrections; 
  }
  return theInstance;
}   
