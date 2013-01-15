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
// $Id$
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
// 

#include "G4ShellCorrection.hh"


G4ThreadLocal G4ShellCorrection* G4ShellCorrection::theInstance = 0;

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
    static G4ThreadLocal G4ShellCorrection *theCorrections_G4MT_TLS_ = 0 ; if (!theCorrections_G4MT_TLS_) theCorrections_G4MT_TLS_ = new  G4ShellCorrection  ;  G4ShellCorrection &theCorrections = *theCorrections_G4MT_TLS_;
    theInstance = &theCorrections; 
  }
  return theInstance;
}   
