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

// Author: Ivana Hrivnacova, 04/07/2012  (ivana@ipno.in2p3.fr)

#ifndef G4HnInformation_h
#define G4HnInformation_h 1

#include "globals.hh"
#include "G4Fcn.hh" 

struct G4HnInformation
{
  G4HnInformation() 
    : fName(""), 
      fXUnitName("none"), 
      fYUnitName("none"), 
      fXFcnName("none"),
      fYFcnName("none"),
      fXUnit(1.0), 
      fYUnit(1.0), 
      fXFcn(G4FcnIdentity),
      fYFcn(G4FcnIdentity),
      fActivation(true),
      fAscii(false)
       {}

  G4HnInformation(const G4String& name, 
                  const G4String& xunitName, const G4String& yunitName,
                  const G4String& fxName, const G4String& fyName,
                  G4double xunit, G4double yunit,
                  G4Fcn fx, G4Fcn fy) 
    : fName(name),  
      fXUnitName(xunitName), 
      fYUnitName(yunitName), 
      fXFcnName(fxName),
      fYFcnName(fyName),
      fXUnit(xunit), 
      fYUnit(yunit), 
      fXFcn(fx),
      fYFcn(fy),
      fActivation(true),
      fAscii(false)
       {}

  G4String fName;
  G4String fXUnitName;
  G4String fYUnitName;
  G4String fXFcnName;
  G4String fYFcnName;
  G4double fXUnit;  
  G4double fYUnit;
  G4Fcn    fXFcn;
  G4Fcn    fYFcn;
  G4bool   fActivation;
  G4bool   fAscii;
};

#endif  
