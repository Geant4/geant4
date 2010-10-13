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
// $Id: G4tgbDetectorConstruction.hh,v 1.3 2010-10-13 07:56:55 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4tgbDetectorConstruction
//
// Class description:
//
// This class is the link between the user packages and G4RunManager
// It passes to the G4RunManager the top volume in the hierarchy
// constructed in G4tgbVolumeMgr.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgbDetectorConstruction_H
#define G4tgbDetectorConstruction_H 1

#include "globals.hh"

class G4VPhysicalVolume;

class G4tgbDetectorConstruction
{
 public:  // with description

   G4tgbDetectorConstruction();
   ~G4tgbDetectorConstruction();
  
   G4VPhysicalVolume* Construct();
};

#endif
