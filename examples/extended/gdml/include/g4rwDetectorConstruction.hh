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
// -*- C++ -*-
// $Id: g4rwDetectorConstruction.hh,v 1.3 2006/06/29 17:21:41 gunter Exp $
#ifndef g4rwDetectorConstruction_H
#define g4rwDetectorConstruction_H 1

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

#include "Saxana/SAXProcessor.h"
#include "Saxana/ProcessingConfigurator.h"

class gogdmlDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    gogdmlDetectorConstruction();
    ~gogdmlDetectorConstruction();

  public:
    G4VPhysicalVolume* Construct();

  private:
    SAXProcessor sxp;
    ProcessingConfigurator config;
    G4VPhysicalVolume* fWorld;
};

#endif

