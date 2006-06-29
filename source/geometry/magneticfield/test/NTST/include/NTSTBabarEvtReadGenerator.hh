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
//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: NTSTBabarEvtReadGenerator.hh,v 1.4 2006-06-29 18:25:21 gunter Exp $
//
// Description:
//	Class NTSTBabarEvtReadGenerator
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Bill Lockman
//
// Copyright Information:
//	Copyright (C) 2000         SCIPP, U.C. Santa Cruz
//
//------------------------------------------------------------------------

#ifndef NTSTBabarEvtReadGenerator_hh
#define NTSTBabarEvtReadGenerator_hh 1

#include <fstream>
#include "globals.hh"
#include "G4VPrimaryGenerator.hh"

class G4Event;

class NTSTBabarEvtReadGenerator:public G4VPrimaryGenerator
{
public:
  NTSTBabarEvtReadGenerator(const char* evfile);
  ~NTSTBabarEvtReadGenerator();
  
  void GeneratePrimaryVertex(G4Event* evt);
  
private:
  G4String fileName;
  std::ifstream inputFile;
};

#endif



