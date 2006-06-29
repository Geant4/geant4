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
// $Id: Tst33VisApplication.cc,v 1.3 2006-06-29 22:01:38 gunter Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33VisApplication.cc
//
// ----------------------------------------------------------------------

#include "Tst33VisApplication.hh"
#include "G4VisManager.hh"
#include "Tst33VisManager.hh"
#include "Tst33VisEventAction.hh"
#include "Tst33VisRunAction.hh"

#include "G4RunManager.hh"

#include "G4UImanager.hh"


Tst33VisApplication::  Tst33VisApplication() 
{
  fVisManager.Initialize();
}

Tst33VisApplication::~Tst33VisApplication(){
}

Tst33VEventAction *Tst33VisApplication::CreateEventAction() {
  return new Tst33VisEventAction;
}



G4UserRunAction *Tst33VisApplication::CreateRunAction(){
  return new Tst33VisRunAction;
}


