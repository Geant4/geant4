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
// $Id: Tst33VisApplication.cc,v 1.2 2002-10-29 16:37:10 dressel Exp $
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


