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
// $Id: Tst33TimedApplication.cc,v 1.2 2002-10-29 16:37:10 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33TimedApplication.cc
//
// ----------------------------------------------------------------------

#include "Tst33TimedApplication.hh"
#include "G4RunManager.hh"
#include "Tst33TimedEventAction.hh"



Tst33TimedApplication::Tst33TimedApplication(G4int time)
  :
  fTime(time)
{}

Tst33TimedApplication::~Tst33TimedApplication()
{}

Tst33VEventAction *Tst33TimedApplication::CreateEventAction() {
  return new Tst33TimedEventAction(fTime);
}



G4UserRunAction *Tst33TimedApplication::CreateRunAction(){
  return 0;
}

