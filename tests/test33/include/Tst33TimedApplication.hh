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
// $Id: Tst33TimedApplication.hh,v 1.5 2006-06-29 21:59:55 gunter Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class Tst33TimedApplication
//
// Class description:
//
// Provides run and event actions used for running for a certain 
// amount of time.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33TimedApplication_hh
#define Tst33TimedApplication_hh Tst33TimedApplication_hh

#include "G4Types.hh"
#include "Tst33VApplication.hh"

class Tst33TimedEventAction;

class Tst33TimedApplication : public Tst33VApplication
{
  public:
    explicit Tst33TimedApplication(G4int time);
    virtual ~Tst33TimedApplication();
  
    virtual G4UserRunAction *CreateRunAction();
    virtual Tst33VEventAction *CreateEventAction();

private:
    G4int fTime;
};

#endif
