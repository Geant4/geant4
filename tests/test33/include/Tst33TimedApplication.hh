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
// $Id: Tst33TimedApplication.hh,v 1.3 2002-11-20 13:09:16 dressel Exp $
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

#include "Tst33VApplication.hh"
#include "globals.hh"

class Tst33TimedEventAction;


class Tst33TimedApplication : public Tst33VApplication {
public:
  explicit Tst33TimedApplication(G4int time);
  virtual ~Tst33TimedApplication();
  
  virtual G4UserRunAction *CreateRunAction();
  virtual Tst33VEventAction *CreateEventAction();
private:
  G4int fTime;
};

#endif
