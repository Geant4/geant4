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
// $Id: Tst33VisApplication.hh,v 1.2 2002-10-29 16:37:10 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class Tst33VisApplication
//
// Class description:
//
// ...

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33VisApplication_hh
#define Tst33VisApplication_hh Tst33VisApplication_hh

#include "Tst33VApplication.hh"
#include "Tst33VisManager.hh"




class Tst33VisApplication : public Tst33VApplication {
public:
  Tst33VisApplication();
  virtual ~Tst33VisApplication();

  virtual G4UserRunAction *CreateRunAction();
  virtual Tst33VEventAction *CreateEventAction();
  

private:
  Tst33VisManager fVisManager;
};

#endif
