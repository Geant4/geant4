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
// $Id: Tst33VApplication.hh,v 1.4 2006-06-29 22:00:01 gunter Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class Tst33VApplication
//
// Class description:
//
// Base class for an application.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33VApplication_hh
#define Tst33VApplication_hh Tst33VApplication_hh


class G4VCellScorer;
class G4UserRunAction;
class Tst33VEventAction;

class Tst33VApplication {
public:
  Tst33VApplication();
  virtual ~Tst33VApplication();

  virtual G4UserRunAction *CreateRunAction() = 0;
  virtual Tst33VEventAction *CreateEventAction() = 0;
};

#endif
