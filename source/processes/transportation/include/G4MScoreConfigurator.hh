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
// $Id: G4MScoreConfigurator.hh,v 1.6 2006/06/29 21:09:46 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// Class G4MScoreConfigurator
//
// Class description:
// This class builds and places the G4MScoreProcess.
// If the object is deleted the process is removed from the 
// process list.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4MScoreConfigurator_hh
#define G4MScoreConfigurator_hh G4MScoreConfigurator_hh

#include "G4Types.hh"
#include "G4ProcessPlacer.hh"

#include "G4VSamplerConfigurator.hh"

class G4VScorer;
class G4VTrackTerminator;
class G4MScoreProcess;

class G4MScoreConfigurator : public G4VSamplerConfigurator
{

public:  // with description

  G4MScoreConfigurator(const G4String &particlename,
                       G4VScorer &scorer);
  virtual ~G4MScoreConfigurator();
  
  virtual void Configure(G4VSamplerConfigurator *preConf);
  virtual const G4VTrackTerminator *GetTrackTerminator() const;

private:

  G4MScoreConfigurator(const G4MScoreConfigurator&);
  G4MScoreConfigurator &operator=(const G4MScoreConfigurator&);

  G4ProcessPlacer fPlacer;
  G4VScorer &fScorer;
  G4MScoreProcess *fMScoreProcess;
};


#endif
