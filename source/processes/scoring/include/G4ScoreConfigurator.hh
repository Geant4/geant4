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
// $Id: G4ScoreConfigurator.hh,v 1.1 2006/11/20 10:02:06 ahoward Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// Class G4ScoreConfigurator
//
// Class description:
// This class builds and places the G4ScoreProcess.
// If the object is deleted the process is removed from the 
// process list.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4ScoreConfigurator_hh
#define G4ScoreConfigurator_hh G4ScoreConfigurator_hh

#include "G4Types.hh"
#include "G4NewProcessPlacer.hh"

#include "G4VSamplerConfigurator.hh"

class G4VScorer;
class G4VTrackTerminator;
class G4ScoreProcess;

class G4VPhysicalVolume;

class G4ScoreConfigurator : public G4VSamplerConfigurator
{

public:  // with description

  G4ScoreConfigurator(G4VPhysicalVolume* worldvolume, 
		      const G4String &particlename,
		      G4VScorer &scorer, 
		      G4bool paraflag);

  virtual ~G4ScoreConfigurator();
  
  virtual void Configure(G4VSamplerConfigurator *preConf);
  virtual const G4VTrackTerminator *GetTrackTerminator() const;

private:

  G4ScoreConfigurator(const G4ScoreConfigurator&);
  G4ScoreConfigurator &operator=(const G4ScoreConfigurator&);

  G4VPhysicalVolume* fWorld;
  G4NewProcessPlacer fPlacer;
  G4VScorer &fScorer;
  G4ScoreProcess *fScoreProcess;


  G4bool paraflag;

};


#endif
