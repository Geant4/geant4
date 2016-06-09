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
// $Id: G4PScoreConfigurator.hh,v 1.6 2006/06/29 21:10:02 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// ----------------------------------------------------------------------
// Class G4PScoreConfigurator
//
// Class description:
// This class builds and places the G4PScoreProcess.
// If the object is deleted the process is removed from the 
// process list.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4PScoreConfigurator_hh
#define G4PScoreConfigurator_hh G4PScoreConfigurator_hh

#include "G4Types.hh"
#include "G4VSamplerConfigurator.hh"
#include "G4ProcessPlacer.hh"

class G4VParallelStepper;
class G4VScorer;
class G4VTrackTerminator;
class G4PScoreProcess;

class G4PScoreConfigurator : public G4VSamplerConfigurator
{

public:  // with description

  G4PScoreConfigurator(const G4String &particlename,
                       G4VParallelStepper &pstepper,
                       G4VScorer &scorer);
  virtual ~G4PScoreConfigurator();

  virtual void Configure(G4VSamplerConfigurator *preConf);
  virtual const G4VTrackTerminator *GetTrackTerminator() const;  
  
private:

  G4PScoreConfigurator(const G4PScoreConfigurator &);
  G4PScoreConfigurator &operator=(const G4PScoreConfigurator &);
  G4ProcessPlacer fPlacer;
  G4VParallelStepper &fPStepper;
  G4VScorer &fScorer;
  G4PScoreProcess *fPScoreProcess;
};


#endif
