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
// $Id: G4VSamplerConfigurator.hh,v 1.1 2007-06-01 08:07:00 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VNewSamplerConfigurator
//
// Class description:
//
// This is an interface for configurators setting up processes
// needed for importance sampling and scoring. 
// The Configurator may be given a pointer to another Configurator.
// If a configurator will be given a pointer to another configurator
// it may obtain a G4VTrackTerminator from the given Configurator.
// This way it is possible to delegate the killing of a track.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VNewSamplerConfigurator_hh
#define G4VNewSamplerConfigurator_hh G4VNewSamplerConfigurator_hh

#include "G4Types.hh"
#include <vector>

class G4VTrackTerminator;

class G4VNewSamplerConfigurator
{

public:  // with description

  G4VNewSamplerConfigurator();
  virtual ~G4VNewSamplerConfigurator();

  virtual void Configure(G4VNewSamplerConfigurator *preConf) = 0;
    // Do the configuration, if preConf is given a
    // G4VTrackTerminator may be obtained from it.

  virtual const G4VTrackTerminator *GetTrackTerminator() const = 0;
    // Return a G4VTrackTerminator or 0.
};

typedef std::vector<G4VNewSamplerConfigurator *> G4Configurators;

#endif
