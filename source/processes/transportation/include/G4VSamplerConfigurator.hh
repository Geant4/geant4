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
// ********************************************************************
//
//
// $Id: G4VSamplerConfigurator.hh,v 1.3 2002-12-12 19:18:40 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VSamplerConfigurator
//
// Class description:
//
// This is an interface for configurators seeting up processes
// needed for importance sampling and scoring. 
// The Configurator may be given a pointer to another Configurator.
// If a configurator will be given a pointer to another configurator
// it may obtain a G4VTrackTerminator from the given Configurator.
// This way it is possible to delegate the killing of a track.
// 
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4VSamplerConfigurator_hh
#define G4VSamplerConfigurator_hh G4VSamplerConfigurator_hh

#include "globals.hh"
#include "g4std/vector"
class G4VTrackTerminator;


class G4VSamplerConfigurator{
public:
  G4VSamplerConfigurator();
  virtual ~G4VSamplerConfigurator();
  virtual void Configure(G4VSamplerConfigurator *preConf) = 0;
    // do the configuration, if preConf is given a
    //  G4VTrackTerminator may be obtained from it
  virtual const G4VTrackTerminator *GetTrackTerminator() const = 0;
    // return a G4VTrackTerminator or 0
};

typedef G4std::vector<G4VSamplerConfigurator *> G4Configurators;


#endif
