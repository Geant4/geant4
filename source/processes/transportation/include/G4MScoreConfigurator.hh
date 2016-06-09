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
// $Id: G4MScoreConfigurator.hh,v 1.5 2003/11/26 14:51:47 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
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
