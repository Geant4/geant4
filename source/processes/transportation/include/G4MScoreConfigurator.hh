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
// $Id: G4MScoreConfigurator.hh,v 1.2 2002-10-16 16:26:58 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4MScoreConfigurator
//
// Class description:
// This class builds and places  the G4MScoreProcess.
// If the object is deleted the process is removed from the 
// process list.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4MScoreConfigurator_hh
#define G4MScoreConfigurator_hh G4MScoreConfigurator_hh

#include "globals.hh"
#include "G4MScoreProcess.hh"
#include "G4ProcessPlacer.hh"

#include "G4VSamplerConfigurator.hh"

class G4VPScorer;
class G4VTrackTerminator;

class G4MScoreConfigurator : public G4VSamplerConfigurator{
public:
  G4MScoreConfigurator(const G4String &particlename,
		       G4VPScorer &scorer);
  virtual ~G4MScoreConfigurator();
  
  virtual void Configure(G4VSamplerConfigurator *preConf);
  virtual const G4VTrackTerminator *GetTrackTerminator() const;

  
private:
  G4MScoreConfigurator(const G4MScoreConfigurator&);
  const G4MScoreConfigurator &operator=(const G4MScoreConfigurator&);

  G4ProcessPlacer fPlacer;
  G4VPScorer &fScorer;
  G4MScoreProcess *fMScoreProcess;
};


#endif
