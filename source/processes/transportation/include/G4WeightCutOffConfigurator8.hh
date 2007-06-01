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
// $Id: G4WeightCutOffConfigurator8.hh,v 1.2 2007-06-01 07:53:26 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4WeightCutOffConfigurator8
//
// Class description:
// This class builds and places the G4WeightCutOffProcess8.
// If the object is deleted the process is removed from the 
// process list.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4WeightCutOffConfigurator8_hh
#define G4WeightCutOffConfigurator8_hh G4WeightCutOffConfigurator8_hh

#include "G4Types.hh"
#include "G4VSamplerConfigurator8.hh"
#include "G4ProcessPlacer8.hh"

class G4WeightCutOffProcess8;
class G4VGCellFinder;
class G4VIStore;

class G4WeightCutOffConfigurator8 : public G4VSamplerConfigurator8
{

public:  // with description

  G4WeightCutOffConfigurator8(const G4String &particlename,
                             G4double wsurvival,
                             G4double wlimit,
                             G4double isource,
                             G4VIStore *istore,
                             const G4VGCellFinder &aGCellFinder);

  virtual ~G4WeightCutOffConfigurator8();
  virtual void Configure(G4VSamplerConfigurator8 *preConf);
  virtual const G4VTrackTerminator *GetTrackTerminator() const ;
  
private:

  G4WeightCutOffConfigurator8(const G4WeightCutOffConfigurator8&);
  G4WeightCutOffConfigurator8 &
  operator=(const G4WeightCutOffConfigurator8&);
  G4ProcessPlacer8 fPlacer;
  G4WeightCutOffProcess8 *fWeightCutOffProcess8;
  G4bool fPlaced;
};

#endif
