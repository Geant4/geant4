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
/// \file processes/phonon/include/G4PhononScattering.hh
/// \brief Definition of the G4PhononScattering class
//
// $Id: G4PhononScattering.hh 75725 2013-11-05 16:52:30Z mkelsey $
//
#ifndef G4PhononScattering_h
#define G4PhononScattering_h 1

#include "G4VPhononProcess.hh"

class G4PhononScattering : public G4VPhononProcess {
public:
  G4PhononScattering(const G4String& processName ="phononScattering" );
  virtual ~G4PhononScattering();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& );
                           
protected:
  virtual G4double GetMeanFreePath(const G4Track&, G4double,
				   G4ForceCondition* );

private:
  // hide assignment operator as private 
  G4PhononScattering(G4PhononScattering&);
  G4PhononScattering& operator=(const G4PhononScattering& right);
};

#endif










