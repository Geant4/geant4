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
// $Id: G4Vee2hadrons.hh 97391 2016-06-02 10:08:45Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4Vee2hadrons
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 12.08.2004
//
// Modifications: 14 July 2014 N. Chikuma revised interfaces  
//

//
// Class Description: base class to compute partial cross sections
//                    of e+e- annihilation into hadrons and 
//                    sample of final state in the centre mass frame

// -------------------------------------------------------------------
//
#ifndef G4Vee2hadrons_h
#define G4Vee2hadrons_h 1

#include <vector>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4eeCrossSections.hh"
#include "G4PhysicsLinearVector.hh"

class G4DynamicParticle;
class G4PhysicsVector;

class G4Vee2hadrons 
{

public:

  explicit G4Vee2hadrons(G4eeCrossSections* cr,
		G4double vlowEnergy,
		G4double vhighEnergy,
		G4double vdelta) : cross(cr)
  {
	lowEnergy  = vlowEnergy;
	highEnergy = vhighEnergy;
	delta      = vdelta;
  };

  virtual ~G4Vee2hadrons() {};

  virtual G4double PeakEnergy() const = 0;

  virtual G4double ComputeCrossSection(G4double) const = 0;

  G4PhysicsVector* PhysicsVector() const
  {
    G4int nbins = std::max(3, G4int((highEnergy - lowEnergy)/delta) );
    G4PhysicsVector* pp = new G4PhysicsLinearVector(lowEnergy,highEnergy,nbins);
    return pp;
  };

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 G4double, const G4ThreeVector&) = 0;

  G4double LowEnergy() const {return lowEnergy;};

  G4double HighEnergy() const {return highEnergy;};
  
private:

  // hide assignment operator
  G4Vee2hadrons & operator=(const  G4Vee2hadrons &right) = delete;
  G4Vee2hadrons(const  G4Vee2hadrons&) = delete;

  // parameters of the table
  G4double lowEnergy;
  G4double highEnergy;
  G4double delta;       

protected:

   G4eeCrossSections* cross;  // class to compute cross section

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
