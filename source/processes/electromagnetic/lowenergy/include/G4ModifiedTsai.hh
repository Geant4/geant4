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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:  G4ModifiedTsai
//
// Author:     Andreia Trindade (andreia@lip.pt)
//             Pedro Rodrigues  (psilva@lip.pt)
//             Luis Peralta     (luis@lip.pt)
// 
// Creation date: 21 March 2003
//
// Modifications: 
// 21 Mar 2003       A. Trindade    First implementation acording with new design
// 24 Mar 2003                      & Fix in Tsai generator in order to prevent theta generation above pi
//               
//
// Class Description: 
//
// Concrete class for Bremsstrahlung Angular Distribution Generation - Tsai Model
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4ModifiedTsai_h
#define G4ModifiedTsai_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VBremAngularDistribution.hh"

class G4ModifiedTsai : public G4VBremAngularDistribution
{

public:

  G4ModifiedTsai(const G4String& name);

  ~G4ModifiedTsai();

  G4double PolarAngle(const G4double initial_energy,
			      const G4double initial_momentum,
			      const G4double final_energy,
			      const G4double final_momentum,
			      const G4int Z) const;

  void PrintGeneratorInformation() const;

protected:

private:

  // hide assignment operator 
     G4ModifiedTsai & operator=(const  G4ModifiedTsai &right);
     G4ModifiedTsai(const  G4ModifiedTsai&);

};

#endif

