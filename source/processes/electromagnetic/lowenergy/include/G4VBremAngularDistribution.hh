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
// File name:  G4VBremAngularDistribution
//
// Author:     Pedro Rodrigues  (psilva@lip.pt)
//             Andreia Trindade (andreia@lip.pt)
//             Luis Peralta     (luis@lip.pt)
//             Design from Maria Grazia Pia (MariaGrazia.Pia@ge.infn.it) 
//            
// 
// Creation date: 21 March 2003
//
// Modifications: 
// 21 Mar 2003       A. Trindade    First implementation acording with new design
//
// Class Description: 
//
// Abstract class for Bremsstrahlung Angular Distribution Generation
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4VBremAngularDistribution_h
#define G4VBremAngularDistribution_h 1

#include "G4ios.hh"
#include "globals.hh"

class G4VBremAngularDistribution 
{

public:

  G4VBremAngularDistribution(const G4String& name);

  virtual ~G4VBremAngularDistribution();

  virtual G4double PolarAngle(const G4double initial_energy,
			      const G4double final_energy,
                              const G4int Z) = 0;

  virtual void PrintGeneratorInformation() const = 0;

protected:

private:

  // hide assignment operator 
     G4VBremAngularDistribution & operator=(const  G4VBremAngularDistribution &right);
     G4VBremAngularDistribution(const  G4VBremAngularDistribution&);

};

#endif

