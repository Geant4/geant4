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
// File name:  G4Generator2BS
//
// Author:     Andreia Trindade (andreia@lip.pt)
//             Pedro Rodrigues  (psilva@lip.pt)
//             Luis Peralta     (luis@lip.pt)
// 
// Creation date: 2 June 2003
//
// Modifications: 
// 02 Jun 2003                           First implementation acording with new design
//               
//
// Class Description: 
//
// Concrete class for Bremsstrahlung Angular Distribution Generation - 2BS Distribution
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4Generator2BS_h
#define G4Generator2BS_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VBremAngularDistribution.hh"

class G4Generator2BS : public G4VBremAngularDistribution
{

public:

  G4Generator2BS(const G4String& name);

  ~G4Generator2BS();

  G4double PolarAngle(const G4double initial_energy,
		      const G4double final_energy,
		      const G4int Z);

  void PrintGeneratorInformation() const;

protected:

  G4double RejectionFunction(G4double value) const;

private:
  G4double z;
  G4double rejection_argument1, rejection_argument2, rejection_argument3;
  G4double EnergyRatio;

  // hide assignment operator 
     G4Generator2BS & operator=(const  G4Generator2BS &right);
     G4Generator2BS(const  G4Generator2BS&);

};

#endif

