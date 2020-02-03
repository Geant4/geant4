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
/// \file field/BlineTracer/include/G4BlineEquation.hh
/// \brief Definition of the G4BlineEquation class
//
//
//
// 
// --------------------------------------------------------------------
//
// G4BlineEquation
//
// Class description:
//
// This class defines the equation of motion needed to trace magnetic
// field lines in the simulation. 

// --------------------------------------------------------------------
// Author: Laurent Desorgher (desorgher@phim.unibe.ch)
//         Created - 2003-10-06
// --------------------------------------------------------------------
#ifndef G4BlineEquation_h
#define G4BlineEquation_h 1

#include "G4Mag_EqRhs.hh"
#include "G4MagneticField.hh"

class G4BlineEquation : public G4Mag_EqRhs
{
  public:  // with description

    G4BlineEquation( G4MagneticField* MagField );
    virtual ~G4BlineEquation();
      // Constructor and destructor.

    virtual void EvaluateRhsGivenB( const G4double y[],
                            const G4double B[3],
                                  G4double dydx[] ) const;
      // Given the value of the magnetic field B, this function 
      // calculates the value of the derivative dydx.

    void SetBackwardDirectionOfIntegration(G4bool abool);  

  private:

    G4bool fBackward_direction;
    G4double fDirection;
};

#endif 
