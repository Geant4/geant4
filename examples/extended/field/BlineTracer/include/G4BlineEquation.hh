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
// $Id: G4BlineEquation.hh,v 1.1 2003/11/25 09:29:46 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
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
    ~G4BlineEquation();
      // Constructor and destructor.

    void EvaluateRhsGivenB( const G4double y[],
                            const G4double B[3],
                                  G4double dydx[] ) const;
      // Given the value of the magnetic field B, this function 
      // calculates the value of the derivative dydx.

    void SetBackwardDirectionOfIntegration(G4bool abool);  

  private:

    G4bool backward_direction;
    G4double direction;
};

#endif 
