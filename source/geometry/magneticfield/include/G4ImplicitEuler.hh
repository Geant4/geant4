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
// $Id: G4ImplicitEuler.hh,v 1.5 2002-11-29 13:48:47 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// W. Wander <wwc@mit.edu> 12/09/97

#ifndef G4IMPLICITEULER_HH
#define G4IMPLICITEULER_HH
#include "G4MagErrorStepper.hh"

class G4ImplicitEuler : public G4MagErrorStepper
{

  public:
    G4ImplicitEuler(G4Mag_EqRhs *EqRhs, G4int numberOfVariables = 6) ;
   ~G4ImplicitEuler();

    void  DumbStepper(  const G4double y[] ,
		        const G4double dydx[] ,
		              G4double h ,
			      G4double yout[] ) ;

    G4int IntegratorOrder() const { return 2 ; } ;

private: // No description

    G4double*  dydxTemp;
    G4double*  yTemp;    
    // Temporaries, created to avoid new/delete on every call
};

#endif /* G4IMPLICITEULER_HH */
