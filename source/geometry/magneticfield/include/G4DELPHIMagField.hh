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
// Class describing the DELPHI magnetic field. This axial symmetry
// field mainly directed along Z axis. The function MagneticField(yTrack,B)
// calculates the magnetic induction vector B in point corresponding to
// yTrack according to parametrization given in:
// P.Billoir, Precise tracking in a quasi-honogeneous magnetic field, DELPHI 
// 87-6 PROG 65, 3 February 1987
//
// 3.2.97 V. Grichine 

#ifndef G4DELPHIMAGFIELD_HH
#define G4DELPHIMAGFIELD_HH

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4Mag_EqRhs.hh"

class G4DELPHIMagField : public G4MagneticField
{
public:
		       
         G4DELPHIMagField() ;	       
           
        ~G4DELPHIMagField() ;
	     
         void GetFieldValue(const G4double yTrack[] ,
	                    G4double B[]              ) const ;
			   

protected:

private:
              
} ;

#endif
