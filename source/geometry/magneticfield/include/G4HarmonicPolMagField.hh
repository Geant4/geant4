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
// Class describing magnetic field parametrised by harmonic polynom up to
// 3rd order. The function MagneticField(yTrack,B) calculates the magnetic field
// induction vector B for the trajectory point yTrack according to formula given
// in: M.Metcalf, Analysis of the SFM Field, OM Development Note AP-10 (revised)
// 09.01.1974
//
// 3.2.97 V.Grichine

#ifndef G4HARMONICPOLMAGFIELD_HH
#define G4HARMONICPOLMAGFIELD_HH

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4Mag_EqRhs.hh"

class G4HarmonicPolMagField : public G4MagneticField
{
public:
		       
         G4HarmonicPolMagField() ;	       
           
        ~G4HarmonicPolMagField() ;
	     
         void GetFieldValue(const G4double yTrack[] ,
	                    G4double B[]              ) const  ;
			   

protected:

private:
              
} ;

#endif
