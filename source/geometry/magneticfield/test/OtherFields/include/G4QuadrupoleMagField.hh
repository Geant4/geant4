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
// $Id: G4QuadrupoleMagField.hh,v 1.2 2002-06-25 13:12:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//


// Class for creation of quadrupole magnetic field
// fGradient - is gardient value for  quadrupole magnetic lense. Then the magnetic
// field components are:
// Bx = B[0] = fGradient*X ,
// By = B[1] = fGradient*Y ,
// Bz = B[2] = 0 .
// Here X,Y,Z are the coordinates of a space point of interest.
//
// 3.2.97 V.Grichine

#ifndef G4QUADRUPOLEMAGFIELD_HH
#define G4QUADRUPOLEMAGFIELD_HH

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4Mag_EqRhs.hh"

class G4QuadrupoleMagField : public G4MagneticField
{
public:
		       
	 G4QuadrupoleMagField(G4double pGradient) ;	       
           
        ~G4QuadrupoleMagField() ;
	     
         void MagneticField(const G4double yTrack[] ,
	                    G4double B[]              ) ;
			   

protected:

private:
          G4double fGradient ;
              
} ;

#endif
