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
// $Id: G4LineCurrentMagField.hh,v 1.1 2002-03-28 13:45:58 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// Class describing line current magnetic field.
// fFieldConstant determines the coefficient in the field law.
// The line current is directed along Z axis and crosses the XY plane in the point
// (0,0) .
//
// 3.2.97 V. Grichine

#ifndef G4LINECURRENTMAGFIELD_HH
#define G4LINECURRENTMAGFIELD_HH

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4Mag_EqRhs.hh"

class G4LineCurrentMagField : public G4MagneticField
{
public:
		       
         G4LineCurrentMagField(G4double pFieldConstant) ;	       
           
        ~G4LineCurrentMagField() ;
	     
         void MagneticField(const G4double yTrack[] ,
	                    G4double B[]              ) ;
			   

protected:

private:
          G4double fFieldConstant ;
              
} ;

#endif
