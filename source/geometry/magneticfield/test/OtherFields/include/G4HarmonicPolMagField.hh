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
//
// $Id: G4HarmonicPolMagField.hh,v 1.2 2006-06-29 18:26:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
	     
         void GetFieldValue(const G4double yTrack[3] ,
	                          G4double *B       ) const ;
			   

protected:

private:
              
} ;

#endif
