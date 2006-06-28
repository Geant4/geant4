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
// -------------------------------------------------------------------
// $Id: MicrobeamEMField.hh,v 1.4 2006-06-28 13:43:01 gunter Exp $
// -------------------------------------------------------------------

#ifndef MicrobeamEMField_h
#define MicrobeamEMField_h 1

#include "globals.hh"
#include "G4ElectroMagneticField.hh"

class MicrobeamEMField
#ifndef STANDALONE
 : public G4ElectroMagneticField
#endif

{
  
public:
  MicrobeamEMField();
  void  GetFieldValue( const  double Point[4], double *Bfield ) const;
		       
  G4bool DoesFieldChangeEnergy() const {return true;}

};

#endif
