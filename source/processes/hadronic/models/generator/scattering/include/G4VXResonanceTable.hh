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
//      GEANT4 Class file
//
//
//      File name:     G4VXResonanceTable
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 4 June 1999
//
//      Modifications: 
//      
// Hadron Kinetic Model
// Base class for cross section tables
//
// -------------------------------------------------------------------

#ifndef G4VXRESONANCETABLE_HH
#define G4VXRESONANCETABLE_HH

#include "globals.hh"
#include "G4PhysicsVector.hh"

class G4VXResonanceTable 
{

public:

  G4VXResonanceTable() { }

  virtual ~G4VXResonanceTable() { }

  virtual G4PhysicsVector* CrossSectionTable() const = 0;


protected:

private:  

};

#endif

