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
// $Id: G4CutsPerMaterialWarning.hh,v 1.1 2001-11-07 22:38:47 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 7 Nov 2001   MGP        Created
//
// -------------------------------------------------------------------
// Class description:
// Functor class to print a warning if the deprecated, unsupported feature
// of having different range cuts for different materials is used
// Will be removed when properly designed "range cuts by region" 
// will be available in Geant4
// Further documentation available from http://www.ge.infn.it/geant4/lowE/

// -------------------------------------------------------------------

#ifndef G4CUTSPERMATERIALWARNING_HH
#define G4CUTSPERMATERIALWARNING_HH 1

#include "globals.hh"

class G4ParticleDefinition;
class G4Material;

class G4CutsPerMaterialWarning
{ 
public:

  G4CutsPerMaterialWarning() { }

  ~G4CutsPerMaterialWarning() { }
 
  void PrintWarning(const G4ParticleDefinition* particle) const;

private:

};
 
#endif
 










