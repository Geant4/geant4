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
// $Id: G4BaryonConstructor.hh,v 1.3 2001-07-11 10:01:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
#ifndef G4BaryonConstructor_h
#define G4BaryonConstructor_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4BaryonConstructor
{
  //This class is a utility class for construction 
  //short lived particles

  public:
    G4BaryonConstructor();
    ~G4BaryonConstructor();
  
  public:
    void ConstructParticle();

  protected:
    void ConstructNucleons();
    void ConstructStrangeBaryons();
    void ConstructCharmBaryons();
    void ConstructBottomBaryons();
};

#endif
