// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BarionConstructor.hh,v 1.1 1999-01-07 16:09:53 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
#ifndef G4BarionConstructor_h
#define G4BarionConstructor_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4BarionConstructor
{
  //This class is a utility class for construction 
  //short lived particles

  public:
    G4BarionConstructor();
    ~G4BarionConstructor();
  
  public:
    void ConstructParticle();

  protected:
    void ConstructNucleons();
    void ConstructStrangeBarions();
    void ConstructCharmBarions();
    void ConstructBottomBarions();
};

#endif
