// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VBoson.hh,v 1.2 1999-12-15 14:50:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//  Revised, Hisaya Kurashige, 15 Dec 1996
//  Revised, Hisaya Kurashige, 24 Feb 1997
// ------------------------------------------------------------

#ifndef G4VBoson_h
#define G4VBoson_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4Material.hh"

#include "G4PhysicsLogVector.hh"
#include "G4ParticleWithCuts.hh"

class G4VBoson : public G4ParticleWithCuts
{
  //  A virtual class for Bosons particles. It defines
  //  public methods which describe the behavior of a
  //  boson.

  private:

      const G4VBoson & operator=(const G4VBoson &right);

  public:

      G4VBoson(const G4String&  aName,  
               G4double         mass,     
               G4double         width,
               G4double         charge,   
               G4int            iSpin,
               G4int            iParity,
               G4int            iConjugation,
               G4int            iIsospin,   
               G4int            iIsospinZ, 
               G4int            gParity,
               const G4String&  pType,
               G4int            lepton,
               G4int            baryon,
               G4int            encoding,
               G4bool           stable,
               G4double         lifetime,
               G4DecayTable     *decaytable)
	: G4ParticleWithCuts(aName, mass, width, charge, iSpin, iParity,
                               iConjugation, iIsospin, iIsospinZ, gParity,
                               pType, lepton, baryon, encoding, stable,
                               lifetime, decaytable) {};

      virtual ~G4VBoson() {};

      G4int operator==(const G4VBoson &right) const;
      G4int operator!=(const G4VBoson &right) const;

};

#endif
