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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4SampleResonance.hh,v 1.4 2001/10/04 20:00:33 hpw Exp $
// GEANT4 tag $Name: geant4-04-00 $
//

#ifndef G4SampleResonance_h
#define G4SampleResonance_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4SampleResonance ----------------
//             by Henning Weber, March 2001.
//      helper class for sampling resonance masses
// ------------------------------------------------------------

#include "globals.hh"
#include "g4std/map"
#include "G4ParticleDefinition.hh"


class G4SampleResonance
{
public:

  G4double GetMinimumMass(const G4ParticleDefinition* p) const;
  G4double SampleMass(const G4double poleMass, 
		      const G4double gamma,
		      const G4double minMass,
		      const G4double maxMass) const;
  G4double SampleMass(const G4ParticleDefinition* p, const G4double maxMass) const;

private:  

  G4double BrWigInt0(const G4double x, const G4double gamma, const G4double m0) const
    { return 2.0*gamma*atan( 2.0 * (x-m0)/ gamma  ); }

  G4double BrWigInt1(const G4double x, const G4double gamma, const G4double m0) const
    { return 0.5*gamma*gamma*log( (x-m0)*(x-m0)+gamma*gamma/4.0 ) + m0*BrWigInt0(x,gamma,m0); }

  G4double BrWigInv(const G4double x, const G4double gamma, const G4double m0) const
    { return 0.5*gamma*tan( 0.5*x/gamma )+m0; }

public:  

  typedef G4std::map<const G4ParticleDefinition*, G4double, G4std::less<const G4ParticleDefinition*> >::const_iterator minMassMapIterator; 
  typedef G4std::map<const G4ParticleDefinition*, G4double, G4std::less<const G4ParticleDefinition*> > minMassMapType;

private:

  static minMassMapType minMassCache;

};


#endif






