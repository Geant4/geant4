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
#include <map>
#include "G4ParticleDefinition.hh"
#include "G4Log.hh"

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
    { return 2.0*gamma*std::atan( 2.0 * (x-m0)/ gamma  ); }

  G4double BrWigInt1(const G4double x, const G4double gamma, const G4double m0) const
    { return 0.5*gamma*gamma*G4Log( (x-m0)*(x-m0)+gamma*gamma/4.0 ) + m0*BrWigInt0(x,gamma,m0); }

  G4double BrWigInv(const G4double x, const G4double gamma, const G4double m0) const
    { return 0.5*gamma*std::tan( 0.5*x/gamma )+m0; }

public:  

  typedef std::map<const G4ParticleDefinition*, G4double, std::less<const G4ParticleDefinition*> >::const_iterator minMassMapIterator; 
  typedef std::map<const G4ParticleDefinition*, G4double, std::less<const G4ParticleDefinition*> > minMassMapType;

private:

  static G4ThreadLocal minMassMapType *minMassCache_G4MT_TLS_;

};


#endif






