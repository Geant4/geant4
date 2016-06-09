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
// $Id: G4FermiConfiguration.hh,v 1.2 2003/11/03 09:38:59 lara Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4FermiConfiguration_h
#define G4FermiConfiguration_h 1

#include <deque>

#include "globals.hh"
#include "Randomize.hh"
#include "G4VFermiFragment.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Fragment.hh"


class G4FermiConfiguration 
{
public:
  // Constructors
  inline G4FermiConfiguration();
  inline ~G4FermiConfiguration();
  inline G4FermiConfiguration(const std::vector<const G4VFermiFragment*>&);
  inline G4FermiConfiguration(const G4FermiConfiguration &);
  
  // Operators
  inline const G4FermiConfiguration & operator=(const G4FermiConfiguration &);
  inline G4bool operator==(const G4FermiConfiguration &) const;
  inline G4bool operator!=(const G4FermiConfiguration &) const;
  
  
  inline void SetConfiguration(const std::vector<const G4VFermiFragment*>&);

  G4double DecayProbability(const G4int A, const G4double TotalE);

  G4FragmentVector * GetFragments(const G4Fragment & theNucleus);
  
  inline G4int GetNumberOfFragments() const;

private:

  G4double CoulombBarrier(void);

  G4ParticleMomentum IsotropicVector(const G4double Magnitude = 1.0);


private:
  // Kappa = V/V_0 it is used in calculation of Coulomb energy
  static const G4double Kappa;
  
  // Nuclear radius r0 (is a model parameter)
  static const G4double r0;
  
  std::vector<const G4VFermiFragment*> Configuration;

  struct DeleteFragment 
  {
    template<typename T>
    void operator()(const T* ptr) const
    {
      delete ptr;
    }
  };

};

#include "G4FermiConfiguration.icc"

#endif


