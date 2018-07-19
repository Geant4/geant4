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
// $Id: G4StatMFChannel.hh 107060 2017-11-01 15:00:04Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4StatMFChannel_h
#define G4StatMFChannel_h 1

#include <deque>

#include "G4StatMFParameters.hh"
#include "G4StatMFFragment.hh"


class G4StatMFChannel {

public:
    // Default Constructor
    G4StatMFChannel();

    // Destructor
    ~G4StatMFChannel();

private:

    // Copy constructor
    G4StatMFChannel(const G4StatMFChannel & right);

    // operators
    G4StatMFChannel & operator=(const G4StatMFChannel & right);

    G4bool operator==(const G4StatMFChannel & right) const;
    G4bool operator!=(const G4StatMFChannel & right) const;
	
public:

    void CreateFragment(G4int A, G4int Z);
	
    inline size_t GetMultiplicity(void) { return _theFragments.size();}
	
    // Return false if there is some unphysical fragment
    G4bool CheckFragments(void);

    G4double GetFragmentsCoulombEnergy(void);

    G4double GetFragmentsEnergy(G4double T) const;
	
    G4FragmentVector * GetFragments(G4int anA, G4int anZ, G4double T);
	
private:

    // This method calculates asymptotic fragments momenta.
    void CoulombImpulse(G4int anA, G4int anZ, G4double T);
	
    void PlaceFragments(G4int anA);

    void SolveEqOfMotion(G4int anA, G4int anZ, G4double T);

    // Calculates fragments momentum components at the breakup instant.
    // Fragment kinetic energies will be calculated according to the
    // Boltzamann distribution at given temperature.
    void FragmentsMomenta(G4int NF, G4int idx, G4double T);	

    // Rotates a 3-vector P to close momentum triangle Pa + V + P = 0
    G4ThreeVector RotateMomentum(G4ThreeVector Pa, G4ThreeVector V, 
				 G4ThreeVector P);

private:

    std::deque<G4StatMFFragment*> _theFragments;

    G4int _NumOfNeutralFragments;
	
    G4int _NumOfChargedFragments;

  struct DeleteFragment 
  {
    template<typename T>
    void operator()(const T* ptr) const
    {
      delete ptr;
    }
  };

};

#endif







