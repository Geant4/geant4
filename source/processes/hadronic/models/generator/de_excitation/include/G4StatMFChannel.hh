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
// $Id: G4StatMFChannel.hh,v 1.9 2002/12/12 19:17:12 gunter Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4StatMFChannel_h
#define G4StatMFChannel_h 1

#include "g4std/deque"

#include "G4StatMFParameters.hh"
#include "G4StatMFFragment.hh"


class G4StatMFChannel {

public:
    // Default Constructor
    G4StatMFChannel() : 
	_NumOfNeutralFragments(0), 
	_NumOfChargedFragments(0)
	{}


    // Destructor
    ~G4StatMFChannel() { 
	if (!_theFragments.empty()) {
	  G4std::for_each(_theFragments.begin(),_theFragments.end(),
			  DeleteFragment());
	}
    }

private:

    // Copy constructor
    G4StatMFChannel(const G4StatMFChannel & right);

    // operators
    G4StatMFChannel & operator=(const G4StatMFChannel & right);

public:
    G4bool operator==(const G4StatMFChannel & right) const;
    G4bool operator!=(const G4StatMFChannel & right) const;
	
public:

    void CreateFragment(const G4double A, const G4double Z);
	
    G4int GetMultiplicity(void) { return _theFragments.size();}
	
    // Return false if there is some unphysical fragment
    G4bool CheckFragments(void);


    G4double GetFragmentsCoulombEnergy(void);


    G4double GetFragmentsEnergy(const G4double T) const;
	
	
    G4FragmentVector * GetFragments(const G4double anA, const G4double anZ, const G4double T);
	

private:


    // This method calculates asymptotic fragments momenta.
    void CoulombImpulse(const G4double anA, const G4double anZ, const G4double T);
	
    void PlaceFragments(const G4double anA);

    void SolveEqOfMotion(const G4double anA, const G4double anZ, const G4double T);


    // Calculates fragments momentum components at the breakup instant.
    // Fragment kinetic energies will be calculated according to the
    // Boltzamann distribution at given temperature.
    void FragmentsMomenta(const G4int NF, const G4int idx, const G4double T);	


    // Samples a isotropic random vectorwith a magnitud given by Magnitude.
    // By default Magnitude = 1
    G4ThreeVector IsotropicVector(const G4double Magnitude = 1.0);


    // Rotates a 3-vector P to close momentum triangle Pa + V + P = 0
    G4ThreeVector RotateMomentum(G4ThreeVector Pa, G4ThreeVector V, 
				 G4ThreeVector P);

private:

    G4std::deque<G4StatMFFragment*> _theFragments;

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







