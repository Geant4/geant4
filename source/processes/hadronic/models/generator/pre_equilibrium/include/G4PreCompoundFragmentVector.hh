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
// $Id: G4PreCompoundFragmentVector.hh,v 1.5 2002/06/06 17:03:57 larazb Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear Preequilibrium
// by V. Lara 

#ifndef G4PreCompoundFragmentVector_h
#define G4PreCompoundFragmentVector_h 1


#include "G4VPreCompoundFragment.hh"


class G4PreCompoundFragmentVector 
{
  typedef G4std::vector<G4VPreCompoundFragment*>  pcfvector;
public:
  G4PreCompoundFragmentVector();
  ~G4PreCompoundFragmentVector();
	
private:
  G4PreCompoundFragmentVector(const G4PreCompoundFragmentVector &right);
  const G4PreCompoundFragmentVector& 
  operator=(const G4PreCompoundFragmentVector &right);
  G4bool operator==(const G4PreCompoundFragmentVector &right) const;
  G4bool operator!=(const G4PreCompoundFragmentVector &right) const;	
    
public:

  inline void Initialize(const G4Fragment & aFragment);
  G4double CalculateProbabilities(const G4Fragment & aFragment);
	
  G4VPreCompoundFragment * ChooseFragment(void);
		
private:

  pcfvector theChannels;

  G4double TotalEmissionProbability;

  struct DeleteFragment 
  {
    template<typename T>
    void operator()(const T* ptr) const
    {
      delete ptr;
    }
  };

};

#include "G4PreCompoundFragmentVector.icc"

#endif










