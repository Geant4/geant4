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
// $Id: PreCompoundFragmentVector.hh,v 1.1 2003-10-08 12:32:12 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear Preequilibrium
// by V. Lara 

#ifndef PreCompoundFragmentVector_h
#define PreCompoundFragmentVector_h 1


#include "VPreCompoundFragment.hh"


class PreCompoundFragmentVector 
{
  typedef std::vector<VPreCompoundFragment*>  pcfvector;
public:
  PreCompoundFragmentVector();
  ~PreCompoundFragmentVector();
	
private:
  PreCompoundFragmentVector(const PreCompoundFragmentVector &right);
  const PreCompoundFragmentVector& 
  operator=(const PreCompoundFragmentVector &right);
  G4bool operator==(const PreCompoundFragmentVector &right) const;
  G4bool operator!=(const PreCompoundFragmentVector &right) const;	
    
public:

  inline void Initialize(const G4Fragment & aFragment);
  G4double CalculateProbabilities(const G4Fragment & aFragment,G4double dLevelDensity);
	
  VPreCompoundFragment * ChooseFragment(void);
		
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

#include "PreCompoundFragmentVector.icc"

#endif










