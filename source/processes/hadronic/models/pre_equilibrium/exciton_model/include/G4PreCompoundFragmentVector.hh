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
// $Id: G4PreCompoundFragmentVector.hh,v 1.3 2006/06/29 20:58:18 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// Hadronic Process: Nuclear Preequilibrium
// by V. Lara 

#ifndef G4PreCompoundFragmentVector_h
#define G4PreCompoundFragmentVector_h 1


#include "G4VPreCompoundFragment.hh"



class G4PreCompoundFragmentVector 
{
  typedef std::vector<G4VPreCompoundFragment*>  pcfvector;
public:
  inline G4PreCompoundFragmentVector(pcfvector * avector);
  inline ~G4PreCompoundFragmentVector();
  
private:
  G4PreCompoundFragmentVector(const G4PreCompoundFragmentVector &right);
  const G4PreCompoundFragmentVector& 
  operator=(const G4PreCompoundFragmentVector &right);
  G4bool operator==(const G4PreCompoundFragmentVector &right) const;
  G4bool operator!=(const G4PreCompoundFragmentVector &right) const;	
  
public:

  inline void Initialize(const G4Fragment & aFragment);
  inline void ResetStage();
  inline void SetVector(pcfvector * avector);

  G4double CalculateProbabilities(const G4Fragment & aFragment);
	
  G4VPreCompoundFragment * ChooseFragment(void);
		
private:

  pcfvector * theChannels;

  G4double TotalEmissionProbability;

};

#include "G4PreCompoundFragmentVector.icc"

#endif

