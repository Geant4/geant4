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
// $Id: G4StatMFMicroCanonical.hh,v 1.1 2003/08/26 18:47:25 lara Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4StatMFMicroCanonical_h
#define G4StatMFMicroCanonical_h 1

#include <vector>

#include "G4VStatMFEnsemble.hh"
#include "G4StatMFMicroPartition.hh"
#include "G4StatMFMicroManager.hh"
#include "G4StatMFParameters.hh"
#include "G4StatMFChannel.hh"

#include "G4Fragment.hh"
#include "Randomize.hh"


class G4StatMFMicroCanonical : public G4VStatMFEnsemble {

public:

    // G4StatMFMicroCanonical class must be initialized with a G4Fragment.
    G4StatMFMicroCanonical(const G4Fragment & theFragment);

    // destructor
    ~G4StatMFMicroCanonical();

private:
    // default constructor
    G4StatMFMicroCanonical() {};


    // copy constructor
    G4StatMFMicroCanonical(const G4StatMFMicroCanonical &right);


    // operators
    G4StatMFMicroCanonical & operator=(const G4StatMFMicroCanonical & right);
    G4bool operator==(const G4StatMFMicroCanonical & right) const;
    G4bool operator!=(const G4StatMFMicroCanonical & right) const;


public:

    // Choice of fragment atomic numbers and charges.
    G4StatMFChannel * ChooseAandZ(const G4Fragment & theFragment);
	
    enum {MaxAllowedMultiplicity = 4};

private:

    // Initailization method
    void Initialize(const G4Fragment & theFragment);

    // Calculate Entropy of Compound Nucleus
    G4double CalcEntropyOfCompoundNucleus(const G4Fragment & theFragment, G4double & TConf);


    G4double CalcFreeInternalEnergy(const G4Fragment & theFragment, const G4double T);


    G4double CalcInvLevelDensity(const G4int anA);
	
	
// Data members
private:
	
    // This is a vector of partitions managers for partitions of different 
    // multiplicities:
    
    std::vector<G4StatMFMicroManager*> _ThePartitionManagerVector;
	


    // Statistical weight of compound nucleus
    G4double _WCompoundNucleus;


  struct DeleteFragment 
  {
    template<typename T>
    void operator()(const T* ptr) const
    {
      delete ptr;
    }
  };

  class SumProbabilities : public std::binary_function<G4double,G4double,G4double>
  {
  public:
    SumProbabilities() : total(0.0) {}
    G4double operator() (G4double& /* probSoFar*/, G4StatMFMicroManager*& manager)
    { 
      total += manager->GetProbability();
      return total;
    }
    
    G4double GetTotal() { return total; }
  public:
    G4double total;
    
  };


};

#endif
