#ifndef G4StatMFMicrocanonical_h
#define G4StatMFMicrocanonical_h 1

//#include <rw/tvvector.h>
#include <rw/tvordvec.h>

#include "G4Fragment.hh"
#include "G4StatMFFragment.hh"
#include "G4StatMFParameters.hh"
#include "G4VStatMFCanonical.hh"
#include "Randomize.hh"

//class G4StatMF1DVector : public RWTValVector<G4int> {
//public:
//  G4StatMF1DVector() {};
//  G4StatMF1DVector(G4int n):RWTValVector<G4int>(n) {};
//};




class G4StatMFMicrocanonical : public G4VStatMFCanonical {

public:

  // G4StatMFMicrocanonical class must be initialized with a G4Fragment.
  G4StatMFMicrocanonical(const G4Fragment & theFragment);

  // destructor
  ~G4StatMFMicrocanonical();

private:
  // default constructor
  G4StatMFMicrocanonical() {};


  // copy constructor
  G4StatMFMicrocanonical(const G4StatMFMicrocanonical &right) {};


  // operators
  G4StatMFMicrocanonical & operator=(const G4StatMFMicrocanonical & right);
  G4bool operator==(const G4StatMFMicrocanonical & right) const;
  G4bool operator!=(const G4StatMFMicrocanonical & right) const;


public:

  // Choice of fragment atomic numbers and charges.
  void ChooseAandZ(const G4Fragment &theFragment);


private:

  // Initailization method
  void Initialize(const G4Fragment & theFragment);

  // Calculate Entropy of Compound Nucleus
  G4double CalcEntropyOfCompoundNucleus(const G4Fragment & theFragment, G4double & TConf);

  G4bool DistributeNucleonsBetweenFragments(const G4int & k, G4int * ANumbers);

  G4double CalcFragmentsConfigProbability(const G4Fragment & theFragment, const G4int & M,
					  const G4int * ANumbers, const G4double & SCompound);


  G4double CalcFreeInternalEnergy(const G4Fragment & theFragment, const G4double & T);


  // Gives fragments charges
  void ChooseZ(const G4Fragment & theFragment, const G4int & FragmentMultiplicity);




  // -----------
  G4double CalcEnergyConfiguration(const G4double A, const G4double Z, const G4int M, 
				   G4double * ECOLA, G4double * EA, const G4int * Anumbers, 
				   const G4double T);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  // Statistical weights 
  G4double W, WW2, WW3, WW4;
  RWTValOrderedVector<G4double> W2, W3, W4;

  // Number of configurations for breakups with multiplicities 2, 3 and 4
  G4int M2, M3, M4;

  // Atomic numbers of fragments for each configuration with multiplicities 2, 3 and 4
  //  RWTValOrderedVector<G4StatMF1DVector> ANum2;
  //  RWTValOrderedVector< RWTValVector<G4int> > ANum2;
  RWTValOrderedVector< G4int* > ANum2;
  //  RWTValOrderedVector<G4StatMF1DVector> ANum3;
  //  RWTValOrderedVector< RWTValVector<G4int> > ANum3;
  RWTValOrderedVector< G4int* > ANum3;
  //  RWTValOrderedVector<G4StatMF1DVector> ANum4;
  //  RWTValOrderedVector< RWTValVector<G4int> > ANum4;
  RWTValOrderedVector< G4int* > ANum4;

  // Statistical weight of compound nucleus
  G4double WCompoundNucleus;


};

#endif
