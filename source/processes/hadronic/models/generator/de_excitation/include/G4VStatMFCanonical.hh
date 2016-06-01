#ifndef G4VStatMFCanonical_h
#define G4VStatMFCanonical_h 1

#include <rw/tvordvec.h>

#include "G4Fragment.hh"
#include "G4StatMFFragment.hh"
#include "G4StatMFParameters.hh"
#include "Randomize.hh"

class G4VStatMFCanonical 
{

public:
  G4VStatMFCanonical() {};

  virtual ~G4VStatMFCanonical() {};
  
private:

  // copy constructor
  G4VStatMFCanonical(const G4VStatMFCanonical & right) {};
  
  // operators
  G4VStatMFCanonical & operator=(const G4VStatMFCanonical & right);
  G4bool operator==(const G4VStatMFCanonical & right);
  G4bool operator!=(const G4VStatMFCanonical & right);

public:

  // Choice of fragment atomic numbers and charges
  virtual void ChooseAandZ(const G4Fragment & theFragment) = 0;

  G4double GetMeanMultiplicity(void) const {return MeanMultiplicity;}

  G4double GetMeanTemperature(void) const { return MeanTemperature; }
  
  G4double GetMeanEntropy(void) const { return MeanEntropy; }

  G4int GetMultiplicity(void) const { return Multiplicity; }

  G4int GetFragmentA(const G4int & i) const 
    {
      if (i < FragmentsA.entries() && i >= 0) return FragmentsA(i);
      else {
	cout << "G4VStatMFCanonical::GetFragmentA: trying to get access to fragment "
	     << i << " from a total of " 
	     << FragmentsZ.entries() << " fragments" << endl;
	return -1;
      }
    }

  G4int GetFragmentZ(const G4int & i) const 
    {
      if (i < FragmentsZ.entries() && i >= 0) return FragmentsZ(i);
      else {
	cout << "G4VStatMFCanonical::GetFragmentZ: trying to get access to fragment "
	     << i << " from a total of " 
	     << FragmentsZ.entries() << " fragments" << endl;
	return -1;
      }
    }

  G4double GetFragmentInvLevelDensity(const G4int & i) const  
    {
      if (i < theChannels.length() && i >= 0) return theChannels(i)->GetInvLevelDensity();
      else {
	cout << "G4VStatMFCanonical::GetFragmentInvLevelDensity: trying to get access to channel "
	     << i << " from a total of "
	     << theChannels.length() << " channels." << endl;
	return 0;
      }
    }

  void SortFragments(void);

  G4int GetNumOfNeutrons(void) const { return NumOfNeutrons; }

  G4int GetNumOfCharged(void) const { return NumOfCharged; }

  G4int GetOrderedA(const G4int & i) const 
    {
      if (i < OrderedA.entries()) return OrderedA(i);
      else return 0;
    }

  G4int GetOrderedZ(const G4int & i) const 
    {
      if (i < OrderedZ.entries()) return OrderedZ(i);
      else return 0;
    }



  G4double Beta(const G4double & T) const ;
  G4double DBetaDT(const G4double & T) const ;


protected:

  // the possible channels
  RWTPtrOrderedVector<G4StatMFFragment> theChannels;

  // Free internal energy at temperature T = 0
  G4double FreeInternalE0;

  
  // Mean breakup multiplicity
  G4double MeanMultiplicity;
  
  // Mean channel temperature
  G4double MeanTemperature;
  
  // Mean channel entropy
  G4double MeanEntropy;


  // Multiplicity
  G4int Multiplicity;

  // Fragment Atomic Numbers 
  RWTValOrderedVector<G4int> FragmentsA;
  // Fragment Charges
  RWTValOrderedVector<G4int> FragmentsZ;

  G4int NumOfNeutrons;

  G4int NumOfCharged;

  RWTValOrderedVector<G4int> OrderedA;
  
  RWTValOrderedVector<G4int> OrderedZ;


};

#endif
