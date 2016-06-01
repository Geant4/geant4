#ifndef G4StatMFMacrocanonical_h
#define G4StatMFMacrocanonical_h 1

#include <rw/tvordvec.h>

#include "G4Fragment.hh"
#include "G4StatMFFragment.hh"
#include "G4StatMFParameters.hh"
#include "G4VStatMFCanonical.hh"
#include "Randomize.hh"


class G4StatMFMacrocanonical : public G4VStatMFCanonical {

public:

  // G4StatMFMacrocanonical class must be initialized with a G4Fragment.
  G4StatMFMacrocanonical(const G4Fragment & theFragment);

  // destructor
  ~G4StatMFMacrocanonical();

private:
  // default constructor
  G4StatMFMacrocanonical() {};


  // copy constructor
  G4StatMFMacrocanonical(const G4StatMFMacrocanonical &right) {};


  // operators
  G4StatMFMacrocanonical & operator=(const G4StatMFMacrocanonical & right);
  G4bool operator==(const G4StatMFMacrocanonical & right) const;
  G4bool operator!=(const G4StatMFMacrocanonical & right) const;


public:

  // Choice of fragment atomic numbers and charges.
  void ChooseAandZ(const G4Fragment &theFragment);

private:

  // Initailization method
  void Initialize(const G4Fragment & theFragment);

  //
  void CalculateTemperature(const G4Fragment & theFragment);

  // Calculates excitation energy per nucleon and summed fragment multiplicity and entropy
  void FragmentsExcitationEnergyAndEntropy(const G4Fragment & theFragment,
					   const G4double Kappa,
					   G4double & ExcitEnergyPerNucleon, 
					   G4double & TotalMultiplicity);

  // This calculates fragment charges over fragment atomic numbers
  void CalculateZARatio(const G4Fragment & theFragment, const G4double & Kappa);

  // 
  void CalculateMultiplicities(const G4Fragment & theFragment, const G4double & Kappa);

  // Calculates fragment multiplicities
  void MeanFragmentMultiplicities(const G4Fragment & theFragment, const G4double & Kappa);

  // Calculate Fragment energies at actual temperature
  void FragmentEnergies(const G4Fragment & theFragment,const G4double & Kappa);

  // Calculates summed fragments entropy
  G4double TotalFragmentsEntropy(const G4double & A, const G4double & Kappa);

  // Determines fragments multiplicities and compute total fragment multiplicity
  G4double ChooseA(const G4double A, RWTValVector<G4double> & ANumbers);

  //
  void ChooseZ(const G4int & Z, const G4double Multiplicity);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  // Chemical Potential \mu
  G4double ChemPotentialMu;

  // Chemical Potential \nu
  G4double ChemPotentialNu;


  // 
  G4double YN, YP, Y2, Y3, Y4;

};

#endif
