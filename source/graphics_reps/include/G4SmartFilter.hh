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
// $Id: G4SmartFilter.hh 66376 2012-12-18 09:42:59Z gcosmo $
//
// Filter with additional funcionality such as active and inverted states, 
// and filtering statistics
//
// Jane Tinslay, March 2006
//
#ifndef G4SMARTFILTER_HH
#define G4SMARTFILTER_HH

#include "G4VFilter.hh"

template <typename T>
class G4SmartFilter : public G4VFilter<T> {

public: // With description

  // Construct with filter name
  G4SmartFilter(const G4String& name);

  virtual ~G4SmartFilter();

  // Evaluate method implemented in subclass
  virtual G4bool Evaluate(const T&) const = 0;

  // Print subclass configuration
  virtual void Print(std::ostream& ostr) const = 0;

  // Clear filter 
  virtual void Clear() = 0;

  // Filter method
  G4bool Accept(const T&) const;
  
  // Print G4SmartFilter configuration
  virtual void PrintAll(std::ostream& ostr) const;

  //Reset
  virtual void Reset();

  // Activate/deactivate filter
  void SetActive(const G4bool&);
  G4bool GetActive() const;

  // Invert filter
  void SetInvert(const G4bool&);
  G4bool GetInvert() const;

  // Set verbosity
  void SetVerbose(const G4bool&);
  G4bool GetVerbose() const;
  

private:

  // Data members
  G4bool fActive;
  G4bool fInvert;
  G4bool fVerbose;
  mutable size_t fNPassed;
  mutable size_t fNProcessed;

};

template <typename T>
G4SmartFilter<T>::G4SmartFilter(const G4String& name)
  :G4VFilter<T>(name)
  ,fActive(true)
  ,fInvert(false)
  ,fVerbose(false)
  ,fNPassed(0)
  ,fNProcessed(0)
{}

template <typename T>
G4SmartFilter<T>::~G4SmartFilter() {}

template <typename T>
G4bool 
G4SmartFilter<T>::Accept(const T& object) const
{
  if (fVerbose) {
    G4cout<<"Begin verbose printout for filter "<<G4VFilter<T>::Name()<<G4endl;
    G4cout<<"Active ? :   "<<fActive<<G4endl;
  }
  
  fNProcessed++;
  
  // Pass everything if filter is not active
  if (!fActive) {
    fNPassed++;
    return true;
  }
  
  // Do filtering
  G4bool passed = Evaluate(object);
  
  // Apply inversion if applicable
  if (fInvert) passed = !passed;
  
  if (passed) fNPassed++;
  
  if (fVerbose) {
    G4cout<<"Inverted ? : "<<fInvert<<G4endl;
    G4cout<<"Passed ?   : "<<passed<<G4endl;
    G4cout<<"End verbose printout for filter "<<G4VFilter<T>::Name()<<G4endl;
  }
  
  return passed;
}

template <typename T>
void 
G4SmartFilter<T>::PrintAll(std::ostream& ostr) const 
{
  ostr<<"Printing data for filter: "<<G4VFilter<T>::Name()<<G4endl;

  Print(ostr);

  ostr<<"Active ?   : " <<fActive<<G4endl;
  ostr<<"Inverted ? : " <<fInvert<<G4endl;
  ostr<<"#Processed : " <<fNProcessed<<G4endl;
  ostr<<"#Passed    : " <<fNPassed<<G4endl;
}

template <typename T>
void
G4SmartFilter<T>::Reset()
{
  fActive = true;
  fInvert = false;
  fNProcessed = 0;
  fNPassed = 0;
  
  // Clear subclass data
  Clear();
}

template <typename T>
void 
G4SmartFilter<T>::SetActive(const G4bool& active) 
{
  fActive = active;
}

template <typename T>
G4bool 
G4SmartFilter<T>::GetActive() const 
{
  return fActive;
}

template <typename T>
void 
G4SmartFilter<T>::SetInvert(const G4bool& invert) 
{
  fInvert = invert;
}

template <typename T>
G4bool
G4SmartFilter<T>::GetInvert() const
{
  return fInvert;
}

template <typename T>
void 
G4SmartFilter<T>::SetVerbose(const G4bool& verbose) 
{
  fVerbose = verbose;
}

template <typename T>
G4bool
G4SmartFilter<T>::GetVerbose() const
{
  return fVerbose;
}

#endif

