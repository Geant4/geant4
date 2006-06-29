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
#ifndef G4RadioactiveDecayRate_h
#define G4RadioactiveDecayRate_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              RadioactiveDecayRate.hh
//
// Version:             0.b.4
// Date:                14/04/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// 13 April 2000, F Lei, DERA UK
// 0.b.4 release. No change to this file     
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "G4ios.hh"
#include "globals.hh"

#include <vector>
////////////////////////////////////////////////////////////////////////////////
//
class G4RadioactiveDecayRate
{
  // class description
  // This class contains the coefficient and decay times of the
  // progeny of a give isotope (A,Z,E).
  //
  // class description - end
public:
  // Constructors
  G4RadioactiveDecayRate();

  //  Destructor
  virtual ~G4RadioactiveDecayRate();

public:
  //  copy constructor and assignment operator
  G4RadioactiveDecayRate(const G4RadioactiveDecayRate &);
  G4RadioactiveDecayRate & operator=(const G4RadioactiveDecayRate &);

public:
  // equality operators
  G4int operator==(const G4RadioactiveDecayRate &right) const
    {return (this == &right);};
  G4int operator!=(const G4RadioactiveDecayRate &right) const
    {return (this != &right);};

  // less-than operator is defined for G4DecayTable
  // G4int operator<(const G4RadioactiveDecayRate &right) const;

public: // with description
  //
  // the inline member functions are self explanatory.
  //
  inline G4int GetZ() const {return Z;}
  inline G4int GetA() const { return A;}
  inline G4double GetE() const { return E;}
  inline G4int GetGeneration() const { return generation;}
  inline std::vector<G4double> GetDecayRateC() const
   {  return decayRateC; }
  inline std::vector<G4double> GetTaos() const {  return taos; }

  inline void SetZ(G4int value) {Z = value;}
  inline void SetA(G4int value) {A = value;}
  inline void SetE(G4double value) {E = value;}
  inline void SetGeneration(G4int value) {generation = value;}
  inline void SetDecayRateC(std::vector<G4double> value)
    {decayRateC = value;}
  inline void SetTaos(std::vector<G4double> value) {taos = value;}

protected:

  G4int                   Z;
  G4int                   A;
  G4double                E;
  G4int                   generation;
  std::vector<G4double> decayRateC;
  std::vector<G4double> taos;

public:

  inline void  SetVerboseLevel(G4int value)
    { verboseLevel = value; }
  inline G4int GetVerboseLevel() const
    { return verboseLevel; }
  void  DumpInfo();

private:
  G4int verboseLevel;
  // control flag for output message
  // G4int verboseLevel;
  //  0: Silent
  //  1: Warning message
  //  2: More

};
#endif
