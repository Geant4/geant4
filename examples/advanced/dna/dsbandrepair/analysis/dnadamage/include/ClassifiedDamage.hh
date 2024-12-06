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
//
/// \file ClassifiedDamage.hh
/// \brief Definition of the ClassifiedDamage class

#ifndef CLASSIFIEDDAMAGE_HH
#define CLASSIFIEDDAMAGE_HH

#include "Damage.hh"

#include <vector>

/// \brief defines a classified DNA damage
class ClassifiedDamage
{
public:

// DNA Classified Damage type
//   Defines the type of damage ollowing SDD formalism 
  enum ClassifiedDamageType{
    fNA = -1,
    fSSB = 0,
    fDSB = 1
  };

  /// \brief constructor
  ClassifiedDamage();
  /// \brief destructor
  ~ClassifiedDamage() = default;

  // Compute the type of classified damage i.e. SSB or DSB
  void ComputeType();
  ClassifiedDamageType GetClassifiedDamageType() const {return fType;};

  const int GetNumDamage() const {return fDamage.size();};
  // Add a damage to the lst of damage
  void AddDamage(Damage);

  // Compute the position in terms of bp of the classified damage
  void ComputeBp();
  unsigned long int GetBpBegin() const{return fBp_begin;};
  unsigned long int GetBpEnd() const{return fBp_end;};
  unsigned long int GetBpBarycenter() const{return fBp_barycenter;};

  // TODO
  // Compute the coordinates of the classified damage
  void ComputePosition();
  /*
  Point GetPosBegin(){return fPos_begin;};
  Point GetPosEnd(){return fPos_end;};
  Point GetPosBarycenter(){return fPos_barycenter;};
  */

  // Compute the complexity of the classified damage
  void ComputeComplexity();
  int GetComplexity() const{return fComplexity;};

  // Reset the classified damage
  void Reset();

  // Tell if base damages have to be taken into account
  void SetIncludeBase(bool pVal) {fIncludeBase =  pVal;};
  bool GetIncludeBase() const {return fIncludeBase;};

  // Le Tuan Anh:
  bool GetIsThereDirectComponentContribution() const 
  {return fIsThereDirectContribution;} // Return true if there is at least 1 direct SB in this cluster
  bool GetIsThereIndirectComponentContribution() const 
  {return fIsThereIndirectContribution;} // Return true if there is at least 1 indirect SB in this cluster

private:

  // CLASSIFIED DAMAGE MEMBERS

  std::vector<Damage>   fDamage;              // List of damage
  ClassifiedDamageType  fType;                // SSB or DSB or other?
  unsigned long int     fBp_begin{0};            // Position
  unsigned long int     fBp_end{0};
  unsigned long int     fBp_barycenter{0};
  int                   fComplexity{0};          // Complexity
  bool                  fIncludeBase{false};        // Base inclusion in the complexity
  bool fIsThereDirectContribution = false; // check if Direct damage appears in cluster a not?
  bool fIsThereIndirectContribution = false; // check if Indirect damage appears in cluster a not?

  };

#endif // CLASSIFIEDDAMAGE_HH