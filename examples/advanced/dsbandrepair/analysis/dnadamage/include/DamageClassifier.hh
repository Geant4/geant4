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
/// \file DamageClassifier.hh
/// \brief Definition of the DamageClassifier class

#ifndef DAMAGECLASSIFIER_HH
#define DAMAGECLASSIFIER_HH

#include <map>
#include <vector>
#include "Damage.hh"
#include "ClassifiedDamage.hh"

/// \brief defines the tool to make cluster of damage
class DamageClassifier
{
public:
  /// \brief constructor
  DamageClassifier() = default;
  /// \brief destructor
  ~DamageClassifier() = default;

  // Make a cluster of damage
  std::vector<ClassifiedDamage> MakeCluster(std::vector<Damage>&,unsigned int pDSBLength,bool pBase);

  // Return the number of SSB inside a list of classified damage
  unsigned int GetNumSSB(const std::vector<ClassifiedDamage>&) const;
  // Return the number of DSB (simple + complex) inside a list of classified damage
  unsigned int GetNumDSB(const std::vector<ClassifiedDamage>&) const;
  // Return the number of complex DSB inside a list of classified damage
  unsigned int GetNumComplexDSB(const std::vector<ClassifiedDamage>&) const;

  // Le Tuan Anh: Return the number of DSB (simple + complex) with the contribution 
  //of at least 1 direct damage inside a list of classified damage
  unsigned int GetNumDSBwithDirectDamage(const std::vector<ClassifiedDamage>&) const; 
  // Le Tuan Anh: Return the number of DSB (simple + complex) with the contribution 
  //of at least 1 indirect damage inside a list of classified damage
  unsigned int GetNumDSBwithIndirectDamage(const std::vector<ClassifiedDamage>&) const; 
  // Le Tuan Anh: Return the number of DSB (simple + complex) with the contribution of 
  //both direct and indirect damage inside a list of classified damage
  unsigned int GetNumDSBwithBothDirectIndirectDamage(const std::vector<ClassifiedDamage>&) const; 

// Utils to sort a vector of damage in maps to easily have access by event or by chromosome ID
// First key is event, second is chromosome ID
std::map<unsigned int,std::map<unsigned int,std::vector<Damage>>> 
SortDamageByEvent(const std::vector<Damage>&);
// First key is chromosome ID, second is event
std::map<unsigned int,std::map<unsigned int,std::vector<Damage>>> 
SortDamageByChromo(const std::vector<Damage>&);

};

#endif //
