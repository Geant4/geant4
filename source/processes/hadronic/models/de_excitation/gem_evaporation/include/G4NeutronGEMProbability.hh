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
// $Id: G4NeutronGEMProbability.hh,v 1.5 2009-09-15 12:54:16 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999) 
//
// J.M. Quesada (July 2009) alpha's and beta's  values according to Furihata's paper 
// (based on notes added on proof in Dostrovskii's paper)


#ifndef G4NeutronGEMProbability_h
#define G4NeutronGEMProbability_h 1


#include "G4GEMProbability.hh"


class G4NeutronGEMProbability : public G4GEMProbability
{
public:
  // Only available constructor
  G4NeutronGEMProbability();
  
  ~G4NeutronGEMProbability() {}
private:  
  // Copy constructor
  G4NeutronGEMProbability(const G4NeutronGEMProbability &right);
  
  const G4NeutronGEMProbability & operator=(const G4NeutronGEMProbability &right);
  G4bool operator==(const G4NeutronGEMProbability &right) const;
  G4bool operator!=(const G4NeutronGEMProbability &right) const;
  
  
private:
  
  virtual G4double CalcAlphaParam(const G4Fragment & fragment) const 
  {
    //JMQ 190709 values according to Furihata's paper (based on notes added on proof in Dostrovskii's paper)
    //            return 0.76+2.2/std::pow(static_cast<G4double>(fragment.GetA()-GetA()),1.0/3.0);
    return 0.76+1.93/std::pow(static_cast<G4double>(fragment.GetA()-GetA()),1.0/3.0);
  }
  
  virtual G4double CalcBetaParam(const G4Fragment & fragment) const 
  {
    //JMQ 190709 values according to Furihata's paper (based on notes added on proof in Dostrovskii's paper)
    //            return (2.12/std::pow(static_cast<G4double>(fragment.GetA()-GetA()),2.0/3.0)-0.05)*MeV/
    return (1.66/std::pow(static_cast<G4double>(fragment.GetA()-GetA()),2.0/3.0)-0.05)*MeV/
      CalcAlphaParam(fragment);
  }
  
  // Excitation energy levels 
  std::vector<G4double> ExcitEnergies;
  // Spin of excitation energy levels 
  std::vector<G4double> ExcitSpins;
  
  std::vector<G4double> ExcitLifetimes;
  
};


#endif
