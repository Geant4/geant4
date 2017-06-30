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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//      File name:     G4NuclearPolarization
//
//      Author:        Jason Detwiler (jasondet@gmail.com)
// 
//      Creation date: Aug 2015 
//
//      Description:   
//      Stores the statistical tensor describing the nuclear polarization 
//      (see Alder and Winther, "Electromagnetic Excitation" (1975),
//      Appendix F). 
//
//      V.Ivanchenko left only polarization tensor and access methods
//                   in this class, also add operators
//                   this allows future implemention of polarized 
//                   hadronic models with polarized initial and 
//                   final states
//
// -------------------------------------------------------------------

#ifndef G4NUCLEARPOLARIZATION_HH
#define G4NUCLEARPOLARIZATION_HH

#include "globals.hh"
#include <vector>

class G4NuclearPolarization
{
public:

  G4NuclearPolarization();

  ~G4NuclearPolarization();

  inline void Unpolarize() 
  { 
    fPolarization.resize(1); 
    fPolarization[0].assign(1, 1.); 
  }
  
  inline void SetPolarization(std::vector< std::vector<G4complex> >& p) 
  { 
    fPolarization = p; 
  }

  inline std::vector< std::vector<G4complex> >& GetPolarization() 
  { 
    return fPolarization; 
  }

  // ============= OPERATORS ==================

  inline G4NuclearPolarization & operator=(const G4NuclearPolarization &right)
  {
    if (this != &right) { fPolarization = right.fPolarization; }
    return *this;
  }
    
  inline G4NuclearPolarization(const G4NuclearPolarization &right )
  { 
    *this = right; 
  }

  G4bool operator==(const G4NuclearPolarization &right) const;
  G4bool operator!=(const G4NuclearPolarization &right) const;

  friend std::ostream& operator<<(std::ostream&, const G4NuclearPolarization&);

private:

  std::vector< std::vector<G4complex> > fPolarization;
};

#endif
