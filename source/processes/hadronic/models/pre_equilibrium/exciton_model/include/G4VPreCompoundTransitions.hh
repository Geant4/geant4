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
//J. M. Quesada (May 08). New virtual classes have been added 
// JMQ (06 September 2008) Also external choices have been added for:
//                      - "never go back"  hipothesis (useNGB=true) 
//                      - CEM transition probabilities (useCEMtr=true)  

#ifndef G4VPreCompoundTransitions_hh
#define G4VPreCompoundTransitions_hh 1

#include "G4Fragment.hh"

class G4VPreCompoundTransitions
{
public:

  G4VPreCompoundTransitions():useNGB(false),useCEMtr(false) {}
  virtual ~G4VPreCompoundTransitions() {}

  virtual G4double CalculateProbability(const G4Fragment& aFragment) = 0;
  virtual G4Fragment PerformTransition(const G4Fragment&  aFragment) = 0;
//J. M. Quesada (May.08) New virtual classes
  virtual G4double GetTransitionProb1()=0;
  virtual G4double GetTransitionProb2()=0;
  virtual G4double GetTransitionProb3()=0;

  // for never go back hypothesis (if useNGB=true, default=false)
  inline void UseNGB(G4bool use){useNGB=use;}
  //for use of CEM transition probabilities (if useCEMtr=true, defaut false)
  inline void UseCEMtr(G4bool use){useCEMtr=use;}

protected:
  G4bool useNGB;
  G4bool useCEMtr;
};

#endif
