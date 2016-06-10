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
// $Id: G4VPreCompoundTransitions.hh 68028 2013-03-13 13:48:15Z gcosmo $
//
// J. M. Quesada (May 08). New virtual classes have been added Prob1,2,3 
// JMQ (06 September 2008) Also external choices have been added for:
//                      - "never go back"  hipothesis (useNGB=true) 
//                      - CEM transition probabilities (useCEMtr=true)  
// 20.08.2010 V.Ivanchenko move constructor and destructor to the source

#ifndef G4VPreCompoundTransitions_hh
#define G4VPreCompoundTransitions_hh 1

#include "G4Fragment.hh"

class G4VPreCompoundTransitions
{
public:

  G4VPreCompoundTransitions();
  virtual ~G4VPreCompoundTransitions();

  virtual G4double CalculateProbability(const G4Fragment& aFragment) = 0;
  virtual void PerformTransition(G4Fragment&  aFragment) = 0;

  inline G4double GetTransitionProb1() const;

  inline G4double GetTransitionProb2() const;

  inline G4double GetTransitionProb3() const;

  // for never go back hypothesis (if useNGB=true, default=false)
  inline void UseNGB(G4bool use){useNGB=use;}
  //for use of CEM transition probabilities (if useCEMtr=true, defaut false)
  inline void UseCEMtr(G4bool use){useCEMtr=use;}

private:

  G4VPreCompoundTransitions(const G4VPreCompoundTransitions &);
  const G4VPreCompoundTransitions& operator=(const G4VPreCompoundTransitions &right);
  G4bool operator==(const G4VPreCompoundTransitions &right) const;
  G4bool operator!=(const G4VPreCompoundTransitions &right) const;

protected:

  G4bool useNGB;
  G4bool useCEMtr;

  G4double TransitionProb1;
  G4double TransitionProb2;
  G4double TransitionProb3;
};

  //J. M.Quesada (May. 08)
  inline G4double G4VPreCompoundTransitions::GetTransitionProb1() const
  {
    return TransitionProb1;
  }
  inline G4double G4VPreCompoundTransitions::GetTransitionProb2() const
  {
    return TransitionProb2;
  }
  inline G4double G4VPreCompoundTransitions::GetTransitionProb3() const
  {
    return TransitionProb3;
  }


#endif
