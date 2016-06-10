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
// $Id: G4VEvaporation.hh 85443 2014-10-29 14:35:57Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations interface
//
// by V. Lara (Oct 1998) written from G4Evaporation.hh (May 1998)
//
// Modifications:
// 03 September 2008 by J. M. Quesada for external choice of inverse 
//    cross section option
// 06 September 2008 (JMQ) Also external choices have been added for 
//    superimposed Coulomb barrier (if useSICBis set true, by default is false) 
// 23 January 2012 by V.Ivanchenko added pointer of G4VPhotonEvaporation to 
//    the constructor
// 

#ifndef G4VEvaporation_h
#define G4VEvaporation_h 1

#include "globals.hh"
#include "G4Fragment.hh"
#include "G4VEvaporationFactory.hh"
#include "G4VEvaporationChannel.hh"
#include <vector>

class G4VEvaporationChannel;

class G4VEvaporation 
{
public:
  G4VEvaporation();
  virtual ~G4VEvaporation(); 

private:  
  G4VEvaporation(const G4VEvaporation &right);
  const G4VEvaporation & operator=(const G4VEvaporation &right);
  G4bool operator==(const G4VEvaporation &right) const;
  G4bool operator!=(const G4VEvaporation &right) const;
  
public:

  // primary fragment is copied to the new instance, the copy is deleted 
  // or is added to the list of products 
  virtual G4FragmentVector * BreakItUp(const G4Fragment &theNucleus) = 0;

  // new interface - vector of products is added to the provided vector
  // primary fragment is deleted or is modified and added to the list
  // of products 
  virtual void BreakFragment(G4FragmentVector*, G4Fragment* theNucleus);

  // definition of options
  virtual void Initialise();

  virtual void SetPhotonEvaporation(G4VEvaporationChannel* ptr);

  inline G4VEvaporationChannel* GetPhotonEvaporation();
  inline G4VEvaporationChannel* GetFissionChannel();

  // for inverse cross section choice
  inline void SetOPTxs(G4int opt) { OPTxs = opt;} 
  // for superimposed Coulomb Barrier for inverse cross sections 	
  inline void UseSICB(G4bool use) { useSICB = use; }	

protected:

  G4VEvaporationChannel* thePhotonEvaporation;

  G4int OPTxs;
  G4bool useSICB;

  std::vector<G4VEvaporationChannel*> * theChannels;
  G4VEvaporationFactory * theChannelFactory;

};

inline G4VEvaporationChannel* G4VEvaporation::GetPhotonEvaporation()
{
  return thePhotonEvaporation;
}

inline G4VEvaporationChannel* G4VEvaporation::GetFissionChannel()
{
  for(std::vector<G4VEvaporationChannel*>::const_iterator iC = theChannels->begin(),
      eC = theChannels->end(); iC!=eC; ++iC) {
    if((*iC)->GetName() == "fission")
      return *iC;
  }
  return 0;
}

#endif
