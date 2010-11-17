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
// $Id: G4ExcitationHandler.hh,v 1.13 2010-11-17 16:20:31 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (May 1998)
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option
//
// Modif (30 June 1998) by V. Lara:
//      -Using G4ParticleTable and therefore G4IonTable
//       it can return all kind of fragments produced in 
//       deexcitation
//      -It uses default algorithms for:
//              Evaporation: G4StatEvaporation
//              MultiFragmentation: G4DummyMF (a dummy one)
//              Fermi Breakup model: G4StatFermiBreakUp
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option
// JMQ (06 September 2008) Also external choices have been added for 
// superimposed Coulomb barrier (if useSICBis set true, by default is false) 

#ifndef G4ExcitationHandler_h
#define G4ExcitationHandler_h 1

#include "G4VMultiFragmentation.hh"
#include "G4VFermiBreakUp.hh"
#include "G4VEvaporation.hh"
#include "G4VPhotonEvaporation.hh"
#include "G4VEvaporationChannel.hh"
#include "G4Fragment.hh"
#include "G4DynamicParticle.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
// needed for default models
#include "G4Evaporation.hh"
#include "G4StatMF.hh"
#include "G4FermiBreakUp.hh"
#include "G4PhotonEvaporation.hh"

class G4IonTable;

class G4ExcitationHandler 
{
public:
  G4ExcitationHandler(); 
  ~G4ExcitationHandler();

private:

  G4ExcitationHandler(const G4ExcitationHandler &right);
  const G4ExcitationHandler & operator=(const G4ExcitationHandler &right);
  G4bool operator==(const G4ExcitationHandler &right) const;
  G4bool operator!=(const G4ExcitationHandler &right) const;
  
public:

  G4ReactionProductVector * BreakItUp(const G4Fragment &theInitialState) const;

  void SetEvaporation(G4VEvaporation *const  value);

  void SetMultiFragmentation(G4VMultiFragmentation *const  value);

  void SetFermiModel(G4VFermiBreakUp *const  value);

  void SetPhotonEvaporation(G4VEvaporationChannel * const value);

  void SetMaxZForFermiBreakUp(G4int aZ);
  void SetMaxAForFermiBreakUp(G4int anA);
  void SetMaxAandZForFermiBreakUp(G4int anA,G4int aZ);
  void SetMinEForMultiFrag(G4double anE);

  // for inverse cross section choice
  inline void SetOPTxs(G4int opt);
  // for superimposed Coulomb Barrir for inverse cross sections
  inline void UseSICB();

private:

  void SetParameters();

  G4ReactionProductVector * Transform(G4FragmentVector * theFragmentVector) const;

  const G4VEvaporation * GetEvaporation() const;

  const G4VMultiFragmentation * GetMultiFragmentation() const;

  const G4VFermiBreakUp * GetFermiModel() const;

  const G4VEvaporationChannel * GetPhotonEvaporation() const;

  G4int GetMaxZ() const;
  G4int GetMaxA() const;
  G4double GetMinE() const;
  
#ifdef debug
  void CheckConservation(const G4Fragment & aFragment,
			 G4FragmentVector * Result) const;
#endif

private:
  
  G4VEvaporation *theEvaporation;
  
  G4VMultiFragmentation *theMultiFragmentation;
  
  G4VFermiBreakUp *theFermiModel;
 
  G4VEvaporationChannel * thePhotonEvaporation;

  G4int maxZForFermiBreakUp;
  G4int maxAForFermiBreakUp;
  G4double minEForMultiFrag;
  G4double minExcitation;

  G4IonTable* theTableOfIons;

  G4bool MyOwnEvaporationClass;
  G4bool MyOwnMultiFragmentationClass;  
  G4bool MyOwnFermiBreakUpClass;
  G4bool MyOwnPhotonEvaporationClass;

  G4int OPTxs;
  G4bool useSICB;
  
  struct DeleteFragment 
  {
    template<typename T>
    void operator()(const T* ptr) const
    {
      delete ptr;
    }
  };
  
};

inline void G4ExcitationHandler::SetOPTxs(G4int opt) 
{ 
  OPTxs = opt; 
  SetParameters();
}

inline void G4ExcitationHandler::UseSICB()
{ 
  useSICB = true; 
  SetParameters();
}

inline const G4VEvaporation * G4ExcitationHandler::GetEvaporation() const
{
  return theEvaporation;
}

inline void G4ExcitationHandler::SetEvaporation(G4VEvaporation *const  value)
{
  if (theEvaporation != 0 && MyOwnEvaporationClass) delete theEvaporation;
  MyOwnEvaporationClass = false;
  theEvaporation = value;
  SetParameters();
}

inline const G4VMultiFragmentation * G4ExcitationHandler::GetMultiFragmentation() const
{
  return theMultiFragmentation;
}

inline void G4ExcitationHandler::SetMultiFragmentation(G4VMultiFragmentation *const  value)
{
  if (theMultiFragmentation != 0 && MyOwnMultiFragmentationClass) delete theMultiFragmentation;
  MyOwnMultiFragmentationClass = false;
  theMultiFragmentation = value;
}

inline const G4VFermiBreakUp * G4ExcitationHandler::GetFermiModel() const
{
  return theFermiModel;
}

inline void G4ExcitationHandler::SetFermiModel(G4VFermiBreakUp *const  value)
{
  if (theFermiModel != 0 && MyOwnFermiBreakUpClass) delete theFermiModel;
  MyOwnFermiBreakUpClass = false;
  theFermiModel = value;
}


inline const G4VEvaporationChannel * G4ExcitationHandler::GetPhotonEvaporation() const
{
  return thePhotonEvaporation;
}

inline void G4ExcitationHandler::SetPhotonEvaporation(G4VEvaporationChannel *const  value)
{
  if (thePhotonEvaporation != 0 && MyOwnPhotonEvaporationClass) delete thePhotonEvaporation;
  MyOwnPhotonEvaporationClass = false;
  thePhotonEvaporation = value;
}

inline void G4ExcitationHandler::SetMaxZForFermiBreakUp(G4int aZ)
{
  maxZForFermiBreakUp = aZ;
}

inline void G4ExcitationHandler::SetMaxAForFermiBreakUp(G4int anA)
{
  maxAForFermiBreakUp = anA;
}

inline void G4ExcitationHandler::SetMaxAandZForFermiBreakUp(G4int anA, G4int aZ)
{
  maxAForFermiBreakUp = anA;
  maxZForFermiBreakUp = aZ;
}

inline void G4ExcitationHandler::SetMinEForMultiFrag(G4double anE)
{
  minEForMultiFrag = anE;
}

inline G4int G4ExcitationHandler::GetMaxZ() const
{
  return maxZForFermiBreakUp;
}

inline G4int G4ExcitationHandler::GetMaxA() const 
{
  return maxAForFermiBreakUp;
}

inline G4double G4ExcitationHandler::GetMinE() const 
{
  return minEForMultiFrag;
}


#endif
