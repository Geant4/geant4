//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ExcitationHandler.hh,v 1.11 2002/12/12 19:17:05 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (May 1998)
// Modif (30 June 1998) by V. Lara:
//      -Using G4ParticleTable and therefore G4IonTable
//       it can return all kind of fragments produced in 
//       deexcitation
//      -It uses default algorithms for:
//              Evaporation: G4StatEvaporation
//              MultiFragmentation: G4DummyMF (a dummy one)
//              Fermi Breakup model: G4StatFermiBreakUp


#ifndef G4ExcitationHandler_h
#define G4ExcitationHandler_h 1

#include "G4VMultiFragmentation.hh"
#include "G4VFermiBreakUp.hh"
#include "G4VEvaporation.hh"
#include "G4VPhotonEvaporation.hh"
#include "G4Fragment.hh"
#include "G4DynamicParticle.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
// needed for default models
#include "G4Evaporation.hh"
#include "G4StatMF.hh"
#include "G4FermiBreakUp.hh"
#include "G4PhotonEvaporation.hh"
#include "G4IonConstructor.hh" 

//#define debug

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

  void SetPhotonEvaporation(G4VPhotonEvaporation * const value);

  void SetMaxZForFermiBreakUp(G4int aZ);
  void SetMaxAForFermiBreakUp(G4int anA);
  void SetMaxAandZForFermiBreakUp(G4int anA,G4int aZ);
  void SetMinEForMultiFrag(G4double anE);

private:

  G4ReactionProductVector * Transform(G4FragmentVector * theFragmentVector) const;

  const G4VEvaporation * GetEvaporation() const;

  const G4VMultiFragmentation * GetMultiFragmentation() const;

  const G4VFermiBreakUp * GetFermiModel() const;

  const G4VPhotonEvaporation * GetPhotonEvaporation() const;

  const G4int GetMaxZ() const;
  const G4int GetMaxA() const;
  const G4double GetMinE() const;

  
#ifdef debug
  void CheckConservation(const G4Fragment & aFragment,
			 G4FragmentVector * Result) const;
#endif
private:
  
  G4VEvaporation *theEvaporation;
  
  G4VMultiFragmentation *theMultiFragmentation;
  
  G4VFermiBreakUp *theFermiModel;
 
  G4VPhotonEvaporation * thePhotonEvaporation;

  G4int maxZForFermiBreakUp;
  G4int maxAForFermiBreakUp;
  G4double minEForMultiFrag;

  G4ParticleTable *theTableOfParticles;

  G4bool MyOwnEvaporationClass;
  G4bool MyOwnMultiFragmentationClass;  
  G4bool MyOwnFermiBreakUpClass;
  G4bool MyOwnPhotonEvaporationClass;

    struct DeleteFragment 
    {
	template<typename T>
	void operator()(const T* ptr) const
	    {
		delete ptr;
	    }
    };


};



inline const G4VEvaporation * G4ExcitationHandler::GetEvaporation() const
{
  return theEvaporation;
}

inline void G4ExcitationHandler::SetEvaporation(G4VEvaporation *const  value)
{
  if (theEvaporation != 0 && MyOwnEvaporationClass) delete theEvaporation;
  MyOwnEvaporationClass = false;
  theEvaporation = value;
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


inline const G4VPhotonEvaporation * G4ExcitationHandler::GetPhotonEvaporation() const
{
  return thePhotonEvaporation;
}

inline void G4ExcitationHandler::SetPhotonEvaporation(G4VPhotonEvaporation *const  value)
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

inline const G4int G4ExcitationHandler::GetMaxZ() const
{
  return maxZForFermiBreakUp;
}

inline const G4int G4ExcitationHandler::GetMaxA() const 
{
  return maxAForFermiBreakUp;
}

inline const G4double G4ExcitationHandler::GetMinE() const 
{
  return minEForMultiFrag;
}


#endif
