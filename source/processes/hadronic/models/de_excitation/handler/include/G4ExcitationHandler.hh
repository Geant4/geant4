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
// $Id: G4ExcitationHandler.hh 104078 2017-05-10 14:50:23Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (May 1998)
//
// Modifications:
// 30 June 1998 by V. Lara:
//      -Using G4ParticleTable and therefore G4IonTable
//       it can return all kind of fragments produced in 
//       deexcitation
//      -It uses default algorithms for:
//              Evaporation: G4StatEvaporation
//              MultiFragmentation: G4DummyMF (a dummy one)
//              Fermi Breakup model: G4StatFermiBreakUp
//
// 03 September 2008 by J. M. Quesada for external choice of inverse 
//    cross section option
// 06 September 2008 JMQ Also external choices have been added for 
//    superimposed Coulomb barrier (if useSICBis set true, by default is false)  
// 23 January 2012 by V.Ivanchenko remove obsolete data members; added access
//    methods to deexcitation components
//                   

#ifndef G4ExcitationHandler_h
#define G4ExcitationHandler_h 1

#include "globals.hh"
#include "G4Fragment.hh"
#include "G4ReactionProductVector.hh"
#include "G4IonTable.hh"
#include "G4DeexPrecoParameters.hh"

class G4VMultiFragmentation;
class G4VFermiBreakUp;
class G4VEvaporation;
class G4VEvaporationChannel;
class G4NistManager;

class G4ExcitationHandler 
{
public:

  explicit G4ExcitationHandler(); 
  ~G4ExcitationHandler();

  G4ReactionProductVector* BreakItUp(const G4Fragment &theInitialState);

  // short model description used for automatic web documentation
  void ModelDescription(std::ostream& outFile) const;

  void Initialise();

  // user defined sub-models
  // deletion is responsibility of this handler if isLocal=true 
  void SetEvaporation(G4VEvaporation* ptr, G4bool isLocal=false);
  void SetMultiFragmentation(G4VMultiFragmentation* ptr);
  void SetFermiModel(G4VFermiBreakUp* ptr);
  void SetPhotonEvaporation(G4VEvaporationChannel* ptr);
  void SetDeexChannelsType(G4DeexChannelType val);

  //======== Obsolete methods to be removed =====

  // parameters of sub-models
  inline void SetMaxZForFermiBreakUp(G4int aZ);
  inline void SetMaxAForFermiBreakUp(G4int anA);
  inline void SetMaxAandZForFermiBreakUp(G4int anA,G4int aZ);
  inline void SetMinEForMultiFrag(G4double anE);

  // access methods
  inline G4VEvaporation* GetEvaporation();
  inline G4VMultiFragmentation* GetMultiFragmentation();
  inline G4VFermiBreakUp* GetFermiModel();
  inline G4VEvaporationChannel* GetPhotonEvaporation();

  // for inverse cross section choice
  inline void SetOPTxs(G4int opt);
  // for superimposed Coulomb Barrir for inverse cross sections
  inline void UseSICB();

  //==============================================

private:

  void SetParameters();

  G4ExcitationHandler(const G4ExcitationHandler &right) = delete;
  const G4ExcitationHandler & operator
  =(const G4ExcitationHandler &right) = delete;
  G4bool operator==(const G4ExcitationHandler &right) const = delete;
  G4bool operator!=(const G4ExcitationHandler &right) const = delete;
  
  G4VEvaporation* theEvaporation;  
  G4VMultiFragmentation* theMultiFragmentation;
  G4VFermiBreakUp* theFermiModel;
  G4VEvaporationChannel* thePhotonEvaporation;

  const G4ParticleDefinition* electron;
  G4int icID;

  G4int maxZForFermiBreakUp;
  G4int maxAForFermiBreakUp;
  G4double minEForMultiFrag;
  G4double minExcitation;

  G4IonTable* theTableOfIons;
  G4NistManager* nist;

  G4int  fVerbose;
  G4bool isInitialised;
  G4bool isEvapLocal;

  // list of fragments to store final result   
  std::vector<G4Fragment*> theResults;

  // list of fragments to store intermediate result   
  std::vector<G4Fragment*> results;

  // list of fragments to apply PhotonEvaporation 
  std::vector<G4Fragment*> thePhotoEvapList;

  // list of fragments to apply Evaporation or Fermi Break-Up
  std::vector<G4Fragment*> theEvapList;          
};

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
  SetMaxAForFermiBreakUp(anA);
  SetMaxZForFermiBreakUp(aZ);
}

inline void G4ExcitationHandler::SetMinEForMultiFrag(G4double anE)
{
  minEForMultiFrag = anE;
}

inline G4VEvaporation* G4ExcitationHandler::GetEvaporation()
{
  return theEvaporation;
}

inline G4VMultiFragmentation* G4ExcitationHandler::GetMultiFragmentation()
{
  return theMultiFragmentation;
}

inline G4VFermiBreakUp* G4ExcitationHandler::GetFermiModel()
{
  return theFermiModel;
}

inline G4VEvaporationChannel* G4ExcitationHandler::GetPhotonEvaporation()
{
  return thePhotonEvaporation;
}

inline void G4ExcitationHandler::SetOPTxs(G4int) 
{}

inline void G4ExcitationHandler::UseSICB()
{}

#endif
