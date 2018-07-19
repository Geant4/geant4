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
// $Id: G4VIntraNuclearTransportModel.hh 96490 2016-04-19 06:57:04Z gcosmo $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, A. Feliciello, 30th June 1998
//      A.Pavliouk 26.11.98
//          In Set...() methods a pointer is deleted now before new
//          value will be asigned.
//	M.Kelsey 07.03.2011
//	    Add data member and Set method to store original projectile
//      V.Ivanchenko 03.01.2012
//          Added G4VPreCompoundModel pointer to the constructor and cleanup
//      V. Uzhinsky Nov. 2012
//          Added method PropagateNuclNucl for simulation of nucleus-nucleus inter. 
// -----------------------------------------------------------------------------

#ifndef G4VIntraNuclearTransportModel_h
#define G4VIntraNuclearTransportModel_h 1

// Class Description
// Base class for intra-nuclear transport models in geant4. By merit 
// of inheriting from this class a intra-nuclear transport model can 
// be used in conjunction with any precompound, string parton model 
// or other high energy generator in the generation of final states 
// for inelastic scattering.
// Class Description - End

#include "G4V3DNucleus.hh"
#include "G4VPreCompoundModel.hh"
#include "G4HadronicInteraction.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4HadProjectile.hh"
#include "G4HadFinalState.hh"

class G4KineticTrackVector;

class G4VIntraNuclearTransportModel : public G4HadronicInteraction
{
public:

  explicit G4VIntraNuclearTransportModel(const G4String& mName = "CascadeModel",
					 G4VPreCompoundModel* ptr = nullptr);

  virtual ~G4VIntraNuclearTransportModel();

  virtual 
  G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries,
				     G4V3DNucleus* theNucleus) = 0;

  virtual 
  G4ReactionProductVector* PropagateNuclNucl(G4KineticTrackVector* theSecondaries,
				     G4V3DNucleus* theNucleus,
                                     G4V3DNucleus* theProjectileNucleus);    // Uzhi Nov. 2012

  inline void SetDeExcitation(G4VPreCompoundModel* ptr);

  inline void Set3DNucleus(G4V3DNucleus* const value);

  inline void SetPrimaryProjectile(const G4HadProjectile &aPrimary);

  inline const G4String& GetModelName() const;

  virtual void ModelDescription(std::ostream& outFile) const ;
  virtual void PropagateModelDescription(std::ostream& outFile) const ;

private:

  G4VIntraNuclearTransportModel(const G4VIntraNuclearTransportModel& right) = delete;
  const G4VIntraNuclearTransportModel& operator=(const G4VIntraNuclearTransportModel &right) = delete;
  int operator==(const G4VIntraNuclearTransportModel& right) const = delete;
  int operator!=(const G4VIntraNuclearTransportModel& right) const = delete;

protected:

  inline G4V3DNucleus* Get3DNucleus() const;

  inline G4VPreCompoundModel* GetDeExcitation() const;

  inline const G4HadProjectile* GetPrimaryProjectile() const;

  G4String theTransportModelName;

  G4V3DNucleus* the3DNucleus;

  G4VPreCompoundModel* theDeExcitation;

  const G4HadProjectile* thePrimaryProjectile;
};

inline const G4String& G4VIntraNuclearTransportModel::GetModelName() const
{
  return theTransportModelName;
}

inline G4V3DNucleus* G4VIntraNuclearTransportModel::Get3DNucleus() const
{
  return the3DNucleus;
}

inline void G4VIntraNuclearTransportModel::Set3DNucleus(G4V3DNucleus* const value)
{
  delete the3DNucleus;  the3DNucleus = value;
}

inline G4VPreCompoundModel* G4VIntraNuclearTransportModel::GetDeExcitation() const
{
  return theDeExcitation;
}

inline void 
G4VIntraNuclearTransportModel::SetDeExcitation(G4VPreCompoundModel* value)
{
  // previous pre-compound model will be deleted at the end of job 
  theDeExcitation = value;
}

inline const G4HadProjectile* 
G4VIntraNuclearTransportModel::GetPrimaryProjectile() const
{
  return thePrimaryProjectile;
}

inline void  
G4VIntraNuclearTransportModel::SetPrimaryProjectile(const G4HadProjectile &aPrimary)
{
  // NOTE:  Previous pointer is NOT deleted: passed by reference, no ownership
  thePrimaryProjectile = &aPrimary;
}

#endif


