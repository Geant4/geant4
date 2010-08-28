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
// $Id: G4VIntraNuclearTransportModel.hh,v 1.5 2010-08-28 15:53:50 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, A. Feliciello, 30th June 1998
//      A.Pavliouk 26.11.98
//          In Set...() methods a pointer is deleted now before new
//          value will be asigned.
// -----------------------------------------------------------------------------

#ifndef G4VIntraNuclearTransportModel_h
#define G4VIntraNuclearTransportModel_h 1

// Class Description
// Base class for intra-nuclear transport models in geant4. By merit of inheriting
// from this class a intra-nuclear transport model can be used in conjunction with
// any precompound, string parton model or other high energy generator in the
// generation of final states for inelastic scattering.
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

  G4VIntraNuclearTransportModel(const G4String& modelName = "CascadeModel");

  virtual ~G4VIntraNuclearTransportModel();

  virtual 
  G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries,
				     G4V3DNucleus* theNucleus) = 0;

  inline void SetDeExcitation(G4VPreCompoundModel* const  value);

  inline void Set3DNucleus(G4V3DNucleus* const value);

  inline const G4String& GetModelName() const;

private:

  G4VIntraNuclearTransportModel(const G4VIntraNuclearTransportModel& right);
  const G4VIntraNuclearTransportModel& operator=(const G4VIntraNuclearTransportModel &right);
  int operator==(const G4VIntraNuclearTransportModel& right) const;
  int operator!=(const G4VIntraNuclearTransportModel& right) const;

protected:

  inline G4V3DNucleus* Get3DNucleus() const;

  inline G4VPreCompoundModel* GetDeExcitation() const;

  G4String theTransportModelName;

  G4V3DNucleus* the3DNucleus;

  G4VPreCompoundModel* theDeExcitation;
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
G4VIntraNuclearTransportModel::SetDeExcitation(G4VPreCompoundModel* const  value)
{
   delete theDeExcitation; 
   theDeExcitation = value;
}

#endif


