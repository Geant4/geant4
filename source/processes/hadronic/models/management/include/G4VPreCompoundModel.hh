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
// $Id$
//

#ifndef G4VPreCompoundModel_h
#define G4VPreCompoundModel_h 1

// -----------------------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation 1998
//
//      V.Ivanchenko 03.01.2012
//          Added G4ExcitationHandler pointer to the constructor and cleanup
// -----------------------------------------------------------------------------

// Class Description
// Base class for pre-equilibrium decay models in geant4. By merit of 
// inheriting from this class a pre-equilibrium decay model can be used 
// in conjunction with any cascade, string parton model or other high 
// energy generator in the generation of final states for inelastic scattering.
// Class Description - End

#include "G4HadronicInteraction.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"

class G4HadProjectile;
class G4HadFinalState;
class G4Nucleus;
class G4Fragment;
class G4ExcitationHandler;

class G4VPreCompoundModel : public G4HadronicInteraction
{
public:

  G4VPreCompoundModel(G4ExcitationHandler* ptr = 0, 
                      const G4String& modelName = "PrecompoundModel");

  virtual ~G4VPreCompoundModel();
  
private:

  // default constructor
  //G4VPreCompoundModel();
  // copy constructor
  G4VPreCompoundModel(const G4VPreCompoundModel &);
  // operators
  const G4VPreCompoundModel& operator=(const G4VPreCompoundModel &right);
  G4bool operator==(const G4VPreCompoundModel &right) const;
  G4bool operator!=(const G4VPreCompoundModel &right) const;

public:

  virtual G4HadFinalState * 
          ApplyYourself(const G4HadProjectile & thePrimary, G4Nucleus & theNucleus) = 0;
  
  virtual G4ReactionProductVector* DeExcite(G4Fragment& aFragment) = 0;

  inline void SetExcitationHandler(G4ExcitationHandler* ptr);
    
  inline G4ExcitationHandler* GetExcitationHandler() const;
  
private:

  G4ExcitationHandler* theExcitationHandler;
};

inline void G4VPreCompoundModel::SetExcitationHandler(G4ExcitationHandler* ptr)
{
  theExcitationHandler = ptr;
}

inline G4ExcitationHandler* G4VPreCompoundModel::GetExcitationHandler() const
{
  return theExcitationHandler;
}

#endif
