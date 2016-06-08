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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VPreCompoundModel.hh,v 1.6 2001/08/01 17:08:15 hpw Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//

#ifndef G4VPreCompoundModel_h
#define G4VPreCompoundModel_h 1

// Class Description
// Base class for pre-equilibrium decay models in geant4. By merit of inheriting
// from this class a pre-equilibrium decay model can be used in conjunction with
// any cascade, string parton model or other high energy generator in the
// generation of final states for inelastic scattering.
// Class Description - End

#include "G4HadronicInteraction.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
class G4Track;
class G4VParticleChange;
class G4Nucleus;
class G4Fragment;
class G4ExcitationHandler;

class G4VPreCompoundModel : public G4HadronicInteraction
{
public:

  G4VPreCompoundModel(G4ExcitationHandler *const value);
  virtual ~G4VPreCompoundModel() {};
  
protected:
  // default constructor
  G4VPreCompoundModel() {};
private:
  // copy constructor
  G4VPreCompoundModel(const G4VPreCompoundModel &right) {};
  // operators
  const G4VPreCompoundModel& operator=(const G4VPreCompoundModel &right);
  G4bool operator==(const G4VPreCompoundModel &right) const;
  G4bool operator!=(const G4VPreCompoundModel &right) const;

public:
  virtual G4VParticleChange * 
          ApplyYourself(const G4Track & thePrimary, G4Nucleus & theNucleus) = 0;
  
  virtual G4ReactionProductVector* 
          DeExcite(const G4Fragment& aFragment) const = 0;

  void SetExcitationHandler(G4ExcitationHandler *const  value);
    
protected:

  const G4ExcitationHandler * GetExcitationHandler() const;
  
private:
  G4ExcitationHandler *theExcitationHandler;
};



inline const G4ExcitationHandler* G4VPreCompoundModel::GetExcitationHandler() const
{
  return theExcitationHandler;
}

inline void G4VPreCompoundModel::SetExcitationHandler(G4ExcitationHandler *const  value)
{
  theExcitationHandler = value;
}




#endif
