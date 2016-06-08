// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPreCompoundModel.hh,v 1.2.8.1 1999/12/07 20:51:44 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//

#ifndef G4VPreCompoundModel_h
#define G4VPreCompoundModel_h 1

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
