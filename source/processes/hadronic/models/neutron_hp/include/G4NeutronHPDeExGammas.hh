// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPDeExGammas.hh,v 1.4 1999-12-15 14:53:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPDeExGammas_h
#define G4NeutronHPDeExGammas_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "G4ReactionProductVector.hh"
#include "G4Gamma.hh"
#include "G4NeutronHPLevel.hh"
#include "G4NeutronHPGamma.hh"
#include "G4ReactionProduct.hh"

class G4NeutronHPDeExGammas
{
  public:
  
  G4NeutronHPDeExGammas() 
  {
    levelStart = NULL;
    levelSize = NULL;
    nLevels = 0;
    theLevels = NULL;
  }
  ~G4NeutronHPDeExGammas() 
  {
    if(levelStart!=NULL) delete [] levelStart;
    if(levelSize!=NULL) delete [] levelSize;
    if(theLevels!=NULL) delete [] theLevels;
  }
  
  void Init(G4std::ifstream & aDataFile);

  inline G4ReactionProductVector * GetDecayGammas(G4int aLevel)
  {
    if(aLevel>nLevels-1 || aLevel<0) return NULL;
    if(nLevels==0) return new G4ReactionProductVector();
    G4ReactionProductVector * result = new G4ReactionProductVector;
    G4DynamicParticleVector * theResult;

    theResult = theLevels[aLevel]. GetDecayGammas();
    G4ReactionProduct * theCurrent;
    G4int i;
    for(i=0; i<theResult->length(); i++)
    {
      theCurrent = new G4ReactionProduct;
      *theCurrent = *(theResult->at(i));
      delete theResult->at(i);
      G4double costheta = 2.*G4UniformRand()-1;
      G4double theta = acos(costheta);
      G4double phi = twopi*G4UniformRand();
      G4double sinth = sin(theta);
      G4double en = theCurrent->GetTotalMomentum();
      G4ThreeVector temp(en*sinth*cos(phi), en*sinth*sin(phi), en*costheta );
      theCurrent->SetMomentum( temp ) ;
      result->insert(theCurrent);
    }
    delete theResult;
    return result;
  }
  
  inline G4NeutronHPLevel * GetLevel(G4int i)
  {
    if(i>nLevels-1) return NULL;
    return theLevels+i;
  }
  
  inline G4int GetNumberOfLevels() { return nLevels; }
  
  inline G4double GetLevelEnergy(G4int aLevel)
  {
    if(aLevel>nLevels-1 || aLevel<0) return 0;
    G4double result = theLevels[aLevel].GetLevelEnergy();
    return result;
  }
  private:
  
  G4int * levelStart;
  G4int * levelSize;
  G4int nLevels;
  G4NeutronHPLevel * theLevels;
};

#endif
