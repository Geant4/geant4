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
// $Id: G4NeutronHPDeExGammas.hh,v 1.8 2002-12-12 19:18:11 gunter Exp $
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
    unsigned int i;
    for(i=0; i<theResult->size(); i++)
    {
      theCurrent = new G4ReactionProduct;
      *theCurrent = *(theResult->operator[](i));
      delete theResult->operator[](i);
      G4double costheta = 2.*G4UniformRand()-1;
      G4double theta = acos(costheta);
      G4double phi = twopi*G4UniformRand();
      G4double sinth = sin(theta);
      G4double en = theCurrent->GetTotalMomentum();
      G4ThreeVector temp(en*sinth*cos(phi), en*sinth*sin(phi), en*costheta );
      theCurrent->SetMomentum( temp ) ;
      result->push_back(theCurrent);
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
