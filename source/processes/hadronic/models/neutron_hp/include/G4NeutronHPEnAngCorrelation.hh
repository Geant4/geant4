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
// $Id: G4NeutronHPEnAngCorrelation.hh,v 1.5 2001-07-11 10:06:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPEnAngCorrelation_h
#define G4NeutronHPEnAngCorrelation_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "globals.hh"
#include "G4NeutronHPProduct.hh"
#include "G4ReactionProduct.hh"

class G4NeutronHPEnAngCorrelation
{
  public:
  G4NeutronHPEnAngCorrelation()
  {
    theProducts = NULL;
    inCharge = false;
    theTotalMeanEnergy = -1.;
  }
  ~G4NeutronHPEnAngCorrelation()
  {
    if(theProducts!=NULL) delete [] theProducts;
  }
  
  inline void Init(G4std::ifstream & aDataFile)
  {
    inCharge = true;
    aDataFile>>targetMass>>frameFlag>>nProducts;
    theProducts = new G4NeutronHPProduct[nProducts];
    for(G4int i=0; i<nProducts; i++)
    {
      theProducts[i].Init(aDataFile);
    }
  }
  
  G4ReactionProduct * SampleOne(G4double anEnergy);

  G4ReactionProductVector * Sample(G4double anEnergy);
  
  inline void SetTarget(G4ReactionProduct & aTarget)
  {
    theTarget = aTarget;
    for(G4int i=0;i<nProducts;i++)theProducts[i].SetTarget(&theTarget);
  }
  
  inline void SetNeutron(G4ReactionProduct & aNeutron)
  {
    theNeutron = aNeutron;
    for(G4int i=0;i<nProducts;i++)theProducts[i].SetNeutron(&theNeutron);
  }
  
  inline G4bool InCharge()
  {
    return inCharge;
  }
  
  inline G4double GetTargetMass() { return targetMass; }
  
  G4double GetTotalMeanEnergy()
  {
     // cashed in 'sample' call
    return theTotalMeanEnergy; 
  }
  
  private:
   
  // data members
  
  G4double targetMass;
  G4int frameFlag; // =1: Target rest frame; =2: CMS system; incident always in lab
  G4int nProducts;
  G4NeutronHPProduct * theProducts;
  G4bool inCharge;
    
  // Utility quantities
  
  G4ReactionProduct theTarget;
  G4ReactionProduct theNeutron;
  
  // cashed values
  
  G4double theTotalMeanEnergy;
  
};

#endif
