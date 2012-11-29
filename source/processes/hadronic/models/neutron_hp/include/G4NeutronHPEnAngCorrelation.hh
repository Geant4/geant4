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
#ifndef G4NeutronHPEnAngCorrelation_h
#define G4NeutronHPEnAngCorrelation_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <fstream>
#include "globals.hh"
#include "G4NeutronHPProduct.hh"
#include "G4ReactionProduct.hh"

class G4NeutronHPEnAngCorrelation
{
  public:
  G4NeutronHPEnAngCorrelation()
  {
    theProducts = 0;
    inCharge = false;
    theTotalMeanEnergy = -1.;
  }
  ~G4NeutronHPEnAngCorrelation()
  {
    if(theProducts!=0) delete [] theProducts;
  }
  
  inline void Init(std::ifstream & aDataFile)
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
