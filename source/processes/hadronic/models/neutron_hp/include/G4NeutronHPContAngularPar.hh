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
// $Id: G4NeutronHPContAngularPar.hh,v 1.13 2008-07-21 23:26:29 tkoi Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 080718 Add ClearHistories method and related class member
//
#ifndef G4NeutronHPContAngularPar_h
#define G4NeutronHPContAngularPar_h 1

#include "G4ios.hh"
#include <fstream>
#include "globals.hh"
#include "G4NeutronHPList.hh"
#include "G4ReactionProduct.hh"
#include "G4NeutronHPInterpolator.hh"
#include "G4InterpolationManager.hh"

class G4NeutronHPContAngularPar
{
  public:
  
  G4NeutronHPContAngularPar()
  {
    theAngular = 0;
    currentMeanEnergy = -2;
     fresh = true;
  }
  ~G4NeutronHPContAngularPar()
  {
    if(theAngular!=0) delete [] theAngular;
  }
  
  void Init(std::ifstream & aDataFile);
  
  G4ReactionProduct * Sample(G4double anEnergy, G4double massCode, G4double mass, 
                             G4int angularRep, G4int interpol);
  
  G4double GetEnergy() { return theEnergy; }
  
  void SetPrimary(G4ReactionProduct * aPrimary)
  {
    thePrimary = aPrimary;
  }
  
  void SetTarget(G4ReactionProduct * aTarget)
  {
    theTarget = aTarget;
  }
  
  void SetTargetCode(G4double aTargetCode) { theTargetCode = aTargetCode; }
  
  void SetInterpolation(G4int theInterpolation)
  {
    theManager.Init(theInterpolation, nEnergies); // one range only
  }

  void Merge(G4double anEnergy, G4InterpolationScheme & aScheme, 
             G4NeutronHPContAngularPar & store1, 
             G4NeutronHPContAngularPar & store2) // hmmmm, this interpolates legendre coefficients. Dangerous @@@
  {
    nDiscreteEnergies = store1.nDiscreteEnergies;
    nAngularParameters = store1.nAngularParameters;
    nEnergies = store1.nEnergies;
    theManager = store1.theManager;
    theEnergy = anEnergy;
    if(theAngular != 0) delete [] theAngular;
    theAngular = new G4NeutronHPList[nEnergies];
    G4int i, ii;
    G4double value;
    for(i=0; i<nEnergies; i++)
    {
      theAngular[i].SetLabel(store1.theAngular[i].GetLabel());
      for(ii=0; ii<nAngularParameters; ii++)
      {
//        G4cout <<"test "<<i<<" "<<store1.theEnergy<<" "<<store2.theEnergy<<" "
//             << store1.theAngular[i].GetValue(ii)<<" "<<
//             store2.theAngular[i].GetValue(ii)<<G4endl;
        value = theInt.Interpolate(aScheme, anEnergy, 
                                   store1.theEnergy, store2.theEnergy,
                                   store1.theAngular[i].GetValue(ii),
                                   store2.theAngular[i].GetValue(ii));
        theAngular[i].SetValue(ii, value);
      }
    }
  };
  
  G4double MeanEnergyOfThisInteraction()
  {
    G4double result;
    if(currentMeanEnergy<-1)
    {
      return 0;
      // throw G4HadronicException(__FILE__, __LINE__, "G4NeutronHPContAngularPar: Logical error in Product class");
    }
    else
    {
      result = currentMeanEnergy;
    }
    currentMeanEnergy = -2;
    return result;
  }
  
  private:
  
  // incoming particle
  G4double theEnergy; 
  
  // number of exit channel energies
  G4int nEnergies; 
  // number of discrete exit channels
  G4int nDiscreteEnergies;
  // number of angular paramerers per channel
  G4int nAngularParameters;
  // knows the interpolation between List labels
  G4InterpolationManager theManager; 
  // on per exit-channel energy
  G4NeutronHPList * theAngular; 
  
  private:
  
  G4NeutronHPInterpolator theInt;
  
  G4double theTargetCode;
  G4ReactionProduct * theTarget;
  G4ReactionProduct * thePrimary;
  
  G4double currentMeanEnergy;

//080718
   public:
      void ClearHistories(){ fresh = true; };
   private:
      G4bool fresh; 
      G4double remaining_energy; // represent energy rest of cascade chain
};
#endif
