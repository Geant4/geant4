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
// $Id: G4LowEIonFragmentation.hh 68028 2013-03-13 13:48:15Z gcosmo $
//
//---------------------------------------------------------------------------
//
// $Id: G4LowEIonFragmentation.hh 68028 2013-03-13 13:48:15Z gcosmo $
//
// ClassName:   G4LowEIonFragmentation
//
// Author:  H.P. Wellisch
//
// Modified:
// 01.06.2010 V.Ivanchenko moved constructor and destructor to the source
// 14.02.2013 V.Ivanchenko remove all static declarations
// 

#ifndef G4LowEIonFragmentation_h
#define G4LowEIonFragmentation_h 1


#include "G4PreCompoundModel.hh"
#include "G4HadronicInteraction.hh"
#include "G4ExcitationHandler.hh"

class G4ParticleDefinition;

class G4LowEIonFragmentation : public G4HadronicInteraction
{
public:
  
  G4LowEIonFragmentation(G4ExcitationHandler * const value);

  G4LowEIonFragmentation();

  virtual ~G4LowEIonFragmentation();

  virtual G4HadFinalState * ApplyYourself(const G4HadProjectile & thePrimary, 
					  G4Nucleus & theNucleus);

  inline G4double GetCrossSection() 
  {
    //    G4cout << "area/millibarn = "<<area/millibarn<<G4endl;
    //    G4cout << "hits = "<<hits<<G4endl;
    //    G4cout << "totalTries = "<<totalTries<<G4endl;
    return area*hits/(static_cast<G4double>(totalTries)*CLHEP::millibarn);
  }

private:  

  G4LowEIonFragmentation(const G4LowEIonFragmentation &);  
  const G4LowEIonFragmentation& operator=(const G4LowEIonFragmentation &right);
  G4bool operator==(const G4LowEIonFragmentation &right) const;
  G4bool operator!=(const G4LowEIonFragmentation &right) const;

  // Members
  
  G4HadFinalState theResult;
  const G4ParticleDefinition* proton;
   
  G4PreCompoundModel * theModel;
  G4ExcitationHandler * theHandler;
  
private:

  G4int hits;
  G4int totalTries;
  G4double area;

};


#endif


