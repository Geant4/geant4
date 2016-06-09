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
// $Id: G4LowEIonFragmentation.hh,v 1.3 2006/06/29 20:58:04 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// by H.P. Wellisch

#ifndef G4LowEIonFragmentation_h
#define G4LowEIonFragmentation_h 1


#include "G4PreCompoundModel.hh"
#include "G4HadronicInteraction.hh"
#include "G4ExcitationHandler.hh"
#include "G4Fancy3DNucleus.hh"

class G4LowEIonFragmentation : public G4HadronicInteraction
{
public:
  
  G4LowEIonFragmentation(G4ExcitationHandler * const value) 
  {
    theHandler = value;
    theModel = new G4PreCompoundModel(theHandler);
  }

  G4LowEIonFragmentation() 
  {
    theHandler = new G4ExcitationHandler;
    theModel = new G4PreCompoundModel(theHandler);
  }

  ~G4LowEIonFragmentation() {delete theModel;}

private:
  
  G4LowEIonFragmentation(const G4LowEIonFragmentation &) : G4HadronicInteraction() {};
  
  const G4LowEIonFragmentation& operator=(const G4LowEIonFragmentation &right);

  G4bool operator==(const G4LowEIonFragmentation &right) const;
  
  G4bool operator!=(const G4LowEIonFragmentation &right) const;

public:
  G4HadFinalState * ApplyYourself(const G4HadProjectile & thePrimary, G4Nucleus & theNucleus);
  static G4double GetCrossSection() 
  {
//    clog << "area/millibarn = "<<area/millibarn<<G4endl;
//    clog << "hits = "<<hits<<G4endl;
//    clog << "totalTries = "<<totalTries<<G4endl;
    return area*static_cast<G4double>(hits)/static_cast<G4double>(totalTries)/millibarn;
  }
private:  

  
  G4HadFinalState theResult;
   
  G4PreCompoundModel * theModel;
  G4ExcitationHandler * theHandler;
  
private:
  static G4int hits;
  static G4int totalTries;
  static G4double area;

};


#endif


