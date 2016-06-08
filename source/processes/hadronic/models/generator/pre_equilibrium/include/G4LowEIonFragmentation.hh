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
// $Id: G4LowEIonFragmentation.hh,v 1.6 2001/10/06 08:17:01 hpw Exp $
// GEANT4 tag $Name: geant4-04-01 $
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
  
  G4LowEIonFragmentation(const G4LowEIonFragmentation &right) {};
  
  const G4LowEIonFragmentation& operator=(const G4LowEIonFragmentation &right);

  G4bool operator==(const G4LowEIonFragmentation &right) const;
  
  G4bool operator!=(const G4LowEIonFragmentation &right) const;

public:
  G4VParticleChange * ApplyYourself(const G4Track & thePrimary, G4Nucleus & theNucleus);
  static G4double GetCrossSection() 
  {
//    clog << "area/millibarn = "<<area/millibarn<<G4endl;
//    clog << "hits = "<<hits<<G4endl;
//    clog << "totalTries = "<<totalTries<<G4endl;
    return area*G4double(hits)/G4double(totalTries)/millibarn;
  }
private:  

  
  G4ParticleChange theResult;
   
  G4PreCompoundModel * theModel;
  G4ExcitationHandler * theHandler;
  
private:
  static G4int hits;
  static G4int totalTries;
  static G4double area;

};


#endif


