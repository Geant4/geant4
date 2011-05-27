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
#ifndef DoubleBandFilter_h
#define DoubleBandFilter_h

#include "Analysis/src/VFilter.h"
#include "globals.hh"

class DoubleBandFilter : public TVANAFilter<G4double>
{
  public:
  
  DoubleBandFilter(G4double higherThan, G4double lowerEquals, G4String aName)
   : TVANAFilter<G4double>(aName)
  {
    theLow = higherThan;
    theHigh = lowerEquals;
  }
  
  G4bool Accept(G4double & anInput)
  {
    if(anInput<=theHigh && anInput>theLow) return true;
    return false;
  }
  
  G4double RelativeGeometricalAcceptance() 
  { 
    // Filters in cos(th)
    G4double result = (theHigh - theLow)/2.;
    G4cout << "acceptance "<< result << G4endl;
    return result;
  }
  
  private:
  
  G4double theLow;
  G4double theHigh;
};

#endif
