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
#ifndef VANAFilter_h
#define VANAFilter_h

template <class Type>
class TVANAFilter
{
  public:
  
  TVANAFilter(G4String aName) : theName(aName) {}
  virtual G4bool Accept(Type & anInput) = 0;
  virtual G4double RelativeGeometricalAcceptance() { return 1;}
  G4String GetName() {return theName;}
  
  private:
  
  TVANAFilter() {}
  
  private:
  
  G4String theName;
};

#endif
