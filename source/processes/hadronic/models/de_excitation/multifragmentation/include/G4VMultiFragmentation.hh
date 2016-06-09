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
// $Id: G4VMultiFragmentation.hh,v 1.2 2005/06/04 13:27:48 jwellisc Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara 

#ifndef G4VMultiFragmentation_h
#define G4VMultiFragmentation_h 1

#include "globals.hh"
#include "G4Fragment.hh"

class G4VMultiFragmentation 
{
public:
  G4VMultiFragmentation();
  virtual ~G4VMultiFragmentation();
  
private:
  G4VMultiFragmentation(const G4VMultiFragmentation &right);
  
  const G4VMultiFragmentation & operator=(const G4VMultiFragmentation &right);
  G4bool operator==(const G4VMultiFragmentation &right) const;
  G4bool operator!=(const G4VMultiFragmentation &right) const;
  
public:
  virtual G4FragmentVector * BreakItUp(const G4Fragment &theNucleus) = 0;
};


#endif


