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
#ifndef Nucleus_h
#define Nucleus_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4Nucleus.hh
//
// Version:             0.b.3
// Date:                29/02/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "g4std/iostream"
////////////////////////////////////////////////////////////////////////////////
//
class Nucleus
{
  // class description
  // The G4Nucleus class is used to contain information identifying an
  // isotope (a,z,e)
  //
  // class description - end
public: // with description
  Nucleus ();
  //    Default constructor
  //
  Nucleus (G4int a, G4int z, G4double e);
  //    Constructor defining new isotope with A,Z.E
  //
  ~Nucleus();
  //  Destructor
  
private:
  G4int a;
  G4int z;
  G4double e;

  //
  //
  // INLINE DECLARATIONS/DEFINITIONS:
  //
public: // with description
  inline  G4int GetA () const {return a;}
  //    Returns the value of a
  inline  G4int GetZ () const {return z;}
  //    Returns the value of z
  inline  G4double GetE () const {return e;}
  //    Returns the value of e

  //
  //
  // DECLARATIONS OF FRIENDS OF THE CLASS.
  //
  friend G4std::ostream &operator << (G4std::ostream &s, const Nucleus &q);

};
////////////////////////////////////////////////////////////////////////////////
#endif



