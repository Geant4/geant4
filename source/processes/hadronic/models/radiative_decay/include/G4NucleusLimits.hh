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
#ifndef G4NucleusLimits_h
#define G4NucleusLimits_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4NucleusLimits.hh
//
// Version:             0.b.4
// Date:                14/04/00
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
// 13 April 2000, F Lei, DERA UK
// 0.b.4 release. No change to this file     
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "g4std/iostream"

////////////////////////////////////////////////////////////////////////////////
//
class G4NucleusLimits
{
  // class description
  //   The G4NucleusLimits class is used to contain information identifying a
  //   range bounded by A [aMin, aMax], and Z on [zMin, zMax].
  //
  // class description - end

public: // with description
  G4NucleusLimits ();
  //    Default constructor where: aMin=1, aMax=240, aMin=1, and
  //    zMax=100.
  //
  G4NucleusLimits
  (G4int aMin1, G4int aMax1, G4int zMin1, G4int zMax1);
  //    Constructor defining new values for aMin, aMax, zMin, and zMax 
  //    respectively.
  //
  ~G4NucleusLimits();
  //  Destructor
  
private:
  G4int aMin;
  G4int aMax;
  G4int zMin;
  G4int zMax;

  //
  //
  // INLINE DECLARATIONS/DEFINTIONS:
  //
public: // with description

  inline  G4int GetAMin () const {return aMin;}
  //    Returns the value of aMin.
  //
  inline  G4int GetAMax () const {return aMax;}
  //    Returns the value of aMax.
  //
  inline  G4int GetZMin () const {return zMin;}
  //    Returns the value of zMin.
  //
  inline  G4int GetZMax () const {return zMax;}
  //    Return the value of zmax.
  //


  friend G4std::ostream &operator << (G4std::ostream &s, const G4NucleusLimits &q);

};
////////////////////////////////////////////////////////////////////////////////
#endif

