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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:     G4CrossSectionSourcePtr
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4CROSSSECTIONSOURCEPTR_HH
#define G4CROSSSECTIONSOURCEPTR_HH

#include "globals.hh"

class G4VCrossSectionSource;

class G4CrossSectionSourcePtr
{
public:

  // This is a wrapper class to store pointers to G4VCrossSectionSource in vectors
  // Constructor 
  G4CrossSectionSourcePtr(G4VCrossSectionSource* x = 0);

  //Destructor
  ~G4CrossSectionSourcePtr() { }

  // Copy constructor
  G4CrossSectionSourcePtr(const G4CrossSectionSourcePtr& xw) : x_(xw.x_) { }

  // Operators

  const G4VCrossSectionSource* operator() () const;
  G4VCrossSectionSource* operator() ();

  G4CrossSectionSourcePtr& operator= (const G4CrossSectionSourcePtr& xw);

  G4bool operator== (const G4CrossSectionSourcePtr& right) const;

  G4bool operator< (const G4CrossSectionSourcePtr& right) { return false; }  

private:  

  G4VCrossSectionSource* x_;

};

#endif


















