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
//      For information related to this code contact:
//
//      File name:     G4CrossSectionSourcePtr
//
//      Author:        
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4VCrossSectionSource.hh"
#include "G4CrossSectionSourcePtr.hh"

G4CrossSectionSourcePtr::G4CrossSectionSourcePtr(G4VCrossSectionSource* x): x_(x)
{ }

G4CrossSectionSourcePtr& G4CrossSectionSourcePtr::operator= (const G4CrossSectionSourcePtr& xw)
{
  x_ = xw.x_; 
  return *this; 
}

G4bool G4CrossSectionSourcePtr::operator==(const G4CrossSectionSourcePtr& right) const
{
  return *(this->operator()()) == *right(); 
}

const G4VCrossSectionSource* G4CrossSectionSourcePtr::operator() () const 
{ return x_; }

G4VCrossSectionSource* G4CrossSectionSourcePtr::operator() ()
{ return x_; }








