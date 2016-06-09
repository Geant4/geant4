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
  if (this != &xw) {
	x_ = xw.x_;
  }
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








