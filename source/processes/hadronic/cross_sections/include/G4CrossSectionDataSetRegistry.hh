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
// $Id: G4CrossSectionDataSetRegistry.hh,v 1.3 2009-08-08 16:21:31 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4CrossSectionDataSetRegistry
//
// Author  V.Ivanchenko  24.01.2009
//
// Modifications:
//
 
//
// Class Description
// This is a singleton keeping pointers to all cross section data sets
// Class Description - End

#ifndef G4CrossSectionDataSetRegistry_h
#define G4CrossSectionDataSetRegistry_h 1

#include <vector>
#include "globals.hh"

class G4VCrossSectionDataSet;

class G4CrossSectionDataSetRegistry
{
public:

  static G4CrossSectionDataSetRegistry* Instance();
  // access 
  
  ~G4CrossSectionDataSetRegistry();
  
  void Register(G4VCrossSectionDataSet*);
  //register new cross section

  void DeRegister(G4VCrossSectionDataSet*);
  //deregister cross section

  void Clean();
  //clean the store
    
private:

  G4CrossSectionDataSetRegistry();

  static G4CrossSectionDataSetRegistry* theInstance;
  
  std::vector <G4VCrossSectionDataSet*> xSections;

};

#endif
