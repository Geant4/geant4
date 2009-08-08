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
// $Id: G4CrossSectionDataSetRegistry.cc,v 1.4 2009-08-08 16:21:31 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:    G4CrossSectionDataSetRegistry
//
// Author  V.Ivanchenko  24.01.2009
//
// Modifications:
//

#include "G4CrossSectionDataSetRegistry.hh"
#include "G4VCrossSectionDataSet.hh"

G4CrossSectionDataSetRegistry* G4CrossSectionDataSetRegistry::theInstance = 0;

G4CrossSectionDataSetRegistry* G4CrossSectionDataSetRegistry::Instance()
{
  if(0 == theInstance) {
    static G4CrossSectionDataSetRegistry manager;
    theInstance = &manager;
  }
  return theInstance;
}

G4CrossSectionDataSetRegistry::G4CrossSectionDataSetRegistry()
{}

G4CrossSectionDataSetRegistry::~G4CrossSectionDataSetRegistry()
{
  Clean();
}

void G4CrossSectionDataSetRegistry::Clean()
{
  std::vector<G4VCrossSectionDataSet*>::iterator it = xSections.begin();
  for (; it != xSections.end(); ++it) {delete *it;}
  xSections.clear();
}

void G4CrossSectionDataSetRegistry::Register(G4VCrossSectionDataSet* p)
{
  if(!p) return;
  std::vector<G4VCrossSectionDataSet*>::iterator it = xSections.begin();
  for (; it != xSections.end(); ++it) {
    if( p == *it ) {return;}
  }
  xSections.push_back(p);
}

void G4CrossSectionDataSetRegistry::DeRegister(G4VCrossSectionDataSet* p)
{
  std::vector<G4VCrossSectionDataSet*>::iterator it = xSections.begin();
  for (; it != xSections.end(); ++it) {
    if( p == *it ) {
      xSections.erase(it);
      return;
    }
  }
}


