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
#ifndef G4CrossSectionFactory_h
#define G4CrossSectionFactory_h 1


#include "globals.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4VCrossSectionDataSet.hh"

class G4VBaseXSFactory
{

public:

  virtual G4VCrossSectionDataSet* Instantiate() = 0;

};


template <typename T> class G4CrossSectionFactory : public G4VBaseXSFactory
{
public:
  
  G4CrossSectionFactory(const G4String& name)
  {
    G4CrossSectionDataSetRegistry::Instance()->AddFactory(name, this);
  }
  
  virtual G4VCrossSectionDataSet* Instantiate() 
  {
    return new T();
  }
};


#define G4_DECLARE_XS_FACTORY(cross_section) \
  const G4CrossSectionFactory<cross_section>& cross_section##Factory = G4CrossSectionFactory<cross_section>(cross_section::Default_Name())

#define G4_REFERENCE_XS_FACTORY(cross_section) \
  class cross_section; \
  extern const G4CrossSectionFactory<cross_section>& cross_section##Factory; \
  const G4CrossSectionFactory<cross_section>& cross_section##FactoryRef = cross_section##Factory


#endif
