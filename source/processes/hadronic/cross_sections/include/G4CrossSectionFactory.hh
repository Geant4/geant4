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
#include "G4VCrossSectionDataSet.hh"
#include "G4CrossSectionFactoryRegistry.hh"
#include "G4Threading.hh"

class G4VBaseXSFactory
{

public:

  virtual G4VCrossSectionDataSet* Instantiate() = 0;

};


//Generic template XS-factory
template <typename T, int mode> class G4CrossSectionFactory : public G4VBaseXSFactory
{
public:
  G4CrossSectionFactory(const G4String& name)
  {
      G4CrossSectionFactoryRegistry::Instance()->Register(name,this);
  }
  
  virtual G4VCrossSectionDataSet* Instantiate() 
  {
      G4ExceptionDescription msg;
      msg<<"Factory mode: "<<mode<<" not supported!";
      G4Exception("G4CrossSectionFactory::Instantiate","CrossSectionFactory001",FatalException,msg);
    return static_cast<T*>(0);
  }
};

//Partial specialized template for non-singleton non-shared factory
// each call to Instantiate creates a new XS
template <typename T> class G4CrossSectionFactory<T,0> : public G4VBaseXSFactory
{
public:
    
    G4CrossSectionFactory(const G4String& name)
    {
        G4CrossSectionFactoryRegistry::Instance()->Register(name,this);
    }
    
    virtual G4VCrossSectionDataSet* Instantiate()
    {
        return new T();
    }
};

//Partial specialized template for singleton, shared factory
// each call to Instantiate returns pointer to static object
template <typename T> class G4CrossSectionFactory<T,1> : public G4VBaseXSFactory
{
public:
    G4CrossSectionFactory(const G4String& name)
    {
        G4CrossSectionFactoryRegistry::Instance()->Register(name,this);
    }
    
    virtual G4VCrossSectionDataSet* Instantiate()
    {
        static T* shared = new T();
        return shared;
    }
};

//Partial specialized template for singleton, shared factory
// each call to Instantiate returns pointer to static thread-local object
template <typename T> class G4CrossSectionFactory<T,2> : public G4VBaseXSFactory
{
    G4CrossSectionFactory(const G4String& name)
    {
        G4CrossSectionFactoryRegistry::Instance()->Register(name,this);
    }
    
    virtual G4VCrossSectionDataSet* Instantiate()
    {
        static G4ThreadLocal T* shared = 0;
        if (!shared) { shared = new T(); }
        return shared;
    }
};


#define G4_BASE_DECLARE_XS_FACTORY(cross_section, flag) \
  const G4CrossSectionFactory<cross_section,flag>& cross_section##Factory = G4CrossSectionFactory<cross_section,flag>(cross_section::Default_Name())

#define G4_BASE_REFERENCE_XS_FACTORY(cross_section,flag) \
  class cross_section; \
  extern const G4CrossSectionFactory<cross_section,flag>& cross_section##Factory; \
  const G4CrossSectionFactory<cross_section,flag>& cross_section##FactoryRef = cross_section##Factory

//Macros to help define and reference factories
#define G4_DECLARE_XS_FACTORY(cross_section) G4_BASE_DECLARE_XS_FACTORY(cross_section,0)
#define G4_DECLARE_SHAREDXS_FACTORY(cross_section) G4_BASE_DECLARE_XS_FACTORY(cross_section,1)
#define G4_DECLARE_SHAREDTLSXS_FACTORY(cross_section) G4_BASE_DECLARE_XS_FACTORY(cross_section,2)


#define G4_REFERENCE_XS_FACTORY(cross_section) G4_BASE_REFERENCE_XS_FACTORY(cross_section,0)
#define G4_REFERENCE_SHAREDXS_FACTORY(cross_section) G4_BASE_REFERENCE_XS_FACTORY(cross_section,1)
#define G4_REFERENCE_SHAREDTLSXS_FACTORY(cross_section) G4_BASE_REFERENCE_XS_FACTORY(cross_section,2)

#endif
