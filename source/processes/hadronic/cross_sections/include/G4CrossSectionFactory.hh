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
  static const G4CrossSectionFactory<cross_section> my_factory(cross_section::Default_Name())
#endif
