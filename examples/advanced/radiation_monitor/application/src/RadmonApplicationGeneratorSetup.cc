//
// File name:     RadmonApplicationGeneratorSetup.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationGeneratorSetup.cc,v 1.1 2005-10-25 16:39:12 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonApplicationGeneratorSetup.hh"
#include "RadmonApplicationOptions.hh"

#include "RadmonGeneratorFixedPosition.hh"
#include "RadmonGeneratorFixedDirection.hh"
#include "RadmonGeneratorFixedEnergy.hh"
#include "RadmonGeneratorFixedParticle.hh"
#include "RadmonGeneratorUniformSphere.hh"
#include "RadmonGeneratorUniformPlane.hh"


#include "RadmonGeneratorsWithLabelFactory.hh"


#define DECLARE_GENERATOR_CONSTRUCTOR(name)     constructor=new name();                                                                  \
                                                if (constructor==0)                                                                      \
                                                {                                                                                        \
                                                 G4cerr << currentOptions.ApplicationName() << ": Cannot allocate " #name "." << G4endl; \
                                                 return false;                                                                           \
                                                }                                                                                        \
                                                factory->AppendGenerator(constructor)

G4bool                                          RadmonApplicationGeneratorSetup :: CreateGenerators(RadmonGeneratorsWithLabelFactory * factory)
{
 RadmonVGeneratorWithLabel * constructor;

 DECLARE_GENERATOR_CONSTRUCTOR(RadmonGeneratorFixedPosition);
 DECLARE_GENERATOR_CONSTRUCTOR(RadmonGeneratorFixedDirection);
 DECLARE_GENERATOR_CONSTRUCTOR(RadmonGeneratorFixedEnergy);
 DECLARE_GENERATOR_CONSTRUCTOR(RadmonGeneratorFixedParticle);
 DECLARE_GENERATOR_CONSTRUCTOR(RadmonGeneratorUniformSphere);
 DECLARE_GENERATOR_CONSTRUCTOR(RadmonGeneratorUniformPlane);
 
 return true;
}
