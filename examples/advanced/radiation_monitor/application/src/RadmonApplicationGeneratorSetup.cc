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
// File name:     RadmonApplicationGeneratorSetup.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationGeneratorSetup.cc,v 1.2 2006-06-28 13:45:52 gunter Exp $
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
