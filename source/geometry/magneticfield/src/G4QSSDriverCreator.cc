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
// G4QSSDriverCreator implementation
//
// Author: J.Apostolakis - October 2021 - February 2023 
// --------------------------------------------------------------------

#include "G4QSSDriverCreator.hh"

#include "G4MagIntegratorStepper.hh"
#include "G4VIntegrationDriver.hh"

#include "G4QSSDriver.hh"
#include "G4QSStepper.hh"
#include "G4QSS2.hh"
#include "G4QSS3.hh"

#include "G4Mag_UsualEqRhs.hh"

#include <cassert>

G4VIntegrationDriver* 
G4QSSDriverCreator::CreateDriver( G4MagIntegratorStepper* pStepper, G4double /*stepMin*/ )
{
  G4VIntegrationDriver* driver = nullptr;
//   pStepper->build_driver(stepMinimum, true); // Original - QSS
  auto qss2stepper = dynamic_cast<G4QSStepper<G4QSS2>*>(pStepper);
  if( qss2stepper != nullptr ) {
    // driver = new G4QSSDriver<G4QSStepper<G4QSS2>>(qss2stepper);
    driver = CreateDriver( qss2stepper );
  }
  auto qss3stepper = dynamic_cast<G4QSStepper<G4QSS3>*>(pStepper);
  if( qss3stepper != nullptr ) {
    // driver = new G4QSSDriver<G4QSStepper<G4QSS3>>(qss3stepper);
    driver= CreateDriver( qss3stepper );          
  }
  return driver;
}

G4QSSDriver<G4QSStepper<G4QSS2>>*
    G4QSSDriverCreator::CreateDriver( G4QSStepper<G4QSS2>* qss2stepper )
{
  G4cout << "---- G4QSSDriver<G4QSS2>* G4QSSDriverCreator::CreateDriver(G4QSStepper<G4QSS2>* ) called.\n";
  return new G4QSSDriver<G4QSStepper<G4QSS2>>(qss2stepper);
}

static constexpr G4int numOfVars= 6;
   
G4QSSDriver<G4QSStepper<G4QSS3>>*
   G4QSSDriverCreator::CreateDriver( G4QSStepper<G4QSS3>* qss3stepper )
{
  G4cout << "---- G4QSSDriver<G4QSS3>* G4QSSDriverCreator::CreateDriver(G4QSStepper<G4QSS3>* ) called.\n";        
  return new G4QSSDriver<G4QSStepper<G4QSS3>>(qss3stepper);
}

G4QSStepper<G4QSS2>* G4QSSDriverCreator::G4QSSDriverCreator::CreateQss2Stepper(G4Mag_EqRhs* Equation)
{
  G4cout << "---- G4QSStepper<G4QSS2>* CreateQss2Stepper(G4Mag_EqRhs* ) CALLED\n";
  return G4QSStepper<G4QSS2>::build_QSS2( Equation, numOfVars, true);
}

G4QSStepper<G4QSS3>* G4QSSDriverCreator::CreateQss3Stepper(G4Mag_EqRhs* Equation)
{
  G4cout << "---- G4QSStepper<G4QSS3>* CreateQss3Stepper(G4Mag_EqRhs* ) CALLED\n";
  return G4QSStepper<G4QSS3>::build_QSS3( Equation, numOfVars, true);
}

G4VIntegrationDriver* G4QSSDriverCreator::CreateQss2Driver(G4Mag_EqRhs* Equation)
{
  assert( dynamic_cast<G4Mag_UsualEqRhs*>(Equation) != nullptr );
  // assert( Equation->GetNumberOfVariables() == numOfVars );

  auto qss2stepper = G4QSStepper<G4QSS2>::build_QSS2( Equation, numOfVars, true);
  return CreateDriver( qss2stepper );
}

G4VIntegrationDriver* G4QSSDriverCreator::
CreateQss3Driver(G4Mag_EqRhs *Equation)
{
  assert( dynamic_cast<G4Mag_UsualEqRhs*>(Equation) != nullptr );        
  // assert( Equation->GetNumberOfVariables() == numOfVars );

  auto qss3stepper = G4QSStepper<G4QSS3>::build_QSS3( Equation, numOfVars, true);
  return CreateDriver( qss3stepper );
}
