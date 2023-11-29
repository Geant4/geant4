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
// G4QSSDriverCreator
//
// QSS driver creator

// Author: J.Apostolakis - October 2021
// -------------------------------------------------------------------
#ifndef QSS_DriverCreator_HH
#define QSS_DriverCreator_HH

#include "G4QSStepper.hh"
#include "G4QSSDriver.hh"

class G4Mag_EqRhs;
class G4QSS2;
class G4QSS3;

class G4QSSDriverCreator
{ 
  public:
         
    static G4VIntegrationDriver* CreateDriver( G4MagIntegratorStepper* pStepper, 
                                               G4double /*stepMinimum*/ );

    static G4QSSDriver<G4QSStepper<G4QSS2>>*
    CreateDriver( G4QSStepper<G4QSS2>* qss2stepper );

    static G4QSSDriver<G4QSStepper<G4QSS3>>*
    CreateDriver( G4QSStepper<G4QSS3>* qss3stepper );

    static G4QSStepper<G4QSS2>* CreateQss2Stepper(G4Mag_EqRhs* Equation);

    static G4QSStepper<G4QSS3>* CreateQss3Stepper(G4Mag_EqRhs* Equation);

    static G4VIntegrationDriver* CreateQss2Driver(G4Mag_EqRhs* Equation);

    static G4VIntegrationDriver* CreateQss3Driver(G4Mag_EqRhs* Equation);
};

#endif
