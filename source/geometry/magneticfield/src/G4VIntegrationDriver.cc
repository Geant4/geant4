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
// $Id: G4VIntegrationDriver.cc 109569 2018-05-02 07:08:33Z gcosmo $
//
// class G4VIntegrationDriver
//
// Class description:
//
// Abstract base class for 'driver' classes which are responsible for
// undertaking integration of an state given an equation of motion and
// within acceptable error bound(s).
//
//  Different integration methods are meant to be provided via this
// common interface, and can span the original type (explicit Runge Kutta
// methods), enhanced RK methods and alternatives such as the
// Bulirsch-Stoer and multi-step methods.
//
//  The drivers' key mission is to insure that the error is below set values. 
//
//    Implementation by Dmitry Sorokin - GSoC 2017
//       Work supported by Google as part of Google Summer of Code 2017.
//    Supervision / code review: John Apostolakis

#include "G4VIntegrationDriver.hh"

void G4VIntegrationDriver::RenewStepperAndAdjust(G4MagIntegratorStepper *)
{
   G4Exception("G4VIntegrationDriver::RenewStepperAndAdjust", "Geometry001", FatalException,
               "This method exists only for the original G4MagIntegratorDriver class.  "
               " Not defined for other classes derived from G4VIntegrationDriver");
}
