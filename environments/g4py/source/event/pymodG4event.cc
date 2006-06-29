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
// $Id: pymodG4event.cc,v 1.4 2006-06-29 15:31:47 gunter Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pymodG4event.cc [Geant4Py module]
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4EventManager();
void export_G4StackManager();
void export_G4Event();
void export_G4UserEventAction();
void export_G4UserStackingAction();
void export_G4ClassificationOfNewTrack();
void export_G4ParticleGun();

BOOST_PYTHON_MODULE(G4event)
{
  export_G4EventManager();
  export_G4StackManager();
  export_G4Event();
  export_G4UserEventAction();
  export_G4UserStackingAction();
  export_G4ClassificationOfNewTrack();
  export_G4ParticleGun();
}

