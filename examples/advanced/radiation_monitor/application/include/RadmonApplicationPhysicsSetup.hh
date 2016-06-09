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
// File name:     RadmonApplicationPhysicsSetup.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationPhysicsSetup.hh,v 1.3 2006/06/29 16:08:21 gunter Exp $
// Tag:           $Name: geant4-09-02 $
//
// Description:   Radmon application generators algoritms setup
//

#ifndef   RADMONAPPLICATIONPHYSICSSETUP_HH
 #define  RADMONAPPLICATIONPHYSICSSETUP_HH

 // Inglude files
 #include "globals.hh"

 // Forward declarations
 class RadmonApplicationOptions;
 class RadmonSubPhysicsListWithLabelFactory;

 class RadmonApplicationPhysicsSetup
 {
  public:
   inline                                       RadmonApplicationPhysicsSetup(const RadmonApplicationOptions & options);
   inline                                      ~RadmonApplicationPhysicsSetup();

  protected:
   G4bool                                       CreateSubPhysicsList(RadmonSubPhysicsListWithLabelFactory * factory);

  // Hidden constructors and operators
                                                RadmonApplicationPhysicsSetup();
                                                RadmonApplicationPhysicsSetup(const RadmonApplicationPhysicsSetup & copy);
   RadmonApplicationPhysicsSetup &              operator=(const RadmonApplicationPhysicsSetup & copy);
   
  // Private attributes
   const RadmonApplicationOptions &             currentOptions;
 };
 
 #include "RadmonApplicationPhysicsSetup.icc"
#endif /* RADMONAPPLICATIONPHYSICSSETUP_HH */
