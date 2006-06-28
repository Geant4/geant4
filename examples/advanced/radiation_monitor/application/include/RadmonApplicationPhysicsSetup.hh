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
// File name:     RadmonApplicationPhysicsSetup.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationPhysicsSetup.hh,v 1.2 2006-06-28 13:45:24 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
