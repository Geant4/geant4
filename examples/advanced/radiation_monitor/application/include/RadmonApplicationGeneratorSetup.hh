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
// File name:     RadmonApplicationGeneratorSetup.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationGeneratorSetup.hh,v 1.2 2006-06-28 13:45:14 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Radmon application generators algoritms setup
//

#ifndef   RADMONAPPLICATIONGENERATORSETUP_HH
 #define  RADMONAPPLICATIONGENERATORSETUP_HH

 // Inglude files
 #include "globals.hh"

 // Forward declarations
 class RadmonApplicationOptions;
 class RadmonGeneratorsWithLabelFactory;

 class RadmonApplicationGeneratorSetup
 {
  public:
   inline                                       RadmonApplicationGeneratorSetup(const RadmonApplicationOptions & options);
   inline                                      ~RadmonApplicationGeneratorSetup();

  protected:
   G4bool                                       CreateGenerators(RadmonGeneratorsWithLabelFactory * factory);

  // Hidden constructors and operators
                                                RadmonApplicationGeneratorSetup();
                                                RadmonApplicationGeneratorSetup(const RadmonApplicationGeneratorSetup & copy);
   RadmonApplicationGeneratorSetup &            operator=(const RadmonApplicationGeneratorSetup & copy);
   
  // Private attributes
   const RadmonApplicationOptions &             currentOptions;
 };
 
 #include "RadmonApplicationGeneratorSetup.icc"
#endif /* RADMONAPPLICATIONGENERATORSETUP_HH */
