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
// File name:     Radmon.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: Radmon.cc,v 1.2 2006-06-28 13:43:45 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Application main file
//

// Include files
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "RadmonApplicationOptions.hh"
#include "RadmonApplication.hh"


// Main
int main(int argc, char * * argv)
{ 
 RadmonApplicationOptions options(argc, argv);

 if (!options.Valid())
  return EXIT_FAILURE;
 
 if (options.Help())
 {
  options.DumpHelp();
  return EXIT_SUCCESS;
 }
 
 RadmonApplication application(options);
 
 if (!application.Valid())
  return EXIT_FAILURE;
  
 G4cout << options.ApplicationName() << ": Application ended successfully." << G4endl;
  
 return EXIT_SUCCESS;
}
