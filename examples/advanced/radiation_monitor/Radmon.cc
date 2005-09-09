//
// File name:     Radmon.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: Radmon.cc,v 1.1 2005-09-09 08:27:34 capra Exp $
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
