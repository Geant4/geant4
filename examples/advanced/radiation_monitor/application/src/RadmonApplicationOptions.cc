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
// File name:     RadmonApplicationOptions.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationOptions.cc,v 1.5 2006/06/29 16:08:43 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//

// Include files
#include "RadmonApplicationOptions.hh"





                                                RadmonApplicationOptions :: RadmonApplicationOptions(int argc, char * * argv)
:
 valid(true),
 interactive(true),
 verbose(false),
 help(false),
 applicationName(argv[0]),
 fileName(0),
 startupFileName(".startup.mac")
{
 int i(1);
 G4bool fileNameGiven(false);
 
 while (i<argc)
 {
  if (argv[i][0]=='-' && argv[i][2]==0)
  {
   switch(argv[i][1])
   {
    case '?':
    case 'H':
    case 'h':
     help=true;
     return;
    
    case 'b':
    case 'B':
     interactive=false;
     break;
    
    case 'v':
    case 'V':
     verbose=true;
     break;
    
    default:
     G4cerr << applicationName << ": Command " << argv[i] << " unknown. Run \"" << applicationName << " -h\" for help." << G4endl;
     valid=false;
     return;
   }
  }
  else
  {
   if (fileNameGiven)
   {
    G4cerr << applicationName << ": Bad syntax. Run \"" << applicationName << " -h\" for help." << G4endl;
    valid=false;
    return;
   }
   
   fileNameGiven=true;
   fileName=argv[i];
  }  
 
  i++;
 }
}



void                                            RadmonApplicationOptions :: DumpHelp(void) const
{
 G4cout << "Usage: " << applicationName << " [-h|-H|-?] [-b] [-v] [<filename>]" << G4endl;
 G4cout << G4endl;
 G4cout << "-h  -H  -?         Usage help" << G4endl;
 G4cout << "-b                 Force non-interactive mode" << G4endl;
 G4cout << "-v                 Verbose output" << G4endl;
 G4cout << "<filename>         Name of the macro to be run" << G4endl;
 G4cout << G4endl;
 G4cout << "If \"" << startupFileName << "\" is present, it will be run before anything else." << G4endl;
 G4cout << G4endl;
}

