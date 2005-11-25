//
// File name:     RadmonApplicationOptions.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationOptions.cc,v 1.3 2005-11-25 01:56:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
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

