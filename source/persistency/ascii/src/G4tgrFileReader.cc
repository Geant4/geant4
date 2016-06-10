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
// $Id: G4tgrFileReader.cc 66872 2013-01-15 01:25:57Z japost $
//
//
// class G4tgrFileReader

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrFileReader.hh"
#include "G4tgrParameterMgr.hh"
#include "G4tgrFileIn.hh"
#include "G4tgrElementSimple.hh"
#include "G4tgrElementFromIsotopes.hh"
#include "G4tgrVolume.hh"
#include "G4tgrPlaceDivRep.hh"
#include "G4tgrPlaceParameterisation.hh"
#include "G4tgrVolumeDivision.hh"
#include "G4tgrVolumeMgr.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMaterialFactory.hh"
#include "G4tgrRotationMatrixFactory.hh"
#include "G4tgrLineProcessor.hh"
#include "G4tgrMessenger.hh"


G4ThreadLocal G4tgrFileReader* G4tgrFileReader::theInstance = 0;


//---------------------------------------------------------------
G4tgrFileReader::G4tgrFileReader()
{
  theLineProcessor = new G4tgrLineProcessor;
}


//---------------------------------------------------------------
G4tgrFileReader::~G4tgrFileReader()
{
  delete theLineProcessor;
  delete theInstance;
}


//---------------------------------------------------------------
G4tgrFileReader* G4tgrFileReader::GetInstance()
{
  if( !theInstance ) {
    theInstance = new G4tgrFileReader;
  }
  return theInstance;
}

//---------------------------------------------------------------
G4bool G4tgrFileReader::ReadFiles() 
{

  std::vector< G4String > wl,wlnew;
    
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << "   Number of geometry data files = " << theTextFiles.size()
           << G4endl;
  }
#endif

  if( theTextFiles.size() == 0 )
  {
    G4Exception("G4tgrFileReader::ReadFiles()", "InvalidInput",
                FatalException, "No files to read ...");
  }

  for( size_t ii = 0; ii < theTextFiles.size(); ii++ )
  {
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 1 )
    {
      G4cout << "   Reading data file " << theTextFiles[ii] << G4endl;
    }
#endif
    
    G4tgrFileIn fin = G4tgrFileIn::GetInstance( theTextFiles[ii] );
    
    G4int nlines = 0;
    for(;;)
    {
      nlines++;
      if(! fin.GetWordsInLine( wlnew ) )  { break; }
      // Check if it is continuation line or first line
      if( wlnew[0].c_str()[0] != ':' )
      {
        wl.insert( wl.end(), wlnew.begin(), wlnew.end() );
#ifdef G4VERBOSE
        if( G4tgrMessenger::GetVerboseLevel() >= 4 )
        {
          G4tgrUtils::DumpVS(wl, "!!!! adding line");
        }
#endif
        continue;
      }
      else
      {
        //----- Process previous tag
#ifdef G4VERBOSE
        if( G4tgrMessenger::GetVerboseLevel() >= 4 )
        {
          G4tgrUtils::DumpVS(wl, "!!!! line read");
        }
#endif
        if( nlines != 1)   // first line has no previous tag
        {
          if( ! theLineProcessor->ProcessLine( wl ) )
          {
            fin.DumpException( "Tag not found: " + wl[0]);
          }
        }
        wl = wlnew;
      }
    }
    
    if( wl.size() != 0 )
    {
      if( ! theLineProcessor->ProcessLine( wl ) )
      {
        fin.DumpException( "Tag not found: " + wl[0]);
      }
    }
  }  
  return 1;
}
