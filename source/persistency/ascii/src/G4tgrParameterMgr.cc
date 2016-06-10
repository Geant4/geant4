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
// $Id: G4tgrParameterMgr.cc 66872 2013-01-15 01:25:57Z japost $
//
//
// class G4tgrParameterMgr

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrParameterMgr.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMaterialFactory.hh"
#include "G4tgrRotationMatrixFactory.hh"
#include "G4tgrFileReader.hh"
#include "G4tgrMessenger.hh"
#include "G4UIcommand.hh"

G4ThreadLocal G4tgrParameterMgr* G4tgrParameterMgr::theInstance = 0;


//-------------------------------------------------------------
G4tgrParameterMgr::G4tgrParameterMgr()
{
}


//-------------------------------------------------------------
G4tgrParameterMgr::~G4tgrParameterMgr()
{
  delete theInstance;
}


//-------------------------------------------------------------
G4tgrParameterMgr* G4tgrParameterMgr::GetInstance()
{
  if( !theInstance )
  {
    theInstance = new G4tgrParameterMgr;
  }
  return theInstance;
}


//-------------------------------------------------------------
void G4tgrParameterMgr::AddParameterNumber( const std::vector<G4String>& wl,
                                                  G4bool mustBeNew  )
{
  CheckIfNewParameter( wl, mustBeNew );

  //----- Convert third argument to double, but then store it as string
  //      for later use in CLHEP evaluator 
  G4float val = G4tgrUtils::GetDouble( wl[2] );
  theParameterList[ wl[1] ] = G4UIcommand::ConvertToString( val );

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgrParameterMgr::AddParameterNumber() -"
           << " parameter added " <<  wl[1]
           << " = " << theParameterList[ wl[1] ] << G4endl; 
  }
#endif
}


//-------------------------------------------------------------
void G4tgrParameterMgr::AddParameterString( const std::vector<G4String>& wl,
                                                  G4bool mustBeNew  )
{
  CheckIfNewParameter( wl, mustBeNew );

  //----- Store parameter as string
  theParameterList[ wl[1] ] = wl[2];

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgrParameterMgr::AddParameterString() -"
           << " parameter added " <<  wl[1]
           << " = " << theParameterList[ wl[1] ] << G4endl; 
  }
#endif
}

//-------------------------------------------------------------
void G4tgrParameterMgr::CheckIfNewParameter( const std::vector<G4String>& wl,
                                                   G4bool mustBeNew  )
{
  //---------- Find first if it exists already
  G4bool existsAlready;
  G4mapss::iterator sdite = theParameterList.find( wl[1] );
  if( sdite == theParameterList.end() )
  {
    existsAlready = false;
  }
  else
  {
    existsAlready = true;
  }

  if(existsAlready)
  {
    if( mustBeNew )
    {
      G4String ErrMessage = "Parameter already exists... " + wl[1];
      G4Exception("G4tgrParameterMgr::CheckParameter()",
                  "IllegalConstruct", FatalException, ErrMessage);
    }
    else
    {
      G4String WarMessage = "Parameter already exists... " + wl[1];
      G4Exception("G4tgrParameterMgr::CheckParameter()",
                  "NotRecommended", JustWarning, WarMessage);
    }
  }
  
  //---------- Check for miminum number of words read 
  G4tgrUtils::CheckWLsize( wl, 3, WLSIZE_EQ, "Parameter::AddParameter");
}


//-------------------------------------------------------------
G4String G4tgrParameterMgr::FindParameter( const G4String& name, G4bool exists )
{
  G4String par = "";
  
  G4mapss::iterator sdite = theParameterList.find( name );
  if( sdite == theParameterList.end() )
  {
    if( exists )
    {
      DumpList();
      G4String ErrMessage = "Parameter not found in list: " + name;
      G4Exception("G4tgrParameterMgr::FindParameter()",
                  "InvalidSetup", FatalException, ErrMessage);
    }
  }
  else
  {
    exists = 1;
    par = ((*sdite).second);
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 3 )
  {
    G4cout << " G4tgrParameterMgr::FindParameter() -"
           << " parameter found " << name << " = " << par << G4endl; 
  }
#endif
  }

  return par;
}


//-------------------------------------------------------------
void G4tgrParameterMgr::DumpList()
{
  //---------- Dump number of objects of each class
  G4cout << " @@@@@@@@@@@@@@@@@@ Dumping parameter list " << G4endl;
  G4mapss::const_iterator cite;
  for( cite = theParameterList.begin();cite != theParameterList.end(); cite++ )
  {
    G4cout << (*cite).first << " = " << (*cite).second << G4endl;
  }
}
