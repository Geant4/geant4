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
// $Id: G4tgbRotationMatrixMgr.cc 66872 2013-01-15 01:25:57Z japost $
//
//
// class G4tgbRotationMatrixMgr

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgbRotationMatrixMgr.hh"

#include "G4SystemOfUnits.hh"
#include "G4tgrRotationMatrixFactory.hh"
#include "G4tgrMessenger.hh"

// -------------------------------------------------------------------------

G4ThreadLocal G4tgbRotationMatrixMgr * G4tgbRotationMatrixMgr::theInstance = 0;


// -------------------------------------------------------------------------
G4tgbRotationMatrixMgr::G4tgbRotationMatrixMgr()
{
}


// -------------------------------------------------------------------------
G4tgbRotationMatrixMgr* G4tgbRotationMatrixMgr::GetInstance()
{
  if( !theInstance )
  {
    theInstance = new G4tgbRotationMatrixMgr;
    theInstance->CopyRotMats();
  }
  return theInstance;
}


// -------------------------------------------------------------------------
G4tgbRotationMatrixMgr::~G4tgbRotationMatrixMgr()
{
  G4mstgbrotm::const_iterator tgbcite;
  for( tgbcite = theTgbRotMats.begin();
       tgbcite != theTgbRotMats.end(); tgbcite++)
  {
    delete (*tgbcite).second;
  }
  theTgbRotMats.clear();
  delete theInstance;
}


// -------------------------------------------------------------------------
void G4tgbRotationMatrixMgr::CopyRotMats()
{
  G4mstgrrotm tgrRotms =
    G4tgrRotationMatrixFactory::GetInstance()->GetRotMatMap();
  G4mstgrrotm::iterator cite;
  for( cite = tgrRotms.begin(); cite != tgrRotms.end(); cite++ )
  {
    G4tgrRotationMatrix* tgr = (*cite).second;
    G4tgbRotationMatrix* tgb = new G4tgbRotationMatrix( tgr );
    theTgbRotMats[tgb->GetName()] = tgb;
  }
}


// -------------------------------------------------------------------------
G4RotationMatrix*
G4tgbRotationMatrixMgr::FindOrBuildG4RotMatrix(const G4String& name)
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgbRotationMatrixMgr::FindOrBuildG4RotMatrix() - "
           << name << G4endl;
  }
#endif
  G4RotationMatrix* g4rotm = FindG4RotMatrix( name );
  if( g4rotm == 0 )
  {
    G4tgbRotationMatrix* hrotm = FindOrBuildTgbRotMatrix( name );
    // GetRotMatrix() never returns 0, otherwise if not found, it crashes
    g4rotm = hrotm->BuildG4RotMatrix();
  }
  return g4rotm;
}        


// -------------------------------------------------------------------------
G4RotationMatrix* G4tgbRotationMatrixMgr::FindG4RotMatrix(const G4String& name)
{
  G4RotationMatrix* g4rotm = 0;

  G4msg4rotm::const_iterator cite = theG4RotMats.find( name );
  if( cite != theG4RotMats.end() )
  {
    g4rotm = (*cite).second;
  } 

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgbRotationMatrixMgr::FindG4RotMatrix(): " << G4endl
           << "   Name: " << name << " = " << g4rotm << G4endl;
  }
#endif
  
  return g4rotm;
}


// -------------------------------------------------------------------------
G4tgbRotationMatrix*
G4tgbRotationMatrixMgr::FindOrBuildTgbRotMatrix(const G4String& name)
{
  G4tgbRotationMatrix* rotm = FindTgbRotMatrix( name );

  if( rotm == 0 )
  {
    G4String ErrMessage = "Rotation Matrix " + name + " not found !";
    G4Exception("G4tgbRotationMatrixFactory::FindOrBuildRotMatrix()",
                "InvalidSetup", FatalException, ErrMessage); 
  }
  return rotm;
}


// -------------------------------------------------------------------------
G4tgbRotationMatrix*
G4tgbRotationMatrixMgr::FindTgbRotMatrix(const G4String& name)
{
  G4tgbRotationMatrix* rotm = 0;

  G4mstgbrotm::const_iterator cite = theTgbRotMats.find( name );
  if( cite != theTgbRotMats.end() )
  {
    rotm = (*cite).second;
  }
  return rotm;
}


// -------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os , const G4RotationMatrix & rot)
{
  os << "[ " 
     << rot.thetaX()/deg << '\t' << rot.phiX()/deg << '\t'
     << rot.thetaY()/deg << '\t' << rot.phiY()/deg << '\t'
     << rot.thetaZ()/deg << '\t' << rot.phiZ()/deg << " ]"
     << G4endl;
  return os;
}
