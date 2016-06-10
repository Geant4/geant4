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
// $Id: G4tgrRotationMatrixFactory.cc 66872 2013-01-15 01:25:57Z japost $
//
//
// class G4tgrRotationMatrixFactory

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrRotationMatrixFactory.hh"
#include "G4tgrMessenger.hh"
#include "G4tgrUtils.hh"


G4ThreadLocal G4tgrRotationMatrixFactory * G4tgrRotationMatrixFactory::theInstance = 0;


// -------------------------------------------------------------------------
G4tgrRotationMatrixFactory* G4tgrRotationMatrixFactory::GetInstance()
{
  if( !theInstance )
  {
    theInstance = new G4tgrRotationMatrixFactory;
  }
  return theInstance;
}


// -------------------------------------------------------------------------
G4tgrRotationMatrixFactory::G4tgrRotationMatrixFactory()
{
}


// -------------------------------------------------------------------------
G4tgrRotationMatrixFactory::~G4tgrRotationMatrixFactory()
{
  G4mstgrrotm::iterator cite;
  for( cite = theTgrRotMats.begin(); cite != theTgrRotMats.end(); cite++)
  {
    delete (*cite).second;
  }
  theTgrRotMats.clear();
  delete theInstance;
}


// -------------------------------------------------------------------------
G4tgrRotationMatrix*
G4tgrRotationMatrixFactory::AddRotMatrix( const std::vector<G4String>& wl )
{
  //---------- Check for miminum number of words read 
  if( wl.size() != 5 && wl.size() != 8 && wl.size() != 11 )
  {
    G4tgrUtils::DumpVS(wl, "G4tgrRotationMatrixFactory::AddRotMatrix()");
    G4Exception("G4tgrRotationMatrixFactory::AddRotMatrix()", "InvalidMatrix",
                FatalException, "Line should have 5, 8 or 11 words !");
  }

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgrRotationMatrixFactory::AddRotMatrix() - Adding: "
           << wl[1] << G4endl;
  }
#endif
  //---------- Look if rotation matrix exists
  if( FindRotMatrix( G4tgrUtils::GetString(wl[1]) ) != 0 )
  {
    G4String ErrMessage = "Rotation matrix repeated... " + wl[1];
    G4Exception("G4tgrRotationMatrixFactory::AddRotMatrix()",
                "InvalidInput", FatalException, ErrMessage);
  } 
 
  G4tgrRotationMatrix* rotm = new G4tgrRotationMatrix( wl );
  theTgrRotMats[ rotm->GetName() ] =  rotm;
  theTgrRotMatList.push_back( rotm );
 
  return rotm;
}


// -------------------------------------------------------------------------
G4tgrRotationMatrix*
G4tgrRotationMatrixFactory::FindRotMatrix(const G4String& name)
{
  G4tgrRotationMatrix* rotm = 0;

  G4mstgrrotm::const_iterator cite = theTgrRotMats.find( name );
  if( cite != theTgrRotMats.end() )
  { 
    rotm = (*cite).second;
  } 

  return rotm;
}


// -------------------------------------------------------------------------
void G4tgrRotationMatrixFactory::DumpRotmList()
{
  G4cout << " @@@@@@@@@@@@@@@@ DUMPING G4tgrRotationMatrix's List " << G4endl;
  G4mstgrrotm::const_iterator cite;
  for(cite = theTgrRotMats.begin(); cite != theTgrRotMats.end(); cite++)
  {
    G4cout << " ROTM: " << (*cite).second->GetName() << G4endl;
  }
}
