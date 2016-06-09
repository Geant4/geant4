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
/// \file geometry/olap/src/OlapManager.cc
/// \brief Implementation of the OlapManager class
//
//
// $Id$
//
// 
// --------------------------------------------------------------
// OlapManager
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#include <stdlib.h>

#include <vector>

#include "OlapManager.hh"
#include "OlapDetConstr.hh"
#include "OlapManagerMessenger.hh"
#include "OlapGenerator.hh"
#include "OlapEventAction.hh"
#include "G4GeoNav.hh"
#include "G4GeometryTolerance.hh"
#include "G4RunManager.hh"
#include "globals.hh"

OlapManager * OlapManager::theInstance = 0;

OlapManager::OlapManager()
  : olapGenerator(0), eventsPerRun(27), lvPos(0),
    polyMode(false), interrupt(false) 
{
   theMessenger = new OlapManagerMessenger(this);
   theLVStore = G4LogicalVolumeStore::GetInstance();

   delta = G4GeometryTolerance::GetInstance()->GetAngularTolerance();
   theRunManager = G4RunManager::GetRunManager();
   const G4VUserDetectorConstruction * aTmp =
      theRunManager->GetUserDetectorConstruction();  
   
   const OlapDetConstr * aConstDetConstr =
     dynamic_cast<const OlapDetConstr*>(aTmp);
   if (!aConstDetConstr)
   {
     G4Exception("OlapManager::OlapManager()", "InvalidSetup", FatalException,
                 "Can't get the OlapDetConstr instance! Exiting...");
   }  
 
   // need a  non-const instance for building NewWords!          
   theDet = const_cast<OlapDetConstr*>(aConstDetConstr);
   
   // instantiate the logical volume tree & add it to the Gui
   theGeoNav = new G4GeoNav(theDet->GetFullWorld()->GetLogicalVolume());
   
   std::vector<G4LogicalVolume*>::iterator aIt2;
   for (aIt2=theLVStore->begin(); aIt2 != theLVStore->end(); aIt2++)
     NoOlapMap[*aIt2] = false;
}


OlapManager::~OlapManager()
{
   delete theMessenger;
   delete theGeoNav;
   delete olapGenerator;
}


OlapManager * OlapManager::GetOlapManager()
{
   if (!theInstance)
      theInstance = new OlapManager();
   
   return theInstance;   
}   


void OlapManager::SetRotation(G4double theta, G4double phi, G4double alpha)
{
  theDet->SetRotation(theta,phi,alpha);
}



void OlapManager::TriggerRun()
{
   const OlapGenerator * aGen = dynamic_cast<const OlapGenerator *>
         (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
   OlapGenerator * aNonConstGen = 0;
   if ( aGen )
   {
        aNonConstGen = const_cast<OlapGenerator*>(aGen);
   }
   else
   {
      G4cout << "Warning: Primary generator is not OlapGenerator!" << G4endl
             << "         Overlap Detection will not work!" << G4endl;
      return;
   }             
   //while
   aNonConstGen->SetAutoIncrement(true); 
   theRunManager->BeamOn(eventsPerRun);
   aNonConstGen->SetAutoIncrement(false);
}


void OlapManager::TriggerFull(G4int trg)
{
   const OlapGenerator * aGen = dynamic_cast<const OlapGenerator *>
         (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
   OlapGenerator * aNonConstGen = 0;
   if ( aGen )
   {
        aNonConstGen = const_cast<OlapGenerator*>(aGen);
   }
   else
   {
      G4cout << "Warning: Primary generator is not OlapGenerator!" << G4endl
             << "         Overlap Detection will not work!" << G4endl;
      return;
   }             
   //while
   aNonConstGen->SetAutoIncrement(true); 
   
   interrupt = false; // can be changed to 'true' in the following loop
                      // by subclasses of OlapNotify. In that case,
                      // the triggerFull-command will be stopped.
   
   while (Next() && trg)
   {
      if (theDet->GetNewWorld()->GetLogicalVolume()->GetDaughter(0)->
          GetLogicalVolume()->GetNoDaughters())
      {
         theRunManager->BeamOn(eventsPerRun);   
         trg--;
         if (interrupt)
         {
           break;
         }
      }          
   }
   aNonConstGen->SetAutoIncrement(false);
}

G4VPhysicalVolume * OlapManager::GetNewWorld() 
{
  return theDet->GetNewWorld();
}        


G4VPhysicalVolume * OlapManager::GetFullWorld() 
{
  return theDet->GetFullWorld();
}        

G4bool OlapManager::Next()
{  
   G4LogicalVolume * nextlv = theGeoNav->NextLV();
   if(nextlv)
   {
      theDet->SetNewWorld(nextlv);
      G4cout << nextlv->GetName() << G4endl;
      return true;
   }
   return false;
}

void OlapManager::ListLV(const G4String & aRegexStr)
{
    G4cout << "Logical volumes matching " << aRegexStr << ":" <<G4endl;
    
    std::vector<G4LogicalVolume *> aLVVec;
    G4int c = theGeoNav->FilterLV(aRegexStr,aLVVec);
 
    for(G4int i=0; i<c; i++)
    {
        G4cout << " " << (aLVVec[i])->GetName() << G4endl;
    }
}


void OlapManager::LsLV()
{
   std::vector<G4LogicalVolume *> lvs;
   G4int c = theGeoNav->LsLV(lvs);
   for (G4int i = 0; i<c; i++)
      G4cout << "   " << (lvs[i])->GetName() << G4endl;
}

void OlapManager::PwdLV()
{
    std::vector<G4LogicalVolume *> lvs;
    theGeoNav->PwdLV(lvs);
    G4String temp;
    G4cout << "/ = ";
    for (G4int i=lvs.size()-1; i>=0; i--)
    {
      G4cout << temp << (lvs[i])->GetName() << G4endl;
      temp = temp + G4String("   ");
    }
}


G4int OlapManager::GetNrLVs()
{ 
   return theDet->GetNrLVs(); 
}


void OlapManager::ChangeLV(const G4String & aDir)
{
   G4LogicalVolume * aLv = theGeoNav->ChangeLV(aDir);
   if (aLv)
   {
       G4cout << "New world: " << aLv->GetName() << G4endl;
       theDet->SetNewWorld(aLv);
       notifyNewWorld(theDet->GetNewWorld()->GetLogicalVolume());
   }
}

void OlapManager::GotoLV(const G4String & aRegexStr)
{
    std::vector<G4LogicalVolume*> lvs; 
    if (theGeoNav->FilterLV(aRegexStr,lvs,true))
    {
       G4cout << "New world: " << (lvs[0])->GetName() << G4endl;
       theDet->SetNewWorld((lvs[0]));
       notifyNewWorld(theDet->GetNewWorld()->GetLogicalVolume());
    }   
}


void OlapManager::SetNewWorld(G4LogicalVolume * nw)
{
   if (nw)
   {
     theDet->SetNewWorld(nw);
     notifyNewWorld(theDet->GetNewWorld()->GetLogicalVolume());
   }  
}

void OlapManager::SetGrid(G4int x, G4int y, G4int z)
{
   // try to find the OlapGenerator and set the new grid
   const OlapGenerator * aGen = dynamic_cast<const OlapGenerator *>
         (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
   OlapGenerator * aNonConstGen = 0;
   if ( aGen )
   {
        aNonConstGen = const_cast<OlapGenerator*>(aGen);
        aNonConstGen->SetGrid(x,y,z);
   }
   else
   {
      G4cout << "Warning: Primary generator is not OlapGenerator!" << G4endl
             << "         Overlap Detection will not work!" << G4endl;
      return;
   }     
   
   eventsPerRun = x*y + y*z + x*z; 
   G4cout << "/olap/trigger will trigger "
          << eventsPerRun 
          << " events." << G4endl;
}


void OlapManager::notifyNewWorld(G4LogicalVolume* nw)
{
   std::set<OlapNotify*>::iterator i = theNotifs.begin();
   for(;i!=theNotifs.end();++i)
     (*i)->worldChanged(nw);   
}


void OlapManager::notifyOlaps(const std::vector<OlapInfo*> & ov)
{
   std::set<OlapNotify*>::iterator i = theNotifs.begin();
   for(;i!=theNotifs.end();++i) 
     (*i)->overlaps(ov);   
}


G4LogicalVolume* OlapManager::GetOriginalWorld()
{
   if (theDet)
   {
     return theDet->GetOriginalWorld();
   }
   return 0;
}
