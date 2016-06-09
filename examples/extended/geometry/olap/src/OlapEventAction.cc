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
/// \file geometry/olap/src/OlapEventAction.cc
/// \brief Implementation of the OlapEventAction class
//
//
// $Id$
//
// 
// --------------------------------------------------------------
// OlapEventAction
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#include "globals.hh"
#include <vector>
#include <sstream>

#include "G4VVisManager.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4UImanager.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"

#include "OlapManager.hh"
#include "OlapEventAction.hh"
#include "OlapRunAction.hh"
#include "OlapDetConstr.hh"
#include "OlapGenerator.hh"

//#define OLAP_DEBUG

std::ostream & 
   operator<<(std::ostream& flux, OlapInfo & oi)
{
    OlapStepInfo si;
    si.thePoint = oi.v1;
    si.theHist = oi.hist1;
    flux << "vol 1: point=" << oi.v1 << " vol 2: point=" << oi.v2
         << " axis=" << oi.axis << "  [" << oi.probNot << "] "
         << oi.info << G4endl; 
    if (oi.stAB.size())
    {
      flux << oi.stAB << G4endl << G4endl;
      flux << oi.stBA << G4endl;
    }           
    flux << "A -> B:" << G4endl;
    flux << si << G4endl;
    flux << "B -> A:" << G4endl;
    si.thePoint=oi.v2;
    si.theHist=oi.hist2;
    flux << si << G4endl;
    return flux;  
}   


OlapInfo::~OlapInfo()
{
   while( !stAB.empty() )
   {
      delete stAB.back();
      stAB.pop_back();
   }   

   while( !stBA.empty() )
   {
      delete stBA.back();
      stBA.pop_back();
   }   

}

G4bool OlapInfo::operator==(const OlapInfo & rh) 
{
   G4bool result = false;   
   if (
        (
          ( hist1.GetTopVolume()    == rh.hist1.GetTopVolume() ) &&
          ( hist2.GetTopVolume()    == rh.hist2.GetTopVolume() ) 
        ) 
          ||
        (
          ( hist1.GetTopVolume()    == rh.hist2.GetTopVolume() ) &&
          ( hist2.GetTopVolume()    == rh.hist1.GetTopVolume() )         
        ) 
      )         
       result = true;
   return result;        
}


OlapEventAction::OlapEventAction(OlapRunAction * aRunAction)
  : theRunAction(aRunAction), dontDelete(false)
{ 
    const G4VUserDetectorConstruction * aDet =
      G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    const OlapDetConstr * aConstDet =
      dynamic_cast<const OlapDetConstr *>(aDet);
    if (!aDet)
    {
      G4Exception("OlapEventAction::OlapEventAction()",
                  "InvalidSetup", FatalException,
                  "Can't be used together with an OlapDetConstr! Exiting...");
    }
    theDet = const_cast<OlapDetConstr *>(aConstDet);                     
}


OlapEventAction::~OlapEventAction()
{
}


void OlapEventAction::BeginOfEventAction(const G4Event*)
{
    while (!ABSteps.empty())
    {
       if (! dontDelete ) 
         delete ABSteps.back(); 
       ABSteps.pop_back();
    }  
    
    while (!BASteps.empty())
    {
       if (! dontDelete )
          delete BASteps.back();
       BASteps.pop_back();
    }   
    
    dontDelete=false;
}   

   
void OlapEventAction::EndOfEventAction(const G4Event* anEvent)   
{
   
   const OlapGenerator * aGen = dynamic_cast<const OlapGenerator *>
         (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
   OlapGenerator * aGenerator=0;
   if ( aGen )
   {
        aGenerator = const_cast<OlapGenerator*>(aGen);
   }
   else
   {
      G4cout << "Warning: Primary generator is not OlapGenerator!" << G4endl
             << "         Overlap Detection will not work!" << G4endl;
      return;             
   }             
   
   if (! (anEvent->GetEventID() % 1000))
   {
     G4cout << "Event=" << anEvent->GetEventID() << " : "
        << theDet->GetNewWorld()->GetLogicalVolume()->GetDaughter(0)->GetName()
        << G4endl;
   }            
   
   G4int axx = -1;
   if (aGenerator)
     axx = aGenerator->GetAxis();

   G4double tolerance = OlapManager::GetOlapManager()->Delta();
   
   if ( (ABSteps[ABSteps.size()-1])->theHist.GetTopVolume()->GetLogicalVolume() != 
          theDet->GetNewWorld()->GetLogicalVolume()  ||
        (BASteps[BASteps.size()-1])->theHist.GetTopVolume()->GetLogicalVolume() != 
          theDet->GetNewWorld()->GetLogicalVolume()  
      )
    { 
       OlapInfo * oi =
          new OlapInfo( ABSteps[ABSteps.size()-1]->theHist,
                        BASteps[BASteps.size()-1]->theHist,
                        ABSteps[ABSteps.size()-1]->thePoint,
                        BASteps[BASteps.size()-1]->thePoint,
                        axx,
                        OlapManager::GetOlapManager()->GetOriginalWorld() );
       oi->info = "daughter protrudes mother";                                       
       oi->probNot = false; // it's an overlap
       oi->stAB = ABSteps;
       oi->stBA = BASteps;
       dontDelete = true;
                                       
       if ( !InsertOI(oi) )
       {                   // already detected!
            delete oi;
       }
       return;
    }  
   
   // just in case, not the OlapGenerator is used but
   // another particle generator
   if (!BASteps.size())
     return;
    
   // do we have as many steps from AB as from BA?
   if (ABSteps.size() != BASteps.size())
   {
     OlapInfo * oi =
        new OlapInfo( ABSteps[1]->theHist,
                      BASteps[1]->theHist,
                      ABSteps[1]->thePoint,
                      BASteps[1]->thePoint,
                      axx,
                      OlapManager::GetOlapManager()->GetOriginalWorld() );
     // copy all steps and prohibit deletion of their pointers ...
     oi->stAB = ABSteps;
     oi->stBA = BASteps;
     oi->probNot = true; // probably not an overlap, just numerics ....
     dontDelete = true;
          
     std::ostringstream os;
     os << "A,B diff. steps AB:" << ABSteps.size()
                      << " BA:"  << BASteps.size() << '\0';

     oi->info = G4String(os.str());
     if ( !InsertOI(oi) )
     {
        delete oi;
     }                               
     return;
   }
 
   // now forget about the first and the last points of AB and BA &     
   // compare the remaining points of AB with the remaining points
   // of BA where the points of BA are in descending order while
   // the points of BA are in ascending order
   G4int siz = ABSteps.size();
   G4int j = siz-3;
   for (G4int i = 1; i<(siz-1); i++)
   {
       G4double delta =
          (ABSteps[i]->thePoint - BASteps[j]->thePoint).mag() ;
       if (delta>tolerance) 
       {
          OlapInfo * oi =
              new OlapInfo( ABSteps[i]->theHist,
                            BASteps[j]->theHist,
                            ABSteps[i]->thePoint,
                            BASteps[j]->thePoint,
                            axx,
                            OlapManager::GetOlapManager()->GetOriginalWorld() );
          if ( !InsertOI(oi) )
          {                       // already detected!
            delete oi;
          }
       }  
       j--;         
   }
}  

G4bool OlapEventAction::InsertOI(OlapInfo * oi)
{
   std::vector<OlapInfo*>::iterator it = theRunAction->theOlaps.begin();
   G4bool flag(false);
   while (it != theRunAction->theOlaps.end())
   {
     if (*oi== *(*it)) { // already detected overlap
       flag=true;
       if ((*it)->probNot != oi->probNot )
       {                 // probably a real overlap, not numerics
         flag=false; 
         (*it)->probNot = false;
         oi->probNot = false;
       } 
       break;    
     }                
     it++;
   }
   
   if (flag)
   {
     return false; // already detected overlap
   }  
     
   theRunAction->theOlaps.push_back(oi);
   
   #ifdef OLAP_DEBUG
     G4cout << "===================================" << G4endl;
     G4cout << "No of detected Overlaps: " << theRunAction->theOlaps.size()
            << G4endl; 
     G4cout << "===================================" << G4endl;
     G4cout << oi->hist1 << G4endl;
     G4cout << "CpNo: " << oi->hist1.GetTopVolume()->GetCopyNo() << G4endl;
     G4cout << "Pos:  " << oi->hist1.GetTopTransform().NetTranslation()
            << G4endl;
     G4cout << "Rot:  " << oi->hist1.GetTopTransform().NetRotation().phiX()
            << " " << oi->hist1.GetTopTransform().NetRotation().phiY() << " "
                   << oi->hist1.GetTopTransform().NetRotation().phiZ() << " "
                   << oi->hist1.GetTopTransform().NetRotation().thetaX() << " "
                   << oi->hist1.GetTopTransform().NetRotation().thetaY() << " "
                   << oi->hist1.GetTopTransform().NetRotation().thetaZ()
                   << G4endl ;
     G4cout << oi->v1 << G4endl << "----" << G4endl;                        

     G4cout << oi->hist2 << G4endl;
     G4cout << "CpNo: " << oi->hist2.GetTopVolume()->GetCopyNo() << G4endl;
     G4cout << "Pos:  " << oi->hist2.GetTopTransform().NetTranslation()
            << G4endl;
     G4cout << "Rot:  " << oi->hist2.GetTopTransform().NetRotation().phiX()
            << " " << oi->hist2.GetTopTransform().NetRotation().phiY() << " "
                   << oi->hist2.GetTopTransform().NetRotation().phiZ() << " "
                   << oi->hist2.GetTopTransform().NetRotation().thetaX() << " "
                   << oi->hist2.GetTopTransform().NetRotation().thetaY() << " "
                   << oi->hist2.GetTopTransform().NetRotation().thetaZ() 
                   << G4endl ;
     G4cout << oi->v2 << G4endl << "----" << G4endl;                       

   #endif
  
   static char * c = getenv("OLAP_LOG_EVENT");
   if (c)
     G4cerr << *oi << G4endl;
    
   return true; // new overlap detected!  
}
