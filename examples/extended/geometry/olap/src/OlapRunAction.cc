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
// $Id: OlapRunAction.cc,v 1.2 2002-06-04 09:23:45 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// OlapRunAction
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#include "OlapRunAction.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"

#include "OlapManager.hh"

#include "G4Run.hh"

#include "globals.hh"
#include "g4std/fstream"

OlapRunAction::OlapRunAction() 
{ 
}


OlapRunAction::~OlapRunAction()
{
}

            
void OlapRunAction::BeginOfRunAction(const G4Run* aRun)
{
  while ( !theOlaps.empty() )
  {
     delete theOlaps.back();
     theOlaps.pop_back();
  }   
}


void OlapRunAction::EndOfRunAction(const G4Run* aRun)
{
   if (!theOlaps.size())
     return;

//output to screen
   G4cerr << "===== collected overlaps of run [" << aRun->GetRunID() << "] " 
          << "(ol=" << theOlaps.size() << ")"<< G4endl;
   
   OlapManager::GetOlapManager()->notifyOlaps(theOlaps);
   
   G4std::vector<OlapInfo*>::iterator it = theOlaps.begin();
   G4int i=1;
   while (it!=theOlaps.end())
   {
      G4cerr << "--[" << i << "]--------------------------" << G4endl; 
      G4cerr << "delta=" << ((*it)->v1 - (*it)->v2).mag() << G4endl;
      G4cerr << *(*it) << G4endl;
   
      it++; i++;
   }

//output to file 

   OlapLogManager * logManager = OlapLogManager::GetOlapLogManager();
   if (logManager->areWeLogging || logManager->areWeLoggingByVolume)
   {     
     it = theOlaps.begin();
     i=1;
     G4String fname;

     if(logManager->areWeLogging) 
       fname  = logManager->filename;

     else if(logManager->areWeLoggingByVolume)
     {
       G4String volume;
       if((*it)->hist1.GetDepth() >= 1)
	 volume = (*it)->hist1.GetVolume(1)->GetName();
       else if((*it)->hist2.GetDepth() >= 1)
	 volume = (*it)->hist2.GetVolume(1)->GetName();
       else
	 G4cerr << "error: did not get the filename" << G4endl;

       fname = logManager->logPath + volume + ".log";
     }

     FILE.open(fname, G4std::ios::app);
     FILE << "===== collected overlaps of run [" << aRun->GetRunID() << "] " 
          << "(ol=" << theOlaps.size() << ")"<< G4endl;

     while (it!=theOlaps.end()) {

       FILE << "--[" << i << "]--------------------------" << G4endl; 
       FILE << "delta=" << ((*it)->v1 - (*it)->v2).mag() << G4endl;
       FILE << *(*it) << G4endl;

       it++; i++;
     }
     FILE << "-------------------------------------" << G4endl;
     FILE.close();

     G4cerr << "Output has also been put into "<< fname << G4endl;
   }//end if 

   G4cerr << "-------------------------------------" << G4endl;
   DrawOlaps();

}  

#include "g4std/vector"
#include "G4Circle.hh"
#include "G4Polyline.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "OlapManager.hh"

void OlapRunAction::DrawOlaps()
{
  // iguana takes over ...
}
