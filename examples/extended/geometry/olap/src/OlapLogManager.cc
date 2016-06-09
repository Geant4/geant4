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
// $Id: OlapLogManager.cc,v 1.3 2006-06-29 17:22:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// OlapLogManager
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#include "OlapLogManager.hh"

#include "globals.hh"
#include <fstream>

OlapLogManager * OlapLogManager::theInstance = 0;


OlapLogManager::OlapLogManager()
{
  areWeLogging = false;
  areWeLoggingByVolume = false;
}

OlapLogManager::~OlapLogManager() {}

void OlapLogManager::Logging(G4String aValue)
{
  if(aValue == "false")
  {
    areWeLogging = false; 
  }
  else
  {
    G4int lastLetter = aValue.length();
    G4String::size_type c = lastLetter-1;
    if (aValue[c] == '/')
      aValue = aValue + "olap.log";

    areWeLogging = true;
    filename = aValue;
    G4cerr << "Output will be put into " << aValue << G4endl;
  }
}


void OlapLogManager::LogByVolume(G4String aPath)
{
  if(aPath == "false")
  {
    areWeLoggingByVolume = false;
  }
  else
  {
    areWeLoggingByVolume = true;
    logPath = aPath;
    G4cerr << "Output will be put into " << aPath
           << "NameOfVolume.log" << G4endl;
  }
}


OlapLogManager * OlapLogManager::GetOlapLogManager()
{
   if (!theInstance)
      theInstance = new OlapLogManager();
   
   return theInstance;   
}   
