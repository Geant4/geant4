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
// $Id: OlapLogManager.cc,v 1.1 2002-06-04 07:40:21 gcosmo Exp $
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
#include "g4std/fstream"

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
