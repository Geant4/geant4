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
// $Id: StateProcessMap.cc,v 1.2 2002-06-03 12:09:33 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#include "StateProcessMap.hh"
#include "SAXComponentFactoryTable.hh"
#include "SAXComponentFactory.hh"

#include <iostream>

StateProcessMap::StateProcessMap()
{
  std::cout << "State process map created" << std::endl;
}

StateProcessMap::~StateProcessMap()
{
  Reset();
}

void StateProcessMap::Initialize()
{
  SAXComponentFactoryTable* FTable = SAXComponentFactoryTable::GetInstance();
  const SAXComponentFactoryTable::Components* factories;
  
  factories = FTable->GetComponents( SAXComponentObject::eProcess );
  
  SAXComponentFactoryTable::Components::const_iterator pit;
  
  StateProcessMap::Process process;
  for( pit = factories->begin(); pit != factories->end(); pit++ )
  {
    process = dynamic_cast<SAXStateProcess*>( (*pit)->Build() );
    
    std::cout << "Adding state process for tag: " << process->State() << std::endl;
    
    AddProcess( process->State(), process );
  }
  std::cout << "========================================================" << std::endl;
}

void StateProcessMap::Reset()
{
  fMap.clear();
}

void StateProcessMap::AddProcess( const std::string& tag, StateProcessMap::Process process )
{
  fMap[tag] = process;
}

void StateProcessMap::AddProcess( char* tag, StateProcessMap::Process process )
{
  std::string key = tag;
  AddProcess( key, process );
}

bool StateProcessMap::Check( const std::string& tag )
{
  return( ( fMap.find( tag ) != fMap.end() ) ? true : false );
}

bool StateProcessMap::Check( char*        tag )
{
  std::string key = tag;
  return( Check( key ) );
}

StateProcessMap::Process StateProcessMap::GetProcess( const std::string& tag )
{
  return( fMap[tag] );
}

StateProcessMap::Process StateProcessMap::GetProcess( char*        tag )
{
  std::string key = tag;
  return( GetProcess( key ) );
}


