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
// $Id: GDMLProcessor.cc,v 1.3 2002-08-19 07:35:50 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#include "GDMLProcessor.hh"

// Declare the external component system initialization routines
extern "C" void GDMLProcessLibLoad();
extern "C" void GDMLSubscriberLibLoad();

const double defaultTemp  = STP_Temperature;
const double defaultPress = STP_Pressure;

static GDMLProcessor* sProcessor = 0;

GDMLProcessor* GDMLProcessor::GetInstance()
{
  if( sProcessor == 0 )
  {
    sProcessor = new GDMLProcessor();
  }
  
  return sProcessor;
}

GDMLExpressionEvaluator* GDMLProcessor::GetEvaluator()
{
  return fCalc;
}
  
GDMLProcessor::~GDMLProcessor()
{
  if( fCalc != 0 )
  {
    delete fCalc;
    fCalc = 0;
  }
  
  GDMLProcessor::Positions::iterator       pit;
  GDMLProcessor::Rotations::iterator       rit;
  GDMLProcessor::Solids::iterator          sit;
  GDMLProcessor::LogicalVolumes::iterator  lvit;
  GDMLProcessor::AssemblyVolumes::iterator avit;
  GDMLProcessor::PhysicalVolumes::iterator pvit;
  
  for( pit = fPTable.begin(); pit != fPTable.end(); pit++ )
  {
    G4ThreeVector* victim = (*pit).second;
    if( victim != 0 )
    {
      delete victim;
    }
  }
  for( rit = fRTable.begin(); rit != fRTable.end(); rit++ )
  {
    G4RotationMatrix* victim = (*rit).second;
    if( victim != 0 )
    {
      delete victim;
    }
  }
  for( sit = fSolids.begin(); sit != fSolids.end(); sit++ )
  {
    G4VSolid* victim = (*sit).second;
    if( victim != 0 )
    {
      delete victim;
    }
  }
  for( lvit = fLVolumes.begin(); lvit != fLVolumes.end(); lvit++ )
  {
    G4LogicalVolume* victim = (*lvit).second;
    if( victim != 0 )
    {
      delete victim;
    }
  }
  for( avit = fAVolumes.begin(); avit != fAVolumes.end(); avit++ )
  {
    G4AssemblyVolume* victim = (*avit).second;
    if( victim != 0 )
    {
      delete victim;
    }
  }
  for( pvit = fPVolumes.begin(); pvit != fPVolumes.end(); pvit++ )
  {
    G4VPhysicalVolume* victim = (*pvit).second;
    if( victim != 0 )
    {
      delete victim;
    }
  }
}

void GDMLProcessor::AddPosition( const std::string& name, G4ThreeVector* p )
{
  fPTable[name] = p;
}

void GDMLProcessor::AddPosition( const char* name, G4ThreeVector* p )
{
  std::string key = name;
  AddPosition( key, p );
}

void GDMLProcessor::AddRotation( const std::string& name, G4RotationMatrix* p )
{
  fRTable[name] = p;
}

void GDMLProcessor::AddRotation( const char* name, G4RotationMatrix* p )
{
  std::string key = name;
  AddRotation( key, p );
}

void GDMLProcessor::AddSolid( const std::string& name, G4VSolid* p )
{
  fSolids[name] = p;
}
void GDMLProcessor::AddSolid( const char* name, G4VSolid* p )
{
  std::string key = name;
  AddSolid( key, p );
}

void GDMLProcessor::AddLogicalVolume( const std::string& name, G4LogicalVolume* p )
{
  fLVolumes[name] = p;
}
void GDMLProcessor::AddLogicalVolume( const char* name, G4LogicalVolume* p )
{
  std::string key = name;
  AddLogicalVolume( key, p );
}

void GDMLProcessor::AddAssemblyVolume( const std::string& name, G4AssemblyVolume* p )
{
  fAVolumes[name] = p;
}
void GDMLProcessor::AddAssemblyVolume( const char* name, G4AssemblyVolume* p )
{
  std::string key = name;
  AddAssemblyVolume( key, p );
}

void GDMLProcessor::AddPhysicalVolume( const std::string& name, G4VPhysicalVolume* p )
{
  fPVolumes[name] = p;
}
void GDMLProcessor::AddPhysicalVolume( const char* name, G4VPhysicalVolume* p )
{
  std::string key = name;
  AddPhysicalVolume( key, p );
}

const G4ThreeVector*    GDMLProcessor::GetPosition( const std::string& name )
{
  return fPTable[name];
}

const G4ThreeVector*    GDMLProcessor::GetPosition( const char* name )
{
  std::string key = name;
  return GetPosition( key );
}

const G4RotationMatrix* GDMLProcessor::GetRotation( const std::string& name )
{
  return fRTable[name];
}

const G4RotationMatrix* GDMLProcessor::GetRotation( const char* name )
{
  std::string key = name;
  return GetRotation( key );
}

const G4VSolid* GDMLProcessor::GetSolid( const std::string& name )
{
  return fSolids[name];
}

const G4VSolid* GDMLProcessor::GetSolid( const char* name )
{
  std::string key = name;
  return GetSolid( key );
}

const G4LogicalVolume* GDMLProcessor::GetLogicalVolume( const std::string& name )
{
  return fLVolumes[name];
}

const G4LogicalVolume* GDMLProcessor::GetLogicalVolume( const char* name )
{
  std::string key = name;
  return GetLogicalVolume( key );
}

const G4AssemblyVolume* GDMLProcessor::GetAssemblyVolume( const std::string& name )
{
  return fAVolumes[name];
}

const G4AssemblyVolume* GDMLProcessor::GetAssemblyVolume( const char* name )
{
  std::string key = name;
  return GetAssemblyVolume( key );
}

const G4VPhysicalVolume* GDMLProcessor::GetPhysicalVolume( const std::string& name )
{
  return fPVolumes[name];
}

const G4VPhysicalVolume* GDMLProcessor::GetPhysicalVolume( const char* name )
{
  std::string key = name;
  return GetPhysicalVolume( key );
}

GDMLProcessor::GDMLProcessor()
: fCalc( 0 ), fWorldVolume( 0 )
{
  fCalc = new GDMLExpressionEvaluator();
  
  // We need to initialize our component system
  GDMLProcessLibLoad();
  GDMLSubscriberLibLoad();
}


