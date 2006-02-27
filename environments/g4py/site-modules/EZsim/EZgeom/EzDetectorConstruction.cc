// $Id: EzDetectorConstruction.cc,v 1.1 2006-02-27 09:46:31 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   EzDetectorConstruction.cc
//
//                                         2005 Q
// ====================================================================
#include "EzDetectorConstruction.hh"
#include "G4EzWorld.hh"

// ====================================================================
//
// class description
//
// ====================================================================

////////////////////////////////////////////////
EzDetectorConstruction::EzDetectorConstruction()
////////////////////////////////////////////////
{
}

/////////////////////////////////////////////////
EzDetectorConstruction::~EzDetectorConstruction()
/////////////////////////////////////////////////
{
}

//////////////////////////////////////////////////////
G4VPhysicalVolume* EzDetectorConstruction::Construct()
//////////////////////////////////////////////////////
{
  return G4EzWorld::GetWorldVolume();
}

