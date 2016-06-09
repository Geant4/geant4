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
// $Id: AClass.cc,v 1.4 2006-06-29 15:38:44 gunter Exp $
// ====================================================================
//   AClass.cc
//
//                                         2005 Q
// ====================================================================
#include "AClass.hh"
#include <iostream>

// ====================================================================
//
// class description
//
// ====================================================================

////////////////
AClass::AClass()
  : ival(0)
////////////////
{
}

///////////////////////////////
AClass::AClass(int i, double d)
  : ival(i)
///////////////////////////////
{
  std::cout << "*** AClass constructor, d= "<< d << std::endl;
}

/////////////////
AClass::~AClass()
/////////////////
{
}


/////////////////////
int AClass::AMethod()
/////////////////////
{
  return 0;
}

//////////////////////////
int AClass::AMethod(int i)
//////////////////////////
{
  return 1;
}


////////////////////////////////////
int AClass::AMethod(int i, double d)
////////////////////////////////////
{
  return 2;
}



//////////////////////
int  AClass::BMethod()
//////////////////////
{
  return 0;
}

////////////////////////////////
double AClass::BMethod(double d)
////////////////////////////////
{
  return 1.;

}

///////////////////////////////////////////////////
double AClass::CMethod(int i, double d1, double d2)
///////////////////////////////////////////////////
{
  return i*d1*d2;
}

