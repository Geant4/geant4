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
// $Id: AClass.cc,v 1.3 2006-06-04 21:35:59 kmura Exp $
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

