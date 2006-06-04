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
// $Id: BClass.cc,v 1.3 2006-06-04 21:35:59 kmura Exp $
// ====================================================================
//   BClass.cc
//
//                                         2004 Q
// ====================================================================
#include "BClass.hh"
#include <iostream>

// ====================================================================
//
// class description
//
// ====================================================================

////////////////
BClass::BClass()
  : XBase()
////////////////
{
  ival= -3;
  dval= -3.;
}

/////////////////
BClass::~BClass()
/////////////////
{
}

//////////////////////
void BClass::AMethod()
//////////////////////
{
  std::cout << "%%% BClass:::AMethod is called."
	    << " (ival, dval)= (" << ival << "," << dval << ")"
	    << std::endl;
}

/////////////////////////////////////////////
int BClass::VMethod(const XBase* abase) const
/////////////////////////////////////////////
{
  return abase-> GetIVal();
}

