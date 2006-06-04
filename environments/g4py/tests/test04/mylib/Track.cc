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
// $Id: Track.cc,v 1.3 2006-06-04 21:35:59 kmura Exp $
// ====================================================================
//   Track.cc
//
//                                         2005 Q
// ====================================================================
#include "Track.hh"
#include "Step.hh"
#include <iostream>

// ====================================================================
//
// class description
//
// ====================================================================

//////////////
Track::Track()
//////////////
{
  std::cout << "Track is created. @" 
	    << this << std::endl;
  step= new Step();
}

///////////////
Track::~Track()
///////////////
{
  std::cout << "Track is deleted. @" 
	    << this << std::endl;
  delete step;
}


//////////////////////////////////
const Step* Track::GetStep() const
//////////////////////////////////
{
  return step;
}


///////////////////////////////////
const Step* Track::GetStep1() const
///////////////////////////////////
{
  return GetStep();
}

///////////////////////////////////
const Step* Track::GetStep2() const
///////////////////////////////////
{
  return GetStep();
}

///////////////////////////////////
const Step* Track::GetStep3() const
///////////////////////////////////
{
  return GetStep();
}

///////////////////////////////////
const Step* Track::GetStep4() const
///////////////////////////////////
{
  return GetStep();
}

