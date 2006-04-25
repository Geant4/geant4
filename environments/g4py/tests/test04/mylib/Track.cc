// $Id: Track.cc,v 1.2 2006-04-25 07:54:28 kmura Exp $
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

