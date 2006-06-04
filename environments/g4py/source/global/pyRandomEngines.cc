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
// $Id: pyRandomEngines.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyRandomEngines.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/Ranlux64Engine.h"

using namespace boost::python;
using namespace CLHEP;

// ====================================================================
// module definition
// ====================================================================
void export_RandomEngines()
{
  class_<HepRandomEngine, boost::noncopyable>
    ("HepRandomEngine", "base class of random engine", no_init)
    ;

  class_<HepJamesRandom, bases<HepRandomEngine> >
    ("HepJamesRandom", "HepJames random engine")
    ;

  class_<Ranlux64Engine, bases<HepRandomEngine> >
    ("Ranlux64Engine", "Ranlux64 random engine")
    ;
}
