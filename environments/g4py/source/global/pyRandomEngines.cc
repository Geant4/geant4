// $Id: pyRandomEngines.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
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
