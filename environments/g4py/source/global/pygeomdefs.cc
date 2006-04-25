// $Id: pygeomdefs.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pygeomdefs.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "geomdefs.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_geomdefs()
{
 enum_<EAxis>("EAxis")
   .value("kXAxis",     kXAxis)
   .value("kYAxis",     kYAxis)
   .value("kZAxis",     kZAxis)
   .value("kRho",       kRho)
   .value("kRadial3D",  kRadial3D)
   .value("kPhi",       kPhi)
   .value("kUndefined", kUndefined)
   ;

 enum_<EInside>("EInside")
   .value("kOutside",   kOutside)
   .value("kSurface",   kSurface)
   .value("kInside",    kInside)
   ;

 enum_<EVolume>("EVolume")
   .value("kNormal",        kNormal)
   .value("kReplica",       kReplica)
   .value("kParameterised", kParameterised)
   ;
}

