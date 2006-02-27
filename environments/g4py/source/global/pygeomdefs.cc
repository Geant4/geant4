// $Id: pygeomdefs.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
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

