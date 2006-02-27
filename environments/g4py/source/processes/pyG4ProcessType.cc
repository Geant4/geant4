// $Id: pyG4ProcessType.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4ProcessTyoe.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4ProcessType.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4ProcessType()
{
 enum_<G4ProcessType>("G4ProcessType")
   .value("fNotDefined",         fNotDefined)
   .value("fTransportation",     fTransportation)
   .value("fElectromagnetic",    fElectromagnetic)
   .value("fOptical",            fOptical)
   .value("fHadronic",           fHadronic)
   .value("fPhotolepton_hadron", fPhotolepton_hadron)
   .value("fDecay",              fDecay)
   .value("fGeneral",            fGeneral)
   .value("fParameterisation",   fParameterisation)
   .value("fUserDefined",        fUserDefined)
   ;
}

