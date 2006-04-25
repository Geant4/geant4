// $Id: pyG4ProcessType.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
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

