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
// $Id: pyG4RayTracer.cc,v 1.4 2006-06-07 05:03:10 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4RayTracer.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4RayTracer.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4RayTracer()
{
  class_<G4RayTracer, G4RayTracer*, bases<G4VGraphicsSystem> >
    ("G4RayTracer", "RayTracer visualization module")
    ;
}

