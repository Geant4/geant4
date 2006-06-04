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
// $Id: pymodG4visualization.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pymodG4visualization.cc [Geant4Py module]
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================

void export_G4VisManager();
void export_G4VGraphicsSystem();
void export_G4VRML1File();
void export_G4VRML2File();
void export_G4DAWNFILE();
void export_G4HepRep();
void export_G4HepRepFile();
void export_G4ASCIITree();
void export_G4RayTracer();

#ifdef G4VIS_USE_OPENGLX
void export_G4OpenGLStoredX();
void export_G4OpenGLImmediateX();
#endif

#ifdef G4VIS_USE_OPENGLXM
void export_G4OpenGLStoredXm();
void export_G4OpenGLImmediateXm();
#endif

#ifdef G4VIS_USE_RAYTRACERX
void export_G4RayTracerX();
#endif

BOOST_PYTHON_MODULE(G4visualization) 
{
  export_G4VisManager();
  export_G4VGraphicsSystem();
  export_G4VRML1File();
  export_G4VRML2File();
  export_G4DAWNFILE();
  export_G4HepRep();
  export_G4HepRepFile();
  export_G4ASCIITree();
  export_G4RayTracer();

#ifdef G4VIS_USE_OPENGLX
  export_G4OpenGLStoredX();
  export_G4OpenGLImmediateX();
#endif

#ifdef G4VIS_USE_OPENGLXM
  export_G4OpenGLStoredXm();
  export_G4OpenGLImmediateXm();
#endif

#ifdef G4VIS_USE_RAYTRACERX
  export_G4RayTracerX();
#endif

}

