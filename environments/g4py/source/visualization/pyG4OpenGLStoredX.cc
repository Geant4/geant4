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
// $Id: pyG4OpenGLStoredX.cc,v 1.4 2006-06-07 05:03:10 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4OpenGLStoredX.cc
//
//                                         2005 Q
// ====================================================================
#ifdef G4VIS_USE_OPENGLX

#include <boost/python.hpp>
#include "G4OpenGLStoredX.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4OpenGLStoredX()
{
  class_<G4OpenGLStoredX, G4OpenGLStoredX*, bases<G4VGraphicsSystem> >
    ("G4OpenGLStoredX", "OpenGL(Stored X) visualization module")
    ;
}

#endif
