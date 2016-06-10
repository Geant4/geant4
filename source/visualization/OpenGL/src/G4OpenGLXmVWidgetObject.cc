//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4OpenGLXmVWidgetObject.cc 68043 2013-03-13 14:27:49Z gcosmo $
//
//Virtual base class for all Motif widgets.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmViewer.hh"
#include "G4OpenGLXmVWidgetObject.hh"

G4OpenGLXmVWidgetObject::G4OpenGLXmVWidgetObject ()
: pView(0)
, visual(0)
, top(0)
{}

G4OpenGLXmVWidgetObject::~G4OpenGLXmVWidgetObject ()
{}

G4OpenGLXmViewer* G4OpenGLXmVWidgetObject::GetView ()
{
  return pView;
}

void G4OpenGLXmVWidgetObject::ProcesspView () 
{
  cmap = pView->cmap;
  depth = pView->vi->depth;
  visual = pView->vi->visual;
  borcol = pView->borcol;
  bgnd = pView->bgnd;
  top = pView->toplevel;
}

#endif
