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
//
// $Id: G4OpenGLXmBox.hh,v 1.6 2001-07-11 10:08:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//Box container class

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMBOX_HH
#define G4OPENGLXMBOX_HH

#include "G4OpenGLXmVWidgetContainer.hh"
#include "globals.hh"
#include <Xm/Frame.h>
#include <Xm/RowColumn.h>

class G4OpenGLXmVWidgetComponent;
class G4OpenGLXmVWidgetShell;

class G4OpenGLXmBox : public G4OpenGLXmVWidgetContainer
{

public:
  G4OpenGLXmBox (const char* = NULL,
		 G4bool = False);   //constructor
  virtual ~G4OpenGLXmBox ();  //destructor

  void AddChild (G4OpenGLXmVWidgetComponent*);
  void AddYourselfTo (G4OpenGLXmVWidgetShell*);

  Widget* GetPointerToParent ();
  Widget* GetPointerToWidget ();
  
  const char* GetName ();
  void SetName (const char*);

protected:
  const char* name;
  Widget* parent;
  Widget box_row_col;
  G4bool radio;

private:
  G4OpenGLXmBox (const G4OpenGLXmBox&);
  G4OpenGLXmBox& operator = (const G4OpenGLXmBox&);
};

#endif

#endif
