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
// $Id: G4OpenGL2PSAction.hh 88190 2015-02-02 17:24:54Z gcosmo $
//

#ifndef G4OpenGL2PSAction_h
#define G4OpenGL2PSAction_h 


#include <string>
#include "Geant4_gl2ps.h"

class G4OpenGL2PSAction {

public:
 G4OpenGL2PSAction();
 
 void setFileName(const char*);
 void setExportImageFormat(unsigned int);
 bool enableFileWriting();
  // return true if ok, false is an error occured

 bool disableFileWriting();
  // return true when OK, false if errror

 bool fileWritingEnabled() const;
 void setLineWidth(int);
 void setPointSize(int);
 void setViewport(int,int,int,int);
 bool extendBufferSize();
 void resetBufferSizeParameters();
 void setBufferSize(int);

protected:
 bool G4gl2psBegin();
 const char * fFileName;
 FILE* fFile;
 GLint fViewport[4];
 GLint fBufferSize;
 GLint fBufferSizeLimit;
private:
 unsigned int fExportImageFormat;
};

#endif

