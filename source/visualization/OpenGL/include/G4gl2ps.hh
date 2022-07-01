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

#ifndef G4gl2ps_h
#define G4gl2ps_h 

#include "G4String.hh"

#include <tools/gl2ps_def.h>

class G4gl2ps {

public:
 G4gl2ps();
 ~G4gl2ps();

 void setOpenGLFunctions(tools_gl2ps_gl_funcs_t*);
 void setFileName(const char*);
 void setExportImageFormat(unsigned int);
  
 void setExportImageFormat_PS()  {setExportImageFormat(TOOLS_GL2PS_PS);}
 void setExportImageFormat_EPS() {setExportImageFormat(TOOLS_GL2PS_EPS);}
 void setExportImageFormat_TEX() {setExportImageFormat(TOOLS_GL2PS_TEX);}
 void setExportImageFormat_PDF() {setExportImageFormat(TOOLS_GL2PS_PDF);}
 void setExportImageFormat_SVG() {setExportImageFormat(TOOLS_GL2PS_SVG);}
 void setExportImageFormat_PGF() {setExportImageFormat(TOOLS_GL2PS_PGF);}
  
 bool enableFileWriting();
 void disableFileWriting();
 bool fileWritingEnabled() const;
  
 bool beginPage();
 bool endPage();
  
 void setLineWidth(int);
 void setPointSize(int);
 void addTextOpt(const char*,const char*,
                 tools_GLshort,tools_GLint,tools_GLfloat);
 void setViewport(int,int,int,int);
 bool extendBufferSize();
 void resetBufferSizeParameters();
 void setBufferSize(int);

 tools_GL2PScontextPointer context() const {return fContext;}
protected:
 tools_gl2ps_gl_funcs_t fOpenGLFuncs;
 tools_GL2PScontextPointer fContext;
 FILE* fFile;
 G4String fFileName;
 int fViewport[4];
 int fBufferSize;
 int fBufferSizeLimit;
private:
 unsigned int fExportImageFormat;
};

#endif

