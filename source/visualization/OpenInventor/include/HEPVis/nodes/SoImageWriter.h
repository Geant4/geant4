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
////////////////////////////////////////////////////////////
//
//  This node is able to produce PostScript, GIF files. 
//
////////////////////////////////////////////////////////////


#ifndef HEPVis_SoImageWriter_h
#define HEPVis_SoImageWriter_h 

#include <Inventor/nodes/SoSubNode.h>
#include <Inventor/fields/SoSFEnum.h>
#include <Inventor/fields/SoSFString.h>

#define SoImageWriter Geant4_SoImageWriter

class SoImageWriter : public SoNode {
  SO_NODE_HEADER(SoImageWriter);
public:
/*  enum Format {
    POST_SCRIPT,
    GIF
  };
  SoSFEnum format;*/
  SoSFString fileName;
  SoImageWriter();
  void enable();
  void disable();
  SbBool getStatus() const;
protected:
  virtual  ~SoImageWriter();
public: /*SoEXTENDER*/
  virtual void GLRender(SoGLRenderAction*);
public: /*SoINTERNAL*/
  static void initClass();
private:
  SbBool fEnabled;
  SbBool fStatus;
};

#endif

