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

