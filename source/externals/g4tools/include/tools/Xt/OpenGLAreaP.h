#ifndef tools_Xt_OpenGLAreaP_h
#define tools_Xt_OpenGLAreaP_h

#include "OpenGLArea.h"

#include <X11/IntrinsicP.h>
#include <X11/CoreP.h>
#include <X11/CompositeP.h>

#include <GL/glx.h>

typedef struct
{
  void* extension;
} OpenGLAreaClassPart;

typedef struct _OpenGLAreaClassRec
{
  CoreClassPart core_class;
  CompositeClassPart composite_class;
  OpenGLAreaClassPart openGLArea_class;
} OpenGLAreaClassRec;

#ifdef __cplusplus
extern "C"{
#endif
extern OpenGLAreaClassRec openGLAreaClassRec;
#ifdef __cplusplus
}
#endif

typedef struct
{
  /*Resources :*/
  Boolean doubleBufferOn;
  XtCallbackList paintCallback;
  XtCallbackList eventCallback;
  /**/
  Visual* visual;
  Boolean installColormap;
  GLXContext glContext;
} OpenGLAreaPart;

typedef struct _OpenGLAreaRec
{
  CorePart core;
  CompositePart composite;
  OpenGLAreaPart openGLArea;
} OpenGLAreaRec;

#endif
/*
exlib_build_use Xt
*/
