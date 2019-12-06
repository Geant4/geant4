#ifndef tools_Xt_OpenGLArea_h
#define tools_Xt_OpenGLArea_h

#include <X11/Intrinsic.h>

typedef struct _OpenGLAreaClassRec* OpenGLAreaWidgetClass;
typedef struct _OpenGLAreaRec* OpenGLAreaWidget;

typedef struct {
 int reason;
 XEvent* event;
} XoAnyCallbackStruct;

/* OpenGLArea */
#define XoNdoubleBufferOn "doubleBufferOn"
#define XoNpaintCallback  "paintCallback"
#define XoNeventCallback  "eventCallback"

#define XoCR_PAINT 1
#define XoCR_EVENT 2

#ifdef __cplusplus
extern "C"{
#endif
extern WidgetClass openGLAreaWidgetClass;

void OpenGLAreaPaint(Widget);
int OpenGLAreaWrite_gl2ps(Widget,const char*,const char*);

#ifdef __cplusplus
}
#endif

#endif
/*
exlib_build_use Xt
*/
