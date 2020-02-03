/*
#define DEBUG
*/

/*this*/
#include "OpenGLAreaP.h"

#include <X11/StringDefs.h>

#include <GL/glx.h>

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C"{
#endif
static void InitializeClass(void);
static void InitializeWidget(Widget, Widget,ArgList,Cardinal*);
static void RealizeWidget(Widget, XtValueMask*, XSetWindowAttributes*);
static void DestroyWidget(Widget);
static void ChangeWidgetSize(Widget);
static void DrawWidget(Widget, XEvent*, Region);
static Boolean SetValues(Widget,Widget,Widget,ArgList,Cardinal *);

static void EventHandler(Widget,XtPointer,XEvent*,Boolean*);
static int MakeCurrent(Widget);
static void XWidgetInstallColormap(Widget);
static void XWidgetUninstallColormap(Widget);
static Widget XWidgetGetShell(Widget);
#ifdef __cplusplus
}
#endif

#define athis ((OpenGLAreaWidget)This)->openGLArea
#define acur ((OpenGLAreaWidget)a_current)->openGLArea
#define CWarn printf
#define CWarnF printf

static XtResource resources [] = {
  {XoNdoubleBufferOn,"DoubleBufferOn",XtRBoolean,sizeof(Boolean),
   XtOffset(OpenGLAreaWidget,openGLArea.doubleBufferOn),XtRImmediate,(XtPointer)True},
  {XoNpaintCallback,XtCCallback,XtRCallback,sizeof(XtCallbackList),
   XtOffset(OpenGLAreaWidget,openGLArea.paintCallback),XtRImmediate,(XtPointer)NULL},
  {XoNeventCallback,XtCCallback,XtRCallback,sizeof(XtCallbackList),
   XtOffset(OpenGLAreaWidget,openGLArea.eventCallback),XtRImmediate,(XtPointer)NULL}
};

OpenGLAreaClassRec  openGLAreaClassRec = {
/* Core Class Part */
{
   (WidgetClass) &compositeClassRec, /* pointer to superclass ClassRec   */
    "OpenGLArea",                    /* widget resource class name       */
    sizeof(OpenGLAreaRec),           /* size in bytes of widget record   */
    InitializeClass,                 /* class_initialize                 */
    NULL,                            /* dynamic initialization           */
    FALSE,                           /* has class been initialized?      */
    InitializeWidget,                /* initialize                       */
    NULL,                            /* notify that initialize called    */
    RealizeWidget,                   /* XCreateWindow for widget         */
    NULL,                            /* widget semantics name to proc mapWidget*/
    0,                               /* number of entries in actions     */
    resources,                       /* resources for subclass fields    */
    XtNumber(resources),             /* number of entries in resources   */
    NULLQUARK,                       /* resource class quarkified        */
    TRUE,                            /* compress MotionNotify for widget */
    TRUE,                            /* compress Expose events for widget*/
    TRUE,                            /* compress enter and leave events  */
    TRUE,                            /* select for VisibilityNotify      */
    DestroyWidget,                   /* free data for subclass pointers  */
    ChangeWidgetSize,                /* geom manager changed widget size */
    DrawWidget,                      /* rediplay window                  */
    SetValues,                       /* set subclass resource values     */
    NULL,                            /* notify that SetValues called    */
    XtInheritSetValuesAlmost,        /* SetValues got "Almost" geo reply*/
    NULL,                            /* notify that get_values called    */
    XtInheritAcceptFocus,            /* assign input focus to widget     */
    XtVersion,                       /* version of intrinsics used       */
    NULL,                            /* list of callback offsets         */
    XtInheritTranslations,           /* translations                     */        
    XtInheritQueryGeometry,          /* return preferred geometry        */
    XtInheritDisplayAccelerator,     /* display your accelerator         */
    NULL                             /* pointer to extension record      */
},
/* Composite Class Part */
{
    XtInheritGeometryManager,   /* geometry manager for children   */
    XtInheritChangeManaged,     /* change managed state of child   */
    XtInheritInsertChild,       /* physically add child to parent  */
    XtInheritDeleteChild,       /* physically remove child         */
    NULL                        /* pointer to extension record     */
},
/* OpenGLArea */
{
   NULL 
}
   
};

WidgetClass openGLAreaWidgetClass = (WidgetClass) &openGLAreaClassRec;

static void InitializeClass(void) {}

static void InitializeWidget(Widget a_request,Widget This,ArgList a_args,Cardinal* a_argn) {
  if(a_request->core.width<=0)  This->core.width  = 100;
  if(a_request->core.height<=0) This->core.height = 100;

#ifdef DEBUG
  printf ("debug : OpenGLArea : InitializeWidget : %s\n",XtName(This));
#endif

  athis.visual          = CopyFromParent;
  athis.installColormap = False;
  athis.glContext       = NULL;
  
 {Display* display;
  Screen* screen;
  int iscreen;
  XVisualInfo* vinfo;
  display = XtDisplay(This);
  screen = XtScreen(This);
  iscreen = XScreenNumberOfScreen(screen);
  
 {int error,event;
  if(glXQueryExtension(display,&error,&event)==0) {
    CWarn ("X server does not have OpenGL extensions.\n");
  }}

  if(athis.doubleBufferOn==True) {
    int atbs[] = {
      GLX_RGBA,
      GLX_RED_SIZE,   1,
      GLX_GREEN_SIZE, 1,
      GLX_BLUE_SIZE,  1,
      GLX_ALPHA_SIZE, 1,
      GLX_DEPTH_SIZE, 1,
      GLX_DOUBLEBUFFER,
      None
    };
    vinfo = glXChooseVisual  (display,iscreen,atbs);
    if(vinfo==NULL) { /*GLX_ALPHA_SIZE : problem with Xming. Try without it.*/
      int atbs[] = {
        GLX_RGBA,
        GLX_RED_SIZE,   1,
        GLX_GREEN_SIZE, 1,
        GLX_BLUE_SIZE,  1,
        GLX_DEPTH_SIZE, 1,
        GLX_DOUBLEBUFFER,
        None
      };
      vinfo = glXChooseVisual  (display,iscreen,atbs);
    }
  } else {
    int atbs[] = {
      GLX_RGBA,
      GLX_RED_SIZE,   1,
      GLX_GREEN_SIZE, 1,
      GLX_BLUE_SIZE,  1,
      GLX_ALPHA_SIZE, 1,
      GLX_DEPTH_SIZE, 1,
      None
    };
    vinfo = glXChooseVisual  (display,iscreen,atbs);
    if(vinfo==NULL) { /*GLX_ALPHA_SIZE : problem with Xming. Try without it.*/
      int atbs[] = {
        GLX_RGBA,
        GLX_RED_SIZE,   1,
        GLX_GREEN_SIZE, 1,
        GLX_BLUE_SIZE,  1,
        GLX_DEPTH_SIZE, 1,
        None
      };
      vinfo = glXChooseVisual  (display,iscreen,atbs);
    }      
  }    

  if(vinfo==NULL) {
    CWarn ("Can't choose a visual.\n");
  } else {
    This->core.depth = vinfo->depth;
    athis.visual = vinfo->visual;
    if( (vinfo->depth ==DefaultDepth (display,iscreen))  &&
        (vinfo->visual==DefaultVisual(display,iscreen)) ) {
      This->core.colormap = XDefaultColormap (display,iscreen);
      athis.installColormap = False;
    } else {
      This->core.colormap =
	XCreateColormap (display,XRootWindow(display,iscreen),vinfo->visual, AllocNone); 
      athis.installColormap = True;
    }
    if(This->core.colormap==0L) {
      CWarn ("Can't get/create a X colormap.\n");
    }
    athis.glContext = glXCreateContext (display,vinfo,NULL,GL_FALSE);
    if(athis.glContext==NULL) {
      CWarn ("Can't create a GLX context.\n");
    }
    XFree(vinfo);
  }}

  XtAddEventHandler(This,ButtonPressMask|ButtonReleaseMask|ButtonMotionMask,0,EventHandler,NULL);

#ifdef DEBUG
  printf("debug : OpenGLArea : InitializeWidget : end\n");
#endif

  /*avoid C++ warnings*/
  a_request = NULL;
  a_args = NULL;
  a_argn = 0;
}

static void RealizeWidget(Widget This,XtValueMask* a_mask,XSetWindowAttributes* a_watbs) {
#ifdef DEBUG
  printf("debug : OpenGLArea : RealizeWidget : %s\n",XtName(This));
#endif

  /*Have to create window ourselves due to OpenGL that compells it's visual.*/
  /*In principle colormap is correctly set in a_watbs.*/

  XtCreateWindow(This,(unsigned int)InputOutput,athis.visual,*a_mask,a_watbs);

  /* Call the Realize procedure (XtInheritRealize) */
  if(openGLAreaWidgetClass->core_class.superclass->core_class.realize!=NULL)
    (openGLAreaWidgetClass->core_class.superclass->core_class.realize)(This, a_mask, a_watbs);

  /*If window is delete, all seems ok.*/
  if(athis.installColormap==True) XWidgetInstallColormap(This);

  /*MakeCurrent(This);*/

#ifdef DEBUG
  printf("debug : OpenGLArea : RealizeWidget : end\n");
#endif
}

static void DestroyWidget(Widget This) {
  if(athis.installColormap==True) {
    XWidgetUninstallColormap (This);
    athis.installColormap    = False;
    XFreeColormap(XtDisplay(This),This->core.colormap);
  }
  if(athis.glContext!=NULL) {
    glXMakeCurrent(XtDisplay(This),None,NULL);
    glXDestroyContext(XtDisplay(This),athis.glContext);
    athis.glContext = NULL;
  }
}

#define IFMOD(a_field)  if(athis.a_field != acur.a_field)
static Boolean SetValues(Widget a_current,Widget a_request,Widget This,ArgList a_args,Cardinal* a_argn) {
  IFMOD(doubleBufferOn) {
    /* Can't change buffering here if X window is created. 
       With OpenGL, buffering fix parameter of the X window.
       Buffering must be choosen before the execution of the 
       Realize method that create the window. */
    if(XtIsRealized(This) && (athis.installColormap==True)) {
      CWarn("Can't change buffering after \"realization\" of the widget.\n");
      athis.doubleBufferOn = acur.doubleBufferOn;
    }
  }
  /*avoid C++ warnings*/
  a_request = NULL;
  a_args = NULL;
  a_argn = 0;
  return False;
}

static void ChangeWidgetSize(Widget This) {
#ifdef DEBUG
  printf("debug : OpenGLArea : ChangeWidgetSize : %s\n",XtName(This));
#endif

  /* Call the Resize procedure (XtInheritResize) */
  if(openGLAreaWidgetClass->core_class.superclass->core_class.resize!=NULL)
    (openGLAreaWidgetClass->core_class.superclass->core_class.resize)(This);

#ifdef DEBUG
  printf("debug : OpenGLArea : ChangeWidgetSize : end\n");
#endif
}
     
static void DrawWidget(Widget  This,XEvent* a_event,Region a_region) {
#ifdef DEBUG
  printf("debug : OpenGLArea : DrawWidget : %s\n",XtName(This));
#endif

  if(openGLAreaWidgetClass->core_class.superclass->core_class.expose!=NULL)
    (openGLAreaWidgetClass->core_class.superclass->core_class.expose)(This,a_event,a_region);

  if(MakeCurrent(This)==1) {
#ifdef DEBUG
    printf("debug : OpenGLArea : DrawWidget : %s : MakeCurrent ok : call paintCallback...\n",XtName(This));
#endif
    XoAnyCallbackStruct value;
    value.reason = XoCR_PAINT;
    value.event = a_event;
    XtCallCallbacks(This,XoNpaintCallback,(XtPointer)&value);
    glXSwapBuffers(XtDisplay(This),XtWindow(This));
    glXMakeCurrent(XtDisplay(This),None,NULL);
  }

#ifdef DEBUG
  printf("debug : OpenGLArea : DrawWidget : end\n");
#endif

  /*avoid C++ warnings*/
  a_event = NULL;
  a_region = NULL;
}


/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

void OpenGLAreaPaint(Widget This) {
  if(!XtIsRealized(This)) return;
  if(MakeCurrent(This)==1) {
    XoAnyCallbackStruct value;
    value.reason = XoCR_PAINT;
    value.event = 0;
    XtCallCallbacks(This,XoNpaintCallback,(XtPointer)&value);
    glXSwapBuffers(XtDisplay(This),XtWindow(This));
    glXMakeCurrent(XtDisplay(This),None,NULL);
  }
}

#ifdef TOOLS_XT_OPENGLAREA_HAS_GL2PS

#include <tools/c_gl2ps.h>

int OpenGLAreaWrite_gl2ps(Widget This,const char* aFileName,const char* opts) {
  FILE* file;
  int options;
  int sort;
  inlib_GLint vp[4];
  int bufsize = 0;
  static inlib_gl2ps_gl_funcs_t s_OpenGL_funcs = {
    glIsEnabled,
    glBegin,
    glEnd,
    glGetFloatv,
    glVertex3f,
    glGetBooleanv,
    glGetIntegerv,
    glRenderMode,
    glFeedbackBuffer,
    glPassThrough
  };

  file = fopen(aFileName,"w");
  if(!file) return 1;

  inlib_c_gl2ps_set_gl_funcs(&s_OpenGL_funcs);

  options = TOOLS_GL2PS_OCCLUSION_CULL 
    | TOOLS_GL2PS_BEST_ROOT 
    | TOOLS_GL2PS_SILENT
    | TOOLS_GL2PS_DRAW_BACKGROUND;
  sort = TOOLS_GL2PS_BSP_SORT;
  //sort = TOOLS_GL2PS_SIMPLE_SORT;
    
  vp[0] = 0;
  vp[1] = 0;
  vp[2] = This->core.width;
  vp[3] = This->core.height;

  inlib_c_gl2psBeginPage("title","exlib_Xt_OpenGLArea", 
                 vp,TOOLS_GL2PS_EPS,sort,options, 
                 TOOLS_GL_RGBA,0, NULL,0,0,0,bufsize, 
                 file,aFileName);
    
 {XoAnyCallbackStruct value;
  value.reason = XoCR_PAINT;
  value.event = 0;
  XtCallCallbacks(This,XoNpaintCallback,(XtPointer)&value);}
    
  inlib_c_gl2psEndPage();

  inlib_c_gl2ps_reset_gl_funcs();

  fclose(file);

  return 0;
}

#else
#ifdef __cplusplus
int OpenGLAreaWrite_gl2ps(Widget,const char*,const char*) {return 1;}
#else
int OpenGLAreaWrite_gl2ps(Widget w,const char* f,const char* o) {return 1;}
#endif
#endif

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
static int MakeCurrent(Widget This) {
  if(This==NULL)          return 0;
  if(!XtIsRealized(This)) return 0;
  if(athis.glContext==NULL) return 0;
  return (int)glXMakeCurrent(XtDisplay(This),XtWindow(This),athis.glContext);
}

static void EventHandler(Widget This,XtPointer a_tag,XEvent* a_event ,Boolean* a_continue) {
  XoAnyCallbackStruct value;
  value.reason = XoCR_EVENT;
  value.event = a_event;
  XtCallCallbacks(This,XoNeventCallback,(XtPointer)&value);
  a_tag = NULL;
  a_continue = NULL;
}

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
static void XWidgetInstallColormap(Widget This) {
/* 
  From Mesa/widgets/GLwDrawingArea.c/post_colormap routine.
  Could use also XtSetWMColormapWindows.
*/
  Display*          display;
  XWindowAttributes watbs;
  Widget            shell;
  Window*           ws = NULL;
  int               wn = 0,found,count;
  Window            wshell,wthis;
  Colormap          cmapthis;
/*.........................................................................*/
  if(This==NULL) return;
  if( !XtIsWidget(This) || !XtIsRealized(This) ) return;
  shell = XWidgetGetShell (This);
  if(shell==NULL) return;
  display = XtDisplay (This);
  wthis   = XtWindow  (This);
  wshell  = XtWindow  (shell);
  XGetWMColormapWindows (display,wshell, &ws, &wn);
/* Check if colormap of this is a colormap of a Window in list.*/
  XGetWindowAttributes  (display,wthis,&watbs);
  cmapthis              = watbs.colormap;
  found                 = -1;
  for(count=0;count<wn;count++) {
    Colormap             cmap;
    XGetWindowAttributes (display,ws[count],&watbs);
    cmap                 = watbs.colormap;
    if(cmap==cmapthis) {
      XFree (ws);
      return;  /*done*/
    }
    if(ws[count]==wshell) {
      found = count;
    }
  }
  /*Have to add window of this in list.*/
  if(wn==0) {
    if(ws!=NULL) XFree(ws);
    ws = (Window*)malloc ( 2 * sizeof(Window));
  } else {
    ws = (Window*)realloc (ws,(wn + 2) * sizeof(Window));  
  }
  if(ws==NULL) return;
  if(found==-1) {
    /*Window of shell not in list.*/
    ws[wn] = wthis; wn++;
    ws[wn] = wshell;wn++;
  } else {
    ws[found] = wthis;
    ws[wn]    = wshell; wn++;  /*Shell must be last.*/
  }
  if (XSetWMColormapWindows(display,wshell, ws, wn)==0) {
    CWarnF ("XWidgetInstallColormap: can't install colormap of %s in %s.\n",XtName(This),XtName(shell));
  }
  XFree (ws);
}

static void XWidgetUninstallColormap(Widget This) {
  int               count;
  Widget            shell;
  Display*          display;
  Window            wthis,wshell;
  Window*           ws  = NULL;
  int               wn  = 0;
  Window*           nws = NULL;
  int               nwn = 0;
/*.........................................................................*/
  if(This==NULL) return;
  if( !XtIsWidget(This) || !XtIsRealized(This) ) return;
  shell = XWidgetGetShell (This);
  if(shell==NULL) return;
  display               = XtDisplay (This);
  wthis                 = XtWindow  (This);
  wshell                = XtWindow  (shell);
  XGetWMColormapWindows (display,wshell, &ws, &wn);
  if( (wn==0) || (ws==NULL) ) return;
  nws = (Window*)malloc ( wn  * sizeof(Window));
  if(nws==NULL) {
    XFree (ws);
    return;
  }
  nwn = 0;
  for(count=0;count<wn;count++) {
    if(ws[count]!=wthis) {
      nws[nwn] = ws[count];
      nwn++;
    }
  }
  if(wn!=nwn) {
    if (XSetWMColormapWindows(display,wshell, nws, nwn)==0) {
      CWarnF("XWidgetUninstallColormap: can't install colormap of %s in %s.\n",XtName(This),XtName(shell));
    }
  }
  XFree (ws);
  XFree (nws);
}

Widget XWidgetGetShell(Widget This) {
  Widget widget;
  if(This==NULL) return NULL;
  widget = This;
  while(1) {
    if(widget==NULL) return NULL;
    if(XtIsShell(widget)) return widget;
    widget = XtParent(widget);
  }
}

/*
exlib_build_use Xt inlib
*/
