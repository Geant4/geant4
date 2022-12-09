/* Copyright (C) 2010, Guy Barrand. All rights reserved. */
/* See the file tools.license for terms.                 */

#ifndef tools_gl2ps_def_h
#define tools_gl2ps_def_h

typedef int            tools_GLint;
typedef unsigned int   tools_GLuint;
typedef float          tools_GLfloat;
typedef unsigned int   tools_GLenum;
typedef short          tools_GLshort;
typedef unsigned short tools_GLushort;
typedef int            tools_GLsizei;
typedef unsigned char  tools_GLboolean;

/*----------------------------------------------------------*/
/*---  from gl2ps.h : --------------------------------------*/
/*----------------------------------------------------------*/
#define TOOLS_GL2PSDLL_API inline

#define TOOLS_GL2PS_MAJOR_VERSION 1
#define TOOLS_GL2PS_MINOR_VERSION 4
#define TOOLS_GL2PS_PATCH_VERSION 2
#define TOOLS_GL2PS_EXTRA_VERSION ""

#define TOOLS_GL2PS_VERSION (TOOLS_GL2PS_MAJOR_VERSION + \
                       0.01 * TOOLS_GL2PS_MINOR_VERSION + \
                       0.0001 * TOOLS_GL2PS_PATCH_VERSION)

#define TOOLS_GL2PS_COPYRIGHT "(C) 1999-2020 C. Geuzaine"

/* Output file formats (the values and the ordering are important!) */

#define TOOLS_GL2PS_PS  0
#define TOOLS_GL2PS_EPS 1
#define TOOLS_GL2PS_TEX 2
#define TOOLS_GL2PS_PDF 3
#define TOOLS_GL2PS_SVG 4
#define TOOLS_GL2PS_PGF 5

/* Sorting algorithms */

#define TOOLS_GL2PS_NO_SORT     1
#define TOOLS_GL2PS_SIMPLE_SORT 2
#define TOOLS_GL2PS_BSP_SORT    3

/* Message levels and error codes */

#define TOOLS_GL2PS_SUCCESS       0
#define TOOLS_GL2PS_INFO          1
#define TOOLS_GL2PS_WARNING       2
#define TOOLS_GL2PS_ERROR         3
#define TOOLS_GL2PS_NO_FEEDBACK   4
#define TOOLS_GL2PS_OVERFLOW      5
#define TOOLS_GL2PS_UNINITIALIZED 6

/* Options for tools_gl2psBeginPage */

#define TOOLS_GL2PS_NONE                 0
#define TOOLS_GL2PS_DRAW_BACKGROUND      (1<<0)
#define TOOLS_GL2PS_SIMPLE_LINE_OFFSET   (1<<1)
#define TOOLS_GL2PS_SILENT               (1<<2)
#define TOOLS_GL2PS_BEST_ROOT            (1<<3)
#define TOOLS_GL2PS_OCCLUSION_CULL       (1<<4)
#define TOOLS_GL2PS_NO_TEXT              (1<<5)
#define TOOLS_GL2PS_LANDSCAPE            (1<<6)
#define TOOLS_GL2PS_NO_PS3_SHADING       (1<<7)
#define TOOLS_GL2PS_NO_PIXMAP            (1<<8)
#define TOOLS_GL2PS_USE_CURRENT_VIEWPORT (1<<9)
#define TOOLS_GL2PS_COMPRESS             (1<<10)
#define TOOLS_GL2PS_NO_BLENDING          (1<<11)
#define TOOLS_GL2PS_TIGHT_BOUNDING_BOX   (1<<12)
#define TOOLS_GL2PS_NO_OPENGL_CONTEXT    (1<<13)
#define TOOLS_GL2PS_NO_TEX_FONTSIZE      (1<<14)
#define TOOLS_GL2PS_PORTABLE_SORT        (1<<15)

/* Arguments for tools_gl2psEnable/tools_gl2psDisable */

#define TOOLS_GL2PS_POLYGON_OFFSET_FILL 1
#define TOOLS_GL2PS_POLYGON_BOUNDARY    2
#define TOOLS_GL2PS_LINE_STIPPLE        3
#define TOOLS_GL2PS_BLEND               4


/* Arguments for tools_gl2psLineCap/Join */

#define TOOLS_GL2PS_LINE_CAP_BUTT       0
#define TOOLS_GL2PS_LINE_CAP_ROUND      1
#define TOOLS_GL2PS_LINE_CAP_SQUARE     2

#define TOOLS_GL2PS_LINE_JOIN_MITER     0
#define TOOLS_GL2PS_LINE_JOIN_ROUND     1
#define TOOLS_GL2PS_LINE_JOIN_BEVEL     2

/* Text alignment (o=raster position; default mode is BL):
   +---+ +---+ +---+ +---+ +---+ +---+ +-o-+ o---+ +---o
   | o | o   | |   o |   | |   | |   | |   | |   | |   |
   +---+ +---+ +---+ +-o-+ o---+ +---o +---+ +---+ +---+
    C     CL    CR    B     BL    BR    T     TL    TR */

#define TOOLS_GL2PS_TEXT_C  1
#define TOOLS_GL2PS_TEXT_CL 2
#define TOOLS_GL2PS_TEXT_CR 3
#define TOOLS_GL2PS_TEXT_B  4
#define TOOLS_GL2PS_TEXT_BL 5
#define TOOLS_GL2PS_TEXT_BR 6
#define TOOLS_GL2PS_TEXT_T  7
#define TOOLS_GL2PS_TEXT_TL 8
#define TOOLS_GL2PS_TEXT_TR 9

typedef tools_GLfloat tools_GL2PSrgba[4];
typedef tools_GLfloat tools_GL2PSxyz[3];

typedef struct {
  tools_GL2PSxyz xyz;
  tools_GL2PSrgba rgba;
} tools_GL2PSvertex;

/* Primitive types */
#define TOOLS_GL2PS_NO_TYPE          -1
#define TOOLS_GL2PS_TEXT             1
#define TOOLS_GL2PS_POINT            2
#define TOOLS_GL2PS_LINE             3
#define TOOLS_GL2PS_QUADRANGLE       4
#define TOOLS_GL2PS_TRIANGLE         5
#define TOOLS_GL2PS_PIXMAP           6
#define TOOLS_GL2PS_IMAGEMAP         7
#define TOOLS_GL2PS_IMAGEMAP_WRITTEN 8
#define TOOLS_GL2PS_IMAGEMAP_VISIBLE 9
#define TOOLS_GL2PS_SPECIAL          10

/*----------------------------------------------------------*/
/*--- from OpenGL : ----------------------------------------*/
/*----------------------------------------------------------*/
#define TOOLS_GL_TRUE                           1
#define TOOLS_GL_FALSE                          0

#define TOOLS_GL_FLOAT				0x1406
#define TOOLS_GL_BLEND				0x0BE2

#define TOOLS_GL_SRC_ALPHA				0x0302
#define TOOLS_GL_ONE_MINUS_SRC_ALPHA                  0x0303

#define TOOLS_GL_RGB					0x1907
#define TOOLS_GL_RGBA                                 0x1908

#define TOOLS_GL_POINTS                               0x0000

#define TOOLS_GL_CURRENT_RASTER_POSITION_VALID	0x0B08
#define TOOLS_GL_CURRENT_RASTER_POSITION		0x0B07
#define TOOLS_GL_CURRENT_RASTER_COLOR			0x0B04
#define TOOLS_GL_ZERO					0
#define TOOLS_GL_ONE					1
#define TOOLS_GL_COLOR_INDEX				0x1900

#define TOOLS_GL_POINT_TOKEN				0x0701
#define TOOLS_GL_LINE_TOKEN				0x0702
#define TOOLS_GL_LINE_RESET_TOKEN			0x0707
#define TOOLS_GL_POLYGON_TOKEN			0x0703
#define TOOLS_GL_BITMAP_TOKEN				0x0704
#define TOOLS_GL_DRAW_PIXEL_TOKEN			0x0705

#define TOOLS_GL_COPY_PIXEL_TOKEN			0x0706
#define TOOLS_GL_PASS_THROUGH_TOKEN			0x0700

#define TOOLS_GL_FEEDBACK				0x1C01
#define TOOLS_GL_COLOR_CLEAR_VALUE			0x0C22
#define TOOLS_GL_INDEX_CLEAR_VALUE			0x0C20
#define TOOLS_GL_RENDER				0x1C00
#define TOOLS_GL_VIEWPORT				0x0BA2
#define TOOLS_GL_BLEND_SRC				0x0BE1
#define TOOLS_GL_BLEND_DST				0x0BE0
#define TOOLS_GL_3D_COLOR				0x0602

#define TOOLS_GL_POLYGON_OFFSET_FACTOR		0x8038
#define TOOLS_GL_POLYGON_OFFSET_UNITS			0x2A00
#define TOOLS_GL_LINE_STIPPLE_PATTERN			0x0B25
#define TOOLS_GL_LINE_STIPPLE_REPEAT			0x0B26

#define TOOLS_GL_ZOOM_X				0x0D16
#define TOOLS_GL_ZOOM_Y				0x0D17

/*----------------------------------------------------------*/
/*----------------------------------------------------------*/
/*----------------------------------------------------------*/
typedef struct tools_GL2PScontextRec* tools_GL2PScontextPointer;

typedef tools_GLboolean (*tools_glIsEnabled_func)      (tools_GLenum);
typedef void            (*tools_glBegin_func)          (tools_GLenum);
typedef void            (*tools_glEnd_func)            ();
typedef void            (*tools_glGetFloatv_func)      (tools_GLenum,tools_GLfloat*);
typedef void            (*tools_glVertex3f_func)       (tools_GLfloat,tools_GLfloat,tools_GLfloat);
typedef void            (*tools_glGetBooleanv_func)    (tools_GLenum,tools_GLboolean*);
typedef void            (*tools_glGetIntegerv_func)    (tools_GLenum,tools_GLint*);
typedef tools_GLint     (*tools_glRenderMode_func)     (tools_GLenum);
typedef void            (*tools_glFeedbackBuffer_func) (tools_GLsizei,tools_GLenum,tools_GLfloat*);
typedef void            (*tools_glPassThrough_func)    (tools_GLfloat);

typedef struct {
  tools_glIsEnabled_func      m_glIsEnabled;
  tools_glBegin_func          m_glBegin;
  tools_glEnd_func            m_glEnd;
  tools_glGetFloatv_func      m_glGetFloatv;
  tools_glVertex3f_func       m_glVertex3f;
  tools_glGetBooleanv_func    m_glGetBooleanv;
  tools_glGetIntegerv_func    m_glGetIntegerv;
  tools_glRenderMode_func     m_glRenderMode;
  tools_glFeedbackBuffer_func m_glFeedbackBuffer;
  tools_glPassThrough_func    m_glPassThrough;
} tools_gl2ps_gl_funcs_t;

#endif
