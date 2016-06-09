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
/*
 * GL2PS, an OpenGL to PostScript Printing Library
 * Copyright (C) 1999-2003  Christophe Geuzaine
 *
 * $Id: gl2ps.h,v 1.5 2006/06/29 21:23:11 gunter Exp $
 *
 * E-mail: geuz@geuz.org
 * URL: http://www.geuz.org/gl2ps/
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#ifndef __GL2PS_H__
#define __GL2PS_H__

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

/* To generate a Windows dll, you have to define GL2PSDLL at compile
   time */

#ifdef WIN32
#  include <windows.h>
#  ifdef GL2PSDLL
#    ifdef GL2PSDLL_EXPORTS
#      define GL2PSDLL_API __declspec(dllexport)
#    else
#      define GL2PSDLL_API __declspec(dllimport)
#    endif
#  else
#    define GL2PSDLL_API
#  endif
#else
#  define GL2PSDLL_API
#endif

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#else
#  include <GL/gl.h>
#endif


#define GL2PS_VERSION                    0.8
#define GL2PS_NONE                       0

/* Output file format */

#define GL2PS_PS                         1
#define GL2PS_EPS                        2
#define GL2PS_TEX                        3

/* Sorting algorithms */

#define GL2PS_NO_SORT                    1
#define GL2PS_SIMPLE_SORT                2
#define GL2PS_BSP_SORT                   3

/* Options for gl2psBeginPage */

#define GL2PS_DRAW_BACKGROUND            (1<<0)
#define GL2PS_SIMPLE_LINE_OFFSET         (1<<1)
#define GL2PS_SILENT                     (1<<2)
#define GL2PS_BEST_ROOT                  (1<<3)
#define GL2PS_OCCLUSION_CULL             (1<<4)
#define GL2PS_NO_TEXT                    (1<<5)
#define GL2PS_LANDSCAPE                  (1<<6)
#define GL2PS_NO_PS3_SHADING             (1<<7)
#define GL2PS_NO_PIXMAP                  (1<<8)

/* Arguments for gl2psEnable/gl2psDisable */

#define GL2PS_POLYGON_OFFSET_FILL        1
#define GL2PS_POLYGON_BOUNDARY           2
#define GL2PS_LINE_STIPPLE               3

/* Magic numbers */

#define GL2PS_EPSILON                    5.e-3
#define GL2PS_DEPTH_FACT                 1000.0
#define GL2PS_SIMPLE_OFFSET              0.05
#define GL2PS_SIMPLE_OFFSET_LARGE        1.0
#define GL2PS_ZERO(arg)                  (std::fabs(arg)<1.e-20)
/*#define GL2PS_ZERO(arg)                ((arg)==0.0)*/

/* Message levels and error codes */

#define GL2PS_SUCCESS                    0
#define GL2PS_INFO                       1
#define GL2PS_WARNING                    2
#define GL2PS_ERROR                      3
#define GL2PS_NO_FEEDBACK                4
#define GL2PS_OVERFLOW                   5
#define GL2PS_UNINITIALIZED              6

/* Primitive types */

#define GL2PS_TEXT                       1
#define GL2PS_POINT                      2
#define GL2PS_LINE                       3
#define GL2PS_QUADRANGLE                 4
#define GL2PS_TRIANGLE                   5
#define GL2PS_PIXMAP                     6

/* BSP tree primitive comparison */

#define GL2PS_COINCIDENT                 1
#define GL2PS_IN_FRONT_OF                2
#define GL2PS_IN_BACK_OF                 3
#define GL2PS_SPANNING                   4

/* 2D BSP tree primitive comparison */

#define GL2PS_POINT_COINCIDENT           0
#define GL2PS_POINT_INFRONT              1
#define GL2PS_POINT_BACK                 2

/* Pass through options */

#define GL2PS_BEGIN_POLYGON_OFFSET_FILL  1
#define GL2PS_END_POLYGON_OFFSET_FILL    2
#define GL2PS_BEGIN_POLYGON_BOUNDARY     3
#define GL2PS_END_POLYGON_BOUNDARY       4
#define GL2PS_BEGIN_LINE_STIPPLE         5
#define GL2PS_END_LINE_STIPPLE           6
#define GL2PS_SET_POINT_SIZE             7
#define GL2PS_SET_LINE_WIDTH             8

typedef GLfloat GL2PSrgba[4];
typedef GLfloat GL2PSxyz[3];
typedef GLfloat GL2PSplane[4];

typedef struct _GL2PSbsptree2d GL2PSbsptree2d;

struct _GL2PSbsptree2d {
  GL2PSplane plane;
  GL2PSbsptree2d *front, *back;
};

typedef struct {
  GLint nmax, size, incr, n;
  char *array;
} GL2PSlist;

typedef struct _GL2PSbsptree GL2PSbsptree;

struct _GL2PSbsptree {
  GL2PSplane plane;
  GL2PSlist *primitives;
  GL2PSbsptree *front, *back;
};

typedef struct {
  GL2PSxyz xyz;
  GL2PSrgba rgba;
} GL2PSvertex;

typedef struct {
  GLshort fontsize;
  char *str, *fontname;
} GL2PSstring;

typedef struct {
  GLsizei width, height;
  GLenum format, type;
  GLfloat *pixels;
} GL2PSimage;

typedef struct {
  GLshort type, numverts;
  char boundary, dash, culled;
  GLfloat width, depth;
  GL2PSvertex *verts;
  GL2PSstring *text;
  GL2PSimage *image;
} GL2PSprimitive;

typedef struct {
  GLint format, sort, options, colorsize, colormode, buffersize, maxbestroot;
  const char *title, *producer, *filename;
  GLboolean shade, boundary;
  GLfloat *feedback, offset[2];
  GLint viewport[4];
  GL2PSrgba *colormap, lastrgba, threshold;
  float lastlinewidth;
  GL2PSlist *primitives;
  GL2PSbsptree2d *imagetree;
  FILE *stream;
} GL2PScontext;

/* public functions */

#ifdef __cplusplus
extern "C" {
#endif

GL2PSDLL_API GLint gl2psBeginPage(const char *title, const char *producer, 
				  GLint viewport[4], GLint format, GLint sort,
				  GLint options, GLint colormode,
				  GLint colorsize, GL2PSrgba *colormap, 
				  GLint nr, GLint ng, GLint nb, GLint buffersize,
				  FILE *stream, const char *filename);
GL2PSDLL_API GLint gl2psEndPage(void);
GL2PSDLL_API GLint gl2psBeginViewport(GLint viewport[4]);
GL2PSDLL_API GLint gl2psEndViewport(void);
GL2PSDLL_API GLint gl2psText(const char *str, const char *fontname, GLshort fontsize);
GL2PSDLL_API GLint gl2psDrawPixels(GLsizei width, GLsizei height,
				   GLint xorig, GLint yorig,
				   GLenum format, GLenum type, const void *pixels);
GL2PSDLL_API GLint gl2psEnable(GLint mode);
GL2PSDLL_API GLint gl2psDisable(GLint mode);
GL2PSDLL_API GLint gl2psPointSize(GLfloat value);
GL2PSDLL_API GLint gl2psLineWidth(GLfloat value);

#ifdef __cplusplus
}
#endif

#endif /* __GL2PS_H__ */
