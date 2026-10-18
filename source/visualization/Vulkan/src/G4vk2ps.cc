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
//
//
// Alex Vergara Gil  17th October 2026
// Ported to Vulkan from the Geant4 OpenGL visualization driver.

#include "G4vk2ps.hh"

#include <tools/gl2ps>

#include <cstdlib>
#include <cstring>
#include <limits>

G4vk2ps::G4vk2ps()
{
  fContext = 0;
  fVulkanFuncs.m_glIsEnabled = tools_dummy_glIsEnabled;
  fVulkanFuncs.m_glBegin = tools_dummy_glBegin;
  fVulkanFuncs.m_glEnd = tools_dummy_glEnd;
  fVulkanFuncs.m_glGetFloatv = tools_dummy_glGetFloatv;
  fVulkanFuncs.m_glVertex3f = tools_dummy_glVertex3f;
  fVulkanFuncs.m_glGetBooleanv = tools_dummy_glGetBooleanv;
  fVulkanFuncs.m_glGetIntegerv = tools_dummy_glGetIntegerv;
  fVulkanFuncs.m_glRenderMode = tools_dummy_glRenderMode;
  fVulkanFuncs.m_glFeedbackBuffer = tools_dummy_glFeedbackBuffer;
  fVulkanFuncs.m_glPassThrough = tools_dummy_glPassThrough;

  fFile = 0;
  fViewport[0] = 0;
  fViewport[1] = 0;
  fViewport[2] = 0;
  fViewport[3] = 0;
  fBufferSize = 2048;
  fBufferSizeLimit = (std::numeric_limits<int>::max)();
  fExportImageFormat = TOOLS_GL2PS_PDF;
  resetBufferSizeParameters();
}

G4vk2ps::~G4vk2ps()
{
  if (fFile)
  {
    ::fclose(fFile);
    fFile = 0;
  }
  if (fContext)
  {
    ::tools_gl2psDeleteContext(fContext);
    fContext = 0;
  }
}

void G4vk2ps::setVulkanFunctions(tools_gl2ps_gl_funcs_t* a_funcs)
{
  fVulkanFuncs.m_glIsEnabled = a_funcs->m_glIsEnabled;
  fVulkanFuncs.m_glBegin = a_funcs->m_glBegin;
  fVulkanFuncs.m_glEnd = a_funcs->m_glEnd;
  fVulkanFuncs.m_glGetFloatv = a_funcs->m_glGetFloatv;
  fVulkanFuncs.m_glVertex3f = a_funcs->m_glVertex3f;
  fVulkanFuncs.m_glGetBooleanv = a_funcs->m_glGetBooleanv;
  fVulkanFuncs.m_glGetIntegerv = a_funcs->m_glGetIntegerv;
  fVulkanFuncs.m_glRenderMode = a_funcs->m_glRenderMode;
  fVulkanFuncs.m_glFeedbackBuffer = a_funcs->m_glFeedbackBuffer;
  fVulkanFuncs.m_glPassThrough = a_funcs->m_glPassThrough;
}

void G4vk2ps::resetBufferSizeParameters()
{
  fBufferSize = 2048;
}

void G4vk2ps::setLineWidth(int width)
{
  if (!fContext) return;
  ::tools_gl2psLineWidth(fContext, width);
}

void G4vk2ps::setPointSize(int size)
{
  if (!fContext) return;
  ::tools_gl2psPointSize(fContext, size);
}

void G4vk2ps::addTextOpt(const char* str, const char* fontname, tools_GLshort fontsize,
                         tools_GLint alignment, tools_GLfloat angle)
{
  if (!fContext) return;
  ::tools_gl2psTextOpt(fContext, str, fontname, fontsize, alignment, angle);
}

void G4vk2ps::setViewport(int a, int b, int winSizeX, int winSizeY)
{
  fViewport[0] = a;
  fViewport[1] = b;
  fViewport[2] = winSizeX;
  fViewport[3] = winSizeY;
}

void G4vk2ps::setFileName(const char* aFileName)
{
  fFileName = aFileName;
}

bool G4vk2ps::enableFileWriting()
{
  if (fFile)
  {
    ::fclose(fFile);
    fFile = 0;
  }
  if (fContext)
  {
    ::tools_gl2psDeleteContext(fContext);
    fContext = 0;
  }

  fContext = ::tools_gl2psCreateContext();
  if (!fContext) return false;

  ::tools_gl2ps_set_gl_funcs(fContext, &fVulkanFuncs);

  fFile = ::fopen(fFileName, "wb");
  if (!fFile)
  {
    ::tools_gl2psDeleteContext(fContext);
    fContext = 0;
    return false;
  }

  // No buffering for output file
  setvbuf(fFile, NULL, _IONBF, 2048);

  return true;
}

void G4vk2ps::disableFileWriting()
{
  if (fFile)
  {
    ::fclose(fFile);
    fFile = 0;
  }
  if (fContext)
  {
    ::tools_gl2psDeleteContext(fContext);
    fContext = 0;
  }
}

bool G4vk2ps::fileWritingEnabled() const
{
  return (fContext && fFile ? true : false);
}

bool G4vk2ps::extendBufferSize()
{
  // extend buffer size *2
  if (fBufferSize < (fBufferSizeLimit / 2))
  {
    fBufferSize = fBufferSize * 2;
    return true;
  }
  return false;
}

// FWJ
void G4vk2ps::setBufferSize(int newSize)
{
  fBufferSize = (newSize < int(fBufferSizeLimit)) ? newSize : fBufferSizeLimit;
}

bool G4vk2ps::beginPage()
{
  if (!fContext) return false;
  if (!fFile) return false;

  if ((fViewport[2] <= 0) || (fViewport[3] <= 0)) return false;

  int options =
    TOOLS_GL2PS_BEST_ROOT | TOOLS_GL2PS_DRAW_BACKGROUND | TOOLS_GL2PS_USE_CURRENT_VIEWPORT;
  int sort = TOOLS_GL2PS_BSP_SORT;

  tools_GLint res = ::tools_gl2psBeginPage(fContext, "Geant4 output", "Geant4", fViewport,
                                           fExportImageFormat, sort, options, TOOLS_GL_RGBA, 0,
                                           NULL, 0, 0, 0, fBufferSize, fFile, fFileName.c_str());
  if (res == TOOLS_GL2PS_ERROR) return false;

  // enable blending for all
  ::tools_gl2psEnable(fContext, TOOLS_GL2PS_BLEND);

  return true;
}

bool G4vk2ps::endPage()
{
  int _status = 0;
  if (fContext)
  {
    _status = ::tools_gl2psEndPage(fContext);
  }
  if (_status == TOOLS_GL2PS_OVERFLOW) return false;
  return true;
}

void G4vk2ps::setExportImageFormat(unsigned int type)
{
  if (fFile) return;
  fExportImageFormat = type;
}
