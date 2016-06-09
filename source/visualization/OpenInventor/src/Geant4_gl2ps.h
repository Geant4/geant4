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
#ifndef Geant4_gl2ps_h
#define Geant4_gl2ps_h

// gl2ps-0.73.
// The gl2ps code is prefixed by Geant4_ in order
// to avoid clashes with other gl2ps code that may come at link
// time from other channels.

#define gl2psBeginPage Geant4_gl2psBeginPage
#define gl2psEndPage Geant4_gl2psEndPage
#define gl2psText Geant4_gl2psText
#define gl2psEnable Geant4_gl2psEnable
#define gl2psDisable Geant4_gl2psDisable
#define gl2psPointSize Geant4_gl2psPointSize
#define gl2psLineWidth Geant4_gl2psLineWidth
#define gl2psNumShadeColors Geant4_gl2psNumShadeColors
#define gl2psDrawPixels Geant4_gl2psDrawPixels
#define gl2psBeginViewport Geant4_gl2psBeginViewport
#define gl2psEndViewport Geant4_gl2psEndViewport

#define gl2psMsg Geant4_gl2psMsg
#define gl2psMalloc Geant4_gl2psMalloc
#define gl2psRealloc Geant4_gl2psRealloc
#define gl2psFree Geant4_gl2psFree
#define gl2psListRealloc Geant4_gl2psListRealloc
#define gl2psListCreate Geant4_gl2psListCreate
#define gl2psListDelete Geant4_gl2psListDelete
#define gl2psListAdd Geant4_gl2psListAdd
#define gl2psListNbr Geant4_gl2psListNbr
#define gl2psListPointer Geant4_gl2psListPointer
#define gl2psListSort Geant4_gl2psListSort
#define gl2psListAction Geant4_gl2psListAction
#define gl2psListActionInverse Geant4_gl2psListActionInverse
#define gl2psComparePointPlane Geant4_gl2psComparePointPlane
#define gl2psPsca Geant4_gl2psPsca
#define gl2psPvec Geant4_gl2psPvec
#define gl2psNorm Geant4_gl2psNorm
#define gl2psGetNormal Geant4_gl2psGetNormal
#define gl2psGetPlane Geant4_gl2psGetPlane
#define gl2psCutEdge Geant4_gl2psCutEdge
#define gl2psCreateSplitPrimitive Geant4_gl2psCreateSplitPrimitive
#define gl2psAddIndex Geant4_gl2psAddIndex
#define gl2psGetIndex Geant4_gl2psGetIndex
#define gl2psTestSplitPrimitive Geant4_gl2psTestSplitPrimitive
#define gl2psSplitPrimitive Geant4_gl2psSplitPrimitive
#define gl2psDivideQuad Geant4_gl2psDivideQuad
#define gl2psCompareDepth Geant4_gl2psCompareDepth
#define gl2psTrianglesFirst Geant4_gl2psTrianglesFirst
#define gl2psFindRoot Geant4_gl2psFindRoot
#define gl2psFreePrimitive Geant4_gl2psFreePrimitive
#define gl2psAddPrimitiveInList Geant4_gl2psAddPrimitiveInList
#define gl2psFreeBspTree Geant4_gl2psFreeBspTree
#define gl2psGreater Geant4_gl2psGreater
#define gl2psLess Geant4_gl2psLess
#define gl2psBuildBspTree Geant4_gl2psBuildBspTree
#define gl2psTraverseBspTree Geant4_gl2psTraverseBspTree
#define gl2psGetPlaneFromPoints Geant4_gl2psGetPlaneFromPoints
#define gl2psFreeBspImageTree Geant4_gl2psFreeBspImageTree
#define gl2psCheckPoint Geant4_gl2psCheckPoint
#define gl2psAddPlanesInBspTreeImage Geant4_gl2psAddPlanesInBspTreeImage
#define gl2psCheckPrimitive Geant4_gl2psCheckPrimitive
#define gl2psCreateSplitPrimitive2D Geant4_gl2psCreateSplitPrimitive2D
#define gl2psSplitPrimitive2D Geant4_gl2psSplitPrimitive2D
#define gl2psAddInImageTree Geant4_gl2psAddInImageTree
#define gl2psAddInImage Geant4_gl2psAddInImage
#define gl2psAddBoundaryInList Geant4_gl2psAddBoundaryInList
#define gl2psBuildPolygonBoundary Geant4_gl2psBuildPolygonBoundary
#define gl2psAddPolyPrimitive Geant4_gl2psAddPolyPrimitive
#define gl2psGetVertex Geant4_gl2psGetVertex
#define gl2psParseFeedbackBuffer Geant4_gl2psParseFeedbackBuffer
#define gl2psSameColor Geant4_gl2psSameColor
#define gl2psVertsSameColor Geant4_gl2psVertsSameColor
#define gl2psPrintPostScriptHeader Geant4_gl2psPrintPostScriptHeader
#define gl2psPrintPostScriptColor Geant4_gl2psPrintPostScriptColor
#define gl2psResetPostScriptColor Geant4_gl2psResetPostScriptColor
#define gl2psPrintPostScriptPrimitive Geant4_gl2psPrintPostScriptPrimitive
#define gl2psPrintPostScriptFooter Geant4_gl2psPrintPostScriptFooter
#define gl2psPrintTeXHeader Geant4_gl2psPrintTeXHeader
#define gl2psPrintTeXPrimitive Geant4_gl2psPrintTeXPrimitive
#define gl2psPrintTeXFooter Geant4_gl2psPrintTeXFooter

#define gl2psPrintPostScriptBeginViewport Geant4_gl2psPrintPostScriptBeginViewport
#define gl2psPrintPostScriptEndViewport Geant4_gl2psPrintPostScriptEndViewport

#define gl2psWriteByte Geant4_gl2psWriteByte
#define gl2psGetRGB Geant4_gl2psGetRGB
#define gl2psPrintPostScriptPixmap Geant4_gl2psPrintPostScriptPixmap
#define gl2psListReset Geant4_gl2psListReset
#define gl2psAddInBspImageTree Geant4_gl2psAddInBspImageTree
#define gl2psPrintPrimitives Geant4_gl2psPrintPrimitives
#define gl2psPrintTeXBeginViewport Geant4_gl2psPrintTeXBeginViewport
#define gl2psPrintTeXEndViewport Geant4_gl2psPrintTeXEndViewport

#define gl2ps Geant4_gl2ps

#include "gl2ps.h"

#endif
