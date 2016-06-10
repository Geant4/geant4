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

// gl2ps-1.3.5
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
#define gl2psDrawPixels Geant4_gl2psDrawPixels
#define gl2psBeginViewport Geant4_gl2psBeginViewport
#define gl2psEndViewport Geant4_gl2psEndViewport
#define gl2psTextOpt Geant4_gl2psTextOpt
#define gl2psSetOptions Geant4_gl2psSetOptions
#define gl2psGetOptions Geant4_gl2psGetOptions
#define gl2psSpecial Geant4_gl2psSpecial
#define gl2psBlendFunc Geant4_gl2psBlendFunc
#define gl2psDrawImageMap Geant4_gl2psDrawImageMap
#define gl2psGetFileExtension Geant4_gl2psGetFileExtension
#define gl2psGetFormatDescription Geant4_gl2psGetFormatDescription

#define gl2psMsg Geant4_gl2psMsg
#define gl2psMalloc Geant4_gl2psMalloc
#define gl2psRealloc Geant4_gl2psRealloc
#define gl2psFree Geant4_gl2psFree
#define gl2psWriteBigEndian Geant4_gl2psWriteBigEndian
#define gl2psSetupCompress Geant4_gl2psSetupCompress
#define gl2psFreeCompress Geant4_gl2psFreeCompress
#define gl2psAllocCompress Geant4_gl2psAllocCompress
#define gl2psReallocCompress Geant4_gl2psReallocCompress
#define gl2psWriteBigEndianCompress Geant4_gl2psWriteBigEndianCompress
#define gl2psDeflate Geant4_gl2psDeflate
#define gl2psPrintf Geant4_gl2psPrintf
#define gl2psPrintGzipHeader Geant4_gl2psPrintGzipHeader
#define gl2psPrintGzipFooter Geant4_gl2psPrintGzipFooter
#define gl2psListReset Geant4_gl2psListReset
#define gl2psListRealloc Geant4_gl2psListRealloc
#define gl2psListCreate Geant4_gl2psListCreate
#define gl2psListDelete Geant4_gl2psListDelete
#define gl2psListAdd Geant4_gl2psListAdd
#define gl2psListNbr Geant4_gl2psListNbr
#define gl2psListPointer Geant4_gl2psListPointer
#define gl2psListSort Geant4_gl2psListSort
#define gl2psListAction Geant4_gl2psListAction
#define gl2psListActionInverse Geant4_gl2psListActionInverse
#define gl2psListRead Geant4_gl2psListRead
#define gl2psEncodeBase64Block Geant4_gl2psEncodeBase64Block
#define gl2psListEncodeBase64 Geant4_gl2psListEncodeBase64
#define gl2psSameColor Geant4_gl2psSameColor
#define gl2psVertsSameColor Geant4_gl2psVertsSameColor
#define gl2psSameColorThreshold Geant4_gl2psSameColorThreshold
#define gl2psSetLastColor Geant4_gl2psSetLastColor
#define gl2psGetRGB Geant4_gl2psGetRGB
#define gl2psCopyPixmap Geant4_gl2psCopyPixmap
#define gl2psFreePixmap Geant4_gl2psFreePixmap
#define gl2psUserWritePNG Geant4_gl2psUserWritePNG
#define gl2psUserFlushPNG Geant4_gl2psUserFlushPNG
#define gl2psConvertPixmapToPNG Geant4_gl2psConvertPixmapToPNG
#define gl2psAddText Geant4_gl2psAddText
#define gl2psCopyText Geant4_gl2psCopyText
#define gl2psFreeText Geant4_gl2psFreeText
#define gl2psSupportedBlendMode Geant4_gl2psSupportedBlendMode
#define gl2psAdaptVertexForBlending Geant4_gl2psAdaptVertexForBlending
#define gl2psAssignTriangleProperties Geant4_gl2psAssignTriangleProperties
#define gl2psFillTriangleFromPrimitive Geant4_gl2psFillTriangleFromPrimitive
#define gl2psInitTriangle Geant4_gl2psInitTriangle
#define gl2psCopyPrimitive Geant4_gl2psCopyPrimitive
#define gl2psSamePosition Geant4_gl2psSamePosition
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
#define gl2psFreeImagemap Geant4_gl2psFreeImagemap
#define gl2psFreePrimitive Geant4_gl2psFreePrimitive
#define gl2psAddPrimitiveInList Geant4_gl2psAddPrimitiveInList
#define gl2psFreeBspTree Geant4_gl2psFreeBspTree
#define gl2psGreater Geant4_gl2psGreater
#define gl2psLess Geant4_gl2psLess
#define gl2psBuildBspTree Geant4_gl2psBuildBspTree
#define gl2psTraverseBspTree Geant4_gl2psTraverseBspTree
#define gl2psRescaleAndOffset Geant4_gl2psRescaleAndOffset
#define gl2psGetPlaneFromPoints Geant4_gl2psGetPlaneFromPoints
#define gl2psFreeBspImageTree Geant4_gl2psFreeBspImageTree
#define gl2psCheckPoint Geant4_gl2psCheckPoint
#define gl2psAddPlanesInBspTreeImage Geant4_gl2psAddPlanesInBspTreeImage
#define gl2psCheckPrimitive Geant4_gl2psCheckPrimitive
#define gl2psCreateSplitPrimitive2D Geant4_gl2psCreateSplitPrimitive2D
#define gl2psSplitPrimitive2D Geant4_gl2psSplitPrimitive2D
#define gl2psAddInImageTree Geant4_gl2psAddInImageTree
#define gl2psAddInBspImageTree Geant4_gl2psAddInBspImageTree
#define gl2psAddBoundaryInList Geant4_gl2psAddBoundaryInList
#define gl2psBuildPolygonBoundary Geant4_gl2psBuildPolygonBoundary
#define gl2psAddPolyPrimitive Geant4_gl2psAddPolyPrimitive
#define gl2psGetVertex Geant4_gl2psGetVertex
#define gl2psParseFeedbackBuffer Geant4_gl2psParseFeedbackBuffer
#define gl2psWriteByte Geant4_gl2psWriteByte
#define gl2psPrintPostScriptPixmap Geant4_gl2psPrintPostScriptPixmap
#define gl2psPrintPostScriptImagemap Geant4_gl2psPrintPostScriptImagemap
#define gl2psPrintPostScriptHeader Geant4_gl2psPrintPostScriptHeader
#define gl2psPrintPostScriptColor Geant4_gl2psPrintPostScriptColor
#define gl2psResetPostScriptColor Geant4_gl2psResetPostScriptColor
#define gl2psEndPostScriptLine Geant4_gl2psEndPostScriptLine
#define gl2psParseStipplePattern Geant4_gl2psParseStipplePattern
#define gl2psPrintPostScriptDash Geant4_gl2psPrintPostScriptDash
#define gl2psPrintPostScriptPrimitive Geant4_gl2psPrintPostScriptPrimitive
#define gl2psPrintPostScriptFooter Geant4_gl2psPrintPostScriptFooter
#define gl2psPrintTeXHeader Geant4_gl2psPrintTeXHeader
#define gl2psPrintTeXPrimitive Geant4_gl2psPrintTeXPrimitive
#define gl2psPrintTeXFooter Geant4_gl2psPrintTeXFooter

#define gl2psPrintPostScriptBeginViewport Geant4_gl2psPrintPostScriptBeginViewport
#define gl2psPrintPostScriptEndViewport Geant4_gl2psPrintPostScriptEndViewport

#define gl2psPrintPostScriptFinalPrimitive Geant4_gl2psPrintPostScriptFinalPrimitive
#define gl2psPrintPrimitives Geant4_gl2psPrintPrimitives
#define gl2psPrintTeXBeginViewport Geant4_gl2psPrintTeXBeginViewport
#define gl2psPrintTeXEndViewport Geant4_gl2psPrintTeXEndViewport
#define gl2psPrintTeXFinalPrimitive Geant4_gl2psPrintTeXFinalPrimitive
#define gl2psPrintPDFCompressorType Geant4_gl2psPrintPDFCompressorType
#define gl2psPrintPDFStrokeColor Geant4_gl2psPrintPDFStrokeColor
#define gl2psPrintPDFFillColor Geant4_gl2psPrintPDFFillColor
#define gl2psPrintPDFLineWidth Geant4_gl2psPrintPDFLineWidth
#define gl2psPutPDFText Geant4_gl2psPutPDFText
#define gl2psPutPDFImage Geant4_gl2psPutPDFImage
#define gl2psPDFstacksInit Geant4_gl2psPDFstacksInit
#define gl2psPDFgroupObjectInit Geant4_gl2psPDFgroupObjectInit
#define gl2psPDFgroupListInit Geant4_gl2psPDFgroupListInit
#define gl2psSortOutTrianglePDFgroup Geant4_gl2psSortOutTrianglePDFgroup
#define gl2psPDFgroupListWriteMainStream Geant4_gl2psPDFgroupListWriteMainStream
#define gl2psPDFgroupListWriteGStateResources Geant4_gl2psPDFgroupListWriteGStateResources
#define gl2psPDFgroupListWriteShaderResources Geant4_gl2psPDFgroupListWriteShaderResources
#define gl2psPDFgroupListWriteXObjectResources Geant4_gl2psPDFgroupListWriteXObjectResources
#define gl2psPDFgroupListWriteFontResources Geant4_gl2psPDFgroupListWriteFontResources
#define gl2psPDFgroupListDelete Geant4_gl2psPDFgroupListDelete
#define gl2psPrintPDFInfo Geant4_gl2psPrintPDFInfo
#define gl2psPrintPDFCatalog Geant4_gl2psPrintPDFCatalog
#define gl2psPrintPDFPages Geant4_gl2psPrintPDFPages
#define gl2psOpenPDFDataStream Geant4_gl2psOpenPDFDataStream
#define gl2psOpenPDFDataStreamWritePreface Geant4_gl2psOpenPDFDataStreamWritePreface
#define gl2psPrintPDFHeader Geant4_gl2psPrintPDFHeader
#define gl2psPrintPDFPrimitive Geant4_gl2psPrintPDFPrimitive
#define gl2psClosePDFDataStream Geant4_gl2psClosePDFDataStream
#define gl2psPrintPDFDataStreamLength Geant4_gl2psPrintPDFDataStreamLength
#define gl2psPrintPDFOpenPage Geant4_gl2psPrintPDFOpenPage
#define gl2psPDFgroupListWriteVariableResources Geant4_gl2psPDFgroupListWriteVariableResources
#define gl2psPrintPDFGSObject Geant4_gl2psPrintPDFGSObject
#define gl2psPrintPDFShaderStreamDataCoord Geant4_gl2psPrintPDFShaderStreamDataCoord
#define gl2psPrintPDFShaderStreamDataRGB Geant4_gl2psPrintPDFShaderStreamDataRGB
#define gl2psPrintPDFShaderStreamDataAlpha Geant4_gl2psPrintPDFShaderStreamDataAlpha
#define gl2psPrintPDFShaderStreamData Geant4_gl2psPrintPDFShaderStreamData
#define gl2psPDFRectHull Geant4_gl2psPDFRectHull
#define gl2psPrintPDFShader Geant4_gl2psPrintPDFShader
#define gl2psPrintPDFShaderMask Geant4_gl2psPrintPDFShaderMask
#define gl2psPrintPDFShaderExtGS Geant4_gl2psPrintPDFShaderExtGS
#define gl2psPrintPDFShaderSimpleExtGS Geant4_gl2psPrintPDFShaderSimpleExtGS
#define gl2psPrintPDFPixmapStreamData Geant4_gl2psPrintPDFPixmapStreamData
#define gl2psPrintPDFPixmap Geant4_gl2psPrintPDFPixmap
#define gl2psPrintPDFText Geant4_gl2psPrintPDFText
#define gl2psPDFgroupListWriteObjects Geant4_gl2psPDFgroupListWriteObjects
#define gl2psPrintPDFFooter Geant4_gl2psPrintPDFFooter
#define gl2psPrintPDFBeginViewport Geant4_gl2psPrintPDFBeginViewport
#define gl2psPrintPDFEndViewport Geant4_gl2psPrintPDFEndViewport
#define gl2psPrintPDFFinalPrimitive Geant4_gl2psPrintPDFFinalPrimitive
#define gl2psSVGGetCoordsAndColors Geant4_gl2psSVGGetCoordsAndColors
#define gl2psSVGGetColorString Geant4_gl2psSVGGetColorString
#define gl2psPrintSVGHeader Geant4_gl2psPrintSVGHeader
#define gl2psPrintSVGSmoothTriangle Geant4_gl2psPrintSVGSmoothTriangle
#define gl2psPrintSVGDash Geant4_gl2psPrintSVGDash
#define gl2psEndSVGLine Geant4_gl2psEndSVGLine
#define gl2psPrintSVGPixmap Geant4_gl2psPrintSVGPixmap
#define gl2psPrintSVGPrimitive Geant4_gl2psPrintSVGPrimitive
#define gl2psPrintSVGFooter Geant4_gl2psPrintSVGFooter
#define gl2psPrintSVGBeginViewport Geant4_gl2psPrintSVGBeginViewport
#define gl2psPrintSVGEndViewport Geant4_gl2psPrintSVGEndViewport
#define gl2psPrintSVGFinalPrimitive Geant4_gl2psPrintSVGFinalPrimitive
#define gl2psPrintPGFColor Geant4_gl2psPrintPGFColor
#define gl2psPrintPGFHeader Geant4_gl2psPrintPGFHeader
#define gl2psPrintPGFDash Geant4_gl2psPrintPGFDash
#define gl2psPGFTextAlignment Geant4_gl2psPGFTextAlignment
#define gl2psPrintPGFPrimitive Geant4_gl2psPrintPGFPrimitive
#define gl2psPrintPGFFooter Geant4_gl2psPrintPGFFooter
#define gl2psPrintPGFBeginViewport Geant4_gl2psPrintPGFBeginViewport
#define gl2psPrintPGFEndViewport Geant4_gl2psPrintPGFEndViewport
#define gl2psPrintPGFFinalPrimitive Geant4_gl2psPrintPGFFinalPrimitive
#define gl2psComputeTightBoundingBox Geant4_gl2psComputeTightBoundingBox

#define gl2ps Geant4_gl2ps

#ifndef G4OPENGL_VERSION_2
#include "gl2ps.h"
#endif

#ifdef _WIN32
#define _USE_MATH_DEFINES
#endif

#endif
