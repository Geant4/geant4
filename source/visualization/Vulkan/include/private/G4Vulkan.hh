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
// Alex Vergara Gil  17th October 2026
// Ported to Vulkan from the Geant4 OpenGL visualization driver.
// G4Vulkan.hh : Vulkan rendering engine and state backend for Geant4

#ifndef G4Vulkan_hh
#define G4Vulkan_hh

#ifdef WIN32
#  include "G4VCrossPlatformMacros.hh"
#  include <windows.h>
#endif

#include <vulkan/vulkan.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Primitive types used by the G4Vulkan scene handler, viewers, and gl2ps bridge
using VKenum = uint32_t;
using VKboolean = uint8_t;
using VKbitfield = uint32_t;
using VKbyte = int8_t;
using VKshort = int16_t;
using VKint = int32_t;
using VKsizei = int32_t;
using VKubyte = uint8_t;
using VKushort = uint16_t;
using VKuint = uint32_t;
using VKfloat = float;
using VKclampf = float;
using VKdouble = double;
using VKclampd = double;
using VKvoid = void;

// Boolean values
#ifndef VK_FALSE
#  define VK_FALSE 0
#  define VK_TRUE 1
#endif

// Primitive topology constants
#ifndef VK_POINTS
#  define VK_POINTS 0x0000
#  define VK_LINES 0x0001
#  define VK_LINE_LOOP 0x0002
#  define VK_LINE_STRIP 0x0003
#  define VK_TRIANGLES 0x0004
#  define VK_TRIANGLE_STRIP 0x0005
#  define VK_TRIANGLE_FAN 0x0006
#  define VK_QUADS 0x0007
#  define VK_QUAD_STRIP 0x0008
#  define VK_POLYGON 0x0009
#endif

// Buffer bits
#ifndef VK_DEPTH_BUFFER_BIT
#  define VK_DEPTH_BUFFER_BIT 0x00000100
#  define VK_STENCIL_BUFFER_BIT 0x00000400
#  define VK_COLOR_BUFFER_BIT 0x00004000
#  define VK_LIST_BIT 0x00020000
#endif

// Matrix modes
#ifndef VK_MODELVIEW
#  define VK_MODELVIEW 0x1700
#  define VK_PROJECTION 0x1701
#  define VK_TEXTURE 0x1702
#  define VK_MODELVIEW_MATRIX 0x0BA6
#  define VK_PROJECTION_MATRIX 0x0BA7
#  define VK_VIEWPORT 0x0BA2
#endif

// Depth / Alpha functions
#ifndef VK_NEVER
#  define VK_NEVER 0x0200
#  define VK_LESS 0x0201
#  define VK_EQUAL 0x0202
#  define VK_LEQUAL 0x0203
#  define VK_GREATER 0x0204
#  define VK_NOTEQUAL 0x0205
#  define VK_GEQUAL 0x0206
#  define VK_ALWAYS 0x0207
#endif

// Blending
#ifndef VK_ZERO
#  define VK_ZERO 0
#  define VK_ONE 1
#  define VK_SRC_COLOR 0x0300
#  define VK_ONE_MINUS_SRC_COLOR 0x0301
#  define VK_SRC_ALPHA 0x0302
#  define VK_ONE_MINUS_SRC_ALPHA 0x0303
#  define VK_DST_ALPHA 0x0304
#  define VK_ONE_MINUS_DST_ALPHA 0x0305
#  define VK_DST_COLOR 0x0306
#  define VK_ONE_MINUS_DST_COLOR 0x0307
#  define VK_SRC_ALPHA_SATURATE 0x0308
#  define VK_BLEND 0x0BE2
#  define VK_BLEND_SRC 0x0BE1
#  define VK_BLEND_DST 0x0BE0
#endif

// Render modes
#ifndef VK_RENDER
#  define VK_RENDER 0x1C00
#  define VK_FEEDBACK 0x1C01
#  define VK_SELECT 0x1C02
#  define VK_RENDER_MODE 0x0C40
#endif

// Feedback tokens & types
#ifndef VK_2D
#  define VK_2D 0x0600
#  define VK_3D 0x0601
#  define VK_3D_COLOR 0x0602
#  define VK_3D_COLOR_TEXTURE 0x0603
#  define VK_4D_COLOR_TEXTURE 0x0604
#  define VK_PASS_THROUGH_TOKEN 0x0700
#  define VK_POINT_TOKEN 0x0701
#  define VK_LINE_TOKEN 0x0702
#  define VK_POLYGON_TOKEN 0x0703
#  define VK_BITMAP_TOKEN 0x0704
#  define VK_DRAW_PIXEL_TOKEN 0x0705
#  define VK_COPY_PIXEL_TOKEN 0x0706
#  define VK_LINE_RESET_TOKEN 0x0707
#endif

// State enables & queries
#ifndef VK_POINT_SMOOTH
#  define VK_POINT_SMOOTH 0x0B10
#  define VK_POINT_SIZE 0x0B11
#  define VK_LINE_SMOOTH 0x0B20
#  define VK_LINE_WIDTH 0x0B21
#  define VK_LINE_STIPPLE 0x0B24
#  define VK_LINE_STIPPLE_PATTERN 0x0B25
#  define VK_LINE_STIPPLE_REPEAT 0x0B26
#  define VK_POLYGON_MODE 0x0B40
#  define VK_POLYGON_SMOOTH 0x0B41
#  define VK_POLYGON_STIPPLE 0x0B42
#  define VK_EDGE_FLAG 0x0B43
#  define VK_CULL_FACE 0x0B44
#  define VK_CULL_FACE_MODE 0x0B45
#  define VK_FRONT_FACE 0x0B46
#  define VK_LIGHTING 0x0B50
#  define VK_LIGHT_MODEL_TWO_SIDE 0x0B52
#  define VK_LIGHT_MODEL_AMBIENT 0x0B53
#  define VK_SHADE_MODEL 0x0B54
#  define VK_COLOR_MATERIAL 0x0B57
#  define VK_DEPTH_TEST 0x0B71
#  define VK_DEPTH_WRITEMASK 0x0B72
#  define VK_DEPTH_CLEAR_VALUE 0x0B73
#  define VK_DEPTH_FUNC 0x0B74
#  define VK_STENCIL_TEST 0x0B90
#  define VK_NORMALIZE 0x0BA1
#  define VK_DITHER 0x0BD0
#  define VK_SCISSOR_TEST 0x0C11
#  define VK_CURRENT_COLOR 0x0B00
#  define VK_COLOR_CLEAR_VALUE 0x0C22
#  define VK_MAX_CLIP_PLANES 0x0D32
#  define VK_MAX_LIGHTS 0x0D31
#  define VK_MAX_VIEWPORT_DIMS 0x0D3A
#  define VK_POLYGON_OFFSET_FACTOR 0x8038
#  define VK_POLYGON_OFFSET_UNITS 0x2A00
#  define VK_POLYGON_OFFSET_POINT 0x2A01
#  define VK_POLYGON_OFFSET_LINE 0x2A02
#  define VK_POLYGON_OFFSET_FILL 0x8037
#endif

// Polygon / Face / Shading / Material / Light constants
#ifndef VK_POINT
#  define VK_POINT 0x1B00
#  define VK_LINE 0x1B01
#  define VK_FILL 0x1B02
#  define VK_CW 0x0900
#  define VK_CCW 0x0901
#  define VK_FRONT 0x0404
#  define VK_BACK 0x0405
#  define VK_FRONT_AND_BACK 0x0408
#  define VK_FLAT 0x1D00
#  define VK_SMOOTH 0x1D01
#  define VK_AMBIENT 0x1200
#  define VK_DIFFUSE 0x1201
#  define VK_SPECULAR 0x1202
#  define VK_POSITION 0x1203
#  define VK_SPOT_DIRECTION 0x1204
#  define VK_SPOT_EXPONENT 0x1205
#  define VK_SPOT_CUTOFF 0x1206
#  define VK_CONSTANT_ATTENUATION 0x1207
#  define VK_LINEAR_ATTENUATION 0x1208
#  define VK_QUADRATIC_ATTENUATION 0x1209
#  define VK_EMISSION 0x1600
#  define VK_SHININESS 0x1601
#  define VK_AMBIENT_AND_DIFFUSE 0x1602
#  define VK_LIGHT0 0x4000
#  define VK_LIGHT1 0x4001
#  define VK_LIGHT2 0x4002
#  define VK_LIGHT3 0x4003
#  define VK_LIGHT4 0x4004
#  define VK_LIGHT5 0x4005
#  define VK_LIGHT6 0x4006
#  define VK_LIGHT7 0x4007
#endif

// Clip planes
#ifndef VK_CLIP_PLANE0
#  define VK_CLIP_PLANE0 0x3000
#  define VK_CLIP_PLANE1 0x3001
#  define VK_CLIP_PLANE2 0x3002
#  define VK_CLIP_PLANE3 0x3003
#  define VK_CLIP_PLANE4 0x3004
#  define VK_CLIP_PLANE5 0x3005
#endif

// Display list modes
#ifndef VK_COMPILE
#  define VK_COMPILE 0x1300
#  define VK_COMPILE_AND_EXECUTE 0x1301
#endif

// Data types & Pixel formats
#ifndef VK_BYTE
#  define VK_BYTE 0x1400
#  define VK_UNSIGNED_BYTE 0x1401
#  define VK_SHORT 0x1402
#  define VK_UNSIGNED_SHORT 0x1403
#  define VK_INT 0x1404
#  define VK_UNSIGNED_INT 0x1405
#  define VK_FLOAT 0x1406
#  define VK_2_BYTES 0x1407
#  define VK_3_BYTES 0x1408
#  define VK_4_BYTES 0x1409
#  define VK_DOUBLE 0x140A
#  define VK_RGB 0x1907
#  define VK_RGBA 0x1908
#  define VK_PACK_ALIGNMENT 0x0D05
#  define VK_PACK_ALIGNMENT 0x0D05
#  define VK_UNPACK_ALIGNMENT 0x0CF5
#endif

// Hints
#ifndef VK_DONT_CARE
#  define VK_DONT_CARE 0x1100
#  define VK_FASTEST 0x1101
#  define VK_NICEST 0x1102
#  define VK_PERSPECTIVE_CORRECTION_HINT 0x0C50
#  define VK_POINT_SMOOTH_HINT 0x0C51
#  define VK_LINE_SMOOTH_HINT 0x0C52
#  define VK_POLYGON_SMOOTH_HINT 0x0C53
#  define VK_FOG_HINT 0x0C54
#endif

// Error codes
#ifndef VK_NO_ERROR
#  define VK_NO_ERROR 0
#  define VK_INVALID_ENUM 0x0500
#  define VK_INVALID_VALUE 0x0501
#  define VK_INVALID_OPERATION 0x0502
#  define VK_STACK_OVERFLOW 0x0503
#  define VK_STACK_UNDERFLOW 0x0504
#  define VK_OUT_OF_MEMORY 0x0505
#endif

// Precompiled SPIR-V Vertex Shader (compiled from GLSL 450 with glslc)
inline constexpr uint32_t kG4VulkanVertSPV[] = {
  0x07230203, 0x00010000, 0x000d000b, 0x00000056, 0x00000000, 0x00020011, 0x00000001, 0x0006000b,
  0x00000001, 0x4c534c47, 0x6474732e, 0x3035342e, 0x00000000, 0x0003000e, 0x00000000, 0x00000001,
  0x000a000f, 0x00000000, 0x00000004, 0x6e69616d, 0x00000000, 0x0000000d, 0x00000019, 0x00000034,
  0x00000045, 0x00000047, 0x00030003, 0x00000002, 0x000001c2, 0x000a0004, 0x475f4c47, 0x4c474f4f,
  0x70635f45, 0x74735f70, 0x5f656c79, 0x656e696c, 0x7269645f, 0x69746365, 0x00006576, 0x00080004,
  0x475f4c47, 0x4c474f4f, 0x6e695f45, 0x64756c63, 0x69645f65, 0x74636572, 0x00657669, 0x00040005,
  0x00000004, 0x6e69616d, 0x00000000, 0x00060005, 0x0000000b, 0x505f6c67, 0x65567265, 0x78657472,
  0x00000000, 0x00060006, 0x0000000b, 0x00000000, 0x505f6c67, 0x7469736f, 0x006e6f69, 0x00070006,
  0x0000000b, 0x00000001, 0x505f6c67, 0x746e696f, 0x657a6953, 0x00000000, 0x00070006, 0x0000000b,
  0x00000002, 0x435f6c67, 0x4470696c, 0x61747369, 0x0065636e, 0x00070006, 0x0000000b, 0x00000003,
  0x435f6c67, 0x446c6c75, 0x61747369, 0x0065636e, 0x00030005, 0x0000000d, 0x00000000, 0x00060005,
  0x00000011, 0x68737550, 0x736e6f43, 0x746e6174, 0x00000073, 0x00040006, 0x00000011, 0x00000000,
  0x0070766d, 0x00060006, 0x00000011, 0x00000001, 0x6867696c, 0x72694474, 0x00000000, 0x00060006,
  0x00000011, 0x00000002, 0x6e696f70, 0x7a695374, 0x00000065, 0x00060006, 0x00000011, 0x00000003,
  0x6867696c, 0x676e6974, 0x00000000, 0x00030005, 0x00000013, 0x00006370, 0x00040005, 0x00000019,
  0x6f506e69, 0x00000073, 0x00030005, 0x00000033, 0x0000006e, 0x00050005, 0x00000034, 0x6f4e6e69,
  0x6c616d72, 0x00000000, 0x00030005, 0x00000037, 0x0000006c, 0x00040005, 0x0000003e, 0x66666964,
  0x00000000, 0x00050005, 0x00000045, 0x67617266, 0x6f6c6f43, 0x00000072, 0x00040005, 0x00000047,
  0x6f436e69, 0x00726f6c, 0x00030047, 0x0000000b, 0x00000002, 0x00050048, 0x0000000b, 0x00000000,
  0x0000000b, 0x00000000, 0x00050048, 0x0000000b, 0x00000001, 0x0000000b, 0x00000001, 0x00050048,
  0x0000000b, 0x00000002, 0x0000000b, 0x00000003, 0x00050048, 0x0000000b, 0x00000003, 0x0000000b,
  0x00000004, 0x00030047, 0x00000011, 0x00000002, 0x00040048, 0x00000011, 0x00000000, 0x00000005,
  0x00050048, 0x00000011, 0x00000000, 0x00000007, 0x00000010, 0x00050048, 0x00000011, 0x00000000,
  0x00000023, 0x00000000, 0x00050048, 0x00000011, 0x00000001, 0x00000023, 0x00000040, 0x00050048,
  0x00000011, 0x00000002, 0x00000023, 0x00000050, 0x00050048, 0x00000011, 0x00000003, 0x00000023,
  0x00000054, 0x00040047, 0x00000019, 0x0000001e, 0x00000000, 0x00040047, 0x00000034, 0x0000001e,
  0x00000001, 0x00040047, 0x00000045, 0x0000001e, 0x00000000, 0x00040047, 0x00000047, 0x0000001e,
  0x00000002, 0x00020013, 0x00000002, 0x00030021, 0x00000003, 0x00000002, 0x00030016, 0x00000006,
  0x00000020, 0x00040017, 0x00000007, 0x00000006, 0x00000004, 0x00040015, 0x00000008, 0x00000020,
  0x00000000, 0x0004002b, 0x00000008, 0x00000009, 0x00000001, 0x0004001c, 0x0000000a, 0x00000006,
  0x00000009, 0x0006001e, 0x0000000b, 0x00000007, 0x00000006, 0x0000000a, 0x0000000a, 0x00040020,
  0x0000000c, 0x00000003, 0x0000000b, 0x0004003b, 0x0000000c, 0x0000000d, 0x00000003, 0x00040015,
  0x0000000e, 0x00000020, 0x00000001, 0x0004002b, 0x0000000e, 0x0000000f, 0x00000000, 0x00040018,
  0x00000010, 0x00000007, 0x00000004, 0x0006001e, 0x00000011, 0x00000010, 0x00000007, 0x00000006,
  0x0000000e, 0x00040020, 0x00000012, 0x00000009, 0x00000011, 0x0004003b, 0x00000012, 0x00000013,
  0x00000009, 0x00040020, 0x00000014, 0x00000009, 0x00000010, 0x00040017, 0x00000017, 0x00000006,
  0x00000003, 0x00040020, 0x00000018, 0x00000001, 0x00000017, 0x0004003b, 0x00000018, 0x00000019,
  0x00000001, 0x0004002b, 0x00000006, 0x0000001b, 0x3f800000, 0x00040020, 0x00000021, 0x00000003,
  0x00000007, 0x0004002b, 0x0000000e, 0x00000023, 0x00000001, 0x0004002b, 0x0000000e, 0x00000024,
  0x00000002, 0x00040020, 0x00000025, 0x00000009, 0x00000006, 0x00040020, 0x00000028, 0x00000003,
  0x00000006, 0x0004002b, 0x0000000e, 0x0000002a, 0x00000003, 0x00040020, 0x0000002b, 0x00000009,
  0x0000000e, 0x00020014, 0x0000002e, 0x00040020, 0x00000032, 0x00000007, 0x00000017, 0x0004003b,
  0x00000018, 0x00000034, 0x00000001, 0x00040020, 0x00000038, 0x00000009, 0x00000007, 0x00040020,
  0x0000003d, 0x00000007, 0x00000006, 0x0004002b, 0x00000006, 0x00000043, 0x3e800000, 0x0004003b,
  0x00000021, 0x00000045, 0x00000003, 0x00040020, 0x00000046, 0x00000001, 0x00000007, 0x0004003b,
  0x00000046, 0x00000047, 0x00000001, 0x0004002b, 0x00000008, 0x0000004c, 0x00000003, 0x00040020,
  0x0000004d, 0x00000001, 0x00000006, 0x00050036, 0x00000002, 0x00000004, 0x00000000, 0x00000003,
  0x000200f8, 0x00000005, 0x0004003b, 0x00000032, 0x00000033, 0x00000007, 0x0004003b, 0x00000032,
  0x00000037, 0x00000007, 0x0004003b, 0x0000003d, 0x0000003e, 0x00000007, 0x00050041, 0x00000014,
  0x00000015, 0x00000013, 0x0000000f, 0x0004003d, 0x00000010, 0x00000016, 0x00000015, 0x0004003d,
  0x00000017, 0x0000001a, 0x00000019, 0x00050051, 0x00000006, 0x0000001c, 0x0000001a, 0x00000000,
  0x00050051, 0x00000006, 0x0000001d, 0x0000001a, 0x00000001, 0x00050051, 0x00000006, 0x0000001e,
  0x0000001a, 0x00000002, 0x00070050, 0x00000007, 0x0000001f, 0x0000001c, 0x0000001d, 0x0000001e,
  0x0000001b, 0x00050091, 0x00000007, 0x00000020, 0x00000016, 0x0000001f, 0x00050041, 0x00000021,
  0x00000022, 0x0000000d, 0x0000000f, 0x0003003e, 0x00000022, 0x00000020, 0x00050041, 0x00000025,
  0x00000026, 0x00000013, 0x00000024, 0x0004003d, 0x00000006, 0x00000027, 0x00000026, 0x00050041,
  0x00000028, 0x00000029, 0x0000000d, 0x00000023, 0x0003003e, 0x00000029, 0x00000027, 0x00050041,
  0x0000002b, 0x0000002c, 0x00000013, 0x0000002a, 0x0004003d, 0x0000000e, 0x0000002d, 0x0000002c,
  0x000500ab, 0x0000002e, 0x0000002f, 0x0000002d, 0x0000000f, 0x000300f7, 0x00000031, 0x00000000,
  0x000400fa, 0x0000002f, 0x00000030, 0x00000054, 0x000200f8, 0x00000030, 0x0004003d, 0x00000017,
  0x00000035, 0x00000034, 0x0006000c, 0x00000017, 0x00000036, 0x00000001, 0x00000045, 0x00000035,
  0x0003003e, 0x00000033, 0x00000036, 0x00050041, 0x00000038, 0x00000039, 0x00000013, 0x00000023,
  0x0004003d, 0x00000007, 0x0000003a, 0x00000039, 0x0008004f, 0x00000017, 0x0000003b, 0x0000003a,
  0x0000003a, 0x00000000, 0x00000001, 0x00000002, 0x0006000c, 0x00000017, 0x0000003c, 0x00000001,
  0x00000045, 0x0000003b, 0x0003003e, 0x00000037, 0x0000003c, 0x0004003d, 0x00000017, 0x0000003f,
  0x00000033, 0x0004003d, 0x00000017, 0x00000040, 0x00000037, 0x00050094, 0x00000006, 0x00000041,
  0x0000003f, 0x00000040, 0x0006000c, 0x00000006, 0x00000042, 0x00000001, 0x00000004, 0x00000041,
  0x0007000c, 0x00000006, 0x00000044, 0x00000001, 0x00000028, 0x00000042, 0x00000043, 0x0003003e,
  0x0000003e, 0x00000044, 0x0004003d, 0x00000007, 0x00000048, 0x00000047, 0x0008004f, 0x00000017,
  0x00000049, 0x00000048, 0x00000048, 0x00000000, 0x00000001, 0x00000002, 0x0004003d, 0x00000006,
  0x0000004a, 0x0000003e, 0x0005008e, 0x00000017, 0x0000004b, 0x00000049, 0x0000004a, 0x00050041,
  0x0000004d, 0x0000004e, 0x00000047, 0x0000004c, 0x0004003d, 0x00000006, 0x0000004f, 0x0000004e,
  0x00050051, 0x00000006, 0x00000050, 0x0000004b, 0x00000000, 0x00050051, 0x00000006, 0x00000051,
  0x0000004b, 0x00000001, 0x00050051, 0x00000006, 0x00000052, 0x0000004b, 0x00000002, 0x00070050,
  0x00000007, 0x00000053, 0x00000050, 0x00000051, 0x00000052, 0x0000004f, 0x0003003e, 0x00000045,
  0x00000053, 0x000200f9, 0x00000031, 0x000200f8, 0x00000054, 0x0004003d, 0x00000007, 0x00000055,
  0x00000047, 0x0003003e, 0x00000045, 0x00000055, 0x000200f9, 0x00000031, 0x000200f8, 0x00000031,
  0x000100fd, 0x00010038};

// Precompiled SPIR-V Fragment Shader (compiled from GLSL 450 with glslc)
inline constexpr uint32_t kG4VulkanFragSPV[] = {
  0x07230203, 0x00010000, 0x000d000b, 0x0000000d, 0x00000000, 0x00020011, 0x00000001, 0x0006000b,
  0x00000001, 0x4c534c47, 0x6474732e, 0x3035342e, 0x00000000, 0x0003000e, 0x00000000, 0x00000001,
  0x0007000f, 0x00000004, 0x00000004, 0x6e69616d, 0x00000000, 0x00000009, 0x0000000b, 0x00030010,
  0x00000004, 0x00000007, 0x00030003, 0x00000002, 0x000001c2, 0x000a0004, 0x475f4c47, 0x4c474f4f,
  0x70635f45, 0x74735f70, 0x5f656c79, 0x656e696c, 0x7269645f, 0x69746365, 0x00006576, 0x00080004,
  0x475f4c47, 0x4c474f4f, 0x6e695f45, 0x64756c63, 0x69645f65, 0x74636572, 0x00657669, 0x00040005,
  0x00000004, 0x6e69616d, 0x00000000, 0x00050005, 0x00000009, 0x4374756f, 0x726f6c6f, 0x00000000,
  0x00050005, 0x0000000b, 0x67617266, 0x6f6c6f43, 0x00000072, 0x00040047, 0x00000009, 0x0000001e,
  0x00000000, 0x00040047, 0x0000000b, 0x0000001e, 0x00000000, 0x00020013, 0x00000002, 0x00030021,
  0x00000003, 0x00000002, 0x00030016, 0x00000006, 0x00000020, 0x00040017, 0x00000007, 0x00000006,
  0x00000004, 0x00040020, 0x00000008, 0x00000003, 0x00000007, 0x0004003b, 0x00000008, 0x00000009,
  0x00000003, 0x00040020, 0x0000000a, 0x00000001, 0x00000007, 0x0004003b, 0x0000000a, 0x0000000b,
  0x00000001, 0x00050036, 0x00000002, 0x00000004, 0x00000000, 0x00000003, 0x000200f8, 0x00000005,
  0x0004003d, 0x00000007, 0x0000000c, 0x0000000b, 0x0003003e, 0x00000009, 0x0000000c, 0x000100fd,
  0x00010038};

struct G4VulkanVertex
{
  float pos[3];
  float normal[3];
  float color[4];
};

struct G4VulkanPushConstants
{
  float mvp[16];
  float lightDir[4];
  float pointSize;
  int32_t lighting;
};

struct G4VulkanDrawCmd
{
  VkPipeline pipeline{VK_NULL_HANDLE};
  uint32_t firstVertex{0};
  uint32_t vertexCount{0};
  G4VulkanPushConstants pc{};
};

struct G4VulkanDisplayList
{
  std::vector<G4VulkanVertex> triangles;
  std::vector<G4VulkanVertex> lines;
  std::vector<G4VulkanVertex> points;
  uint32_t gpuTriFirst{0};
  uint32_t gpuTriCount{0};
  uint32_t gpuLineFirst{0};
  uint32_t gpuLineCount{0};
  uint32_t gpuPointFirst{0};
  uint32_t gpuPointCount{0};
  float minPos[3]{1e30f, 1e30f, 1e30f};
  float maxPos[3]{-1e30f, -1e30f, -1e30f};
  bool hasBounds{false};
  bool hasLightingState{false};
  bool lighting{true};
  std::vector<std::function<void()>> stateOps;

  void Clear()
  {
    triangles.clear();
    lines.clear();
    points.clear();
    gpuTriFirst = gpuTriCount = 0;
    gpuLineFirst = gpuLineCount = 0;
    gpuPointFirst = gpuPointCount = 0;
    minPos[0] = minPos[1] = minPos[2] = 1e30f;
    maxPos[0] = maxPos[1] = maxPos[2] = -1e30f;
    hasBounds = false;
    hasLightingState = false;
    lighting = true;
    stateOps.clear();
  }

  void AddBounds(const G4VulkanVertex& v)
  {
    hasBounds = true;
    for (int i = 0; i < 3; ++i)
    {
      if (v.pos[i] < minPos[i]) minPos[i] = v.pos[i];
      if (v.pos[i] > maxPos[i]) maxPos[i] = v.pos[i];
    }
  }
};

class G4VulkanContext
{
  public:
    using Matrix4d = std::array<double, 16>;
    static constexpr uint32_t kMaxStaticVertices = 524288;
    static constexpr uint32_t kMaxDynamicVertices = 524288;
    static constexpr uint32_t kCmdRingSize = 64;
    static constexpr VkDeviceSize kStagingBufferSize = 16 * 1024 * 1024;

    static G4VulkanContext& Instance()
    {
      static G4VulkanContext instance;
      return instance;
    }

    G4VulkanContext() { InitState(); }

    ~G4VulkanContext() { DestroyVulkan(); }

    void InitState()
    {
      fMatrixMode = VK_MODELVIEW;
      fModelViewStack.clear();
      fProjectionStack.clear();
      fTextureStack.clear();
      fModelViewStack.push_back(IdentityMatrix());
      fProjectionStack.push_back(IdentityMatrix());
      fTextureStack.push_back(IdentityMatrix());
      fCachedMvpDirty = true;
      ResizeFramebuffer(fViewport[2], fViewport[3]);
    }

    static Matrix4d IdentityMatrix()
    {
      Matrix4d m{};
      m[0] = m[5] = m[10] = m[15] = 1.0;
      return m;
    }

    static Matrix4d MultiplyMatrix(const Matrix4d& a, const Matrix4d& b)
    {
      Matrix4d r{};
      for (int c = 0; c < 4; ++c)
      {
        for (int row = 0; row < 4; ++row)
        {
          double sum = 0.0;
          for (int k = 0; k < 4; ++k)
          {
            sum += a[k * 4 + row] * b[c * 4 + k];
          }
          r[c * 4 + row] = sum;
        }
      }
      return r;
    }

    static void TransformVec4(const Matrix4d& m, const double in[4], double out[4])
    {
      for (int row = 0; row < 4; ++row)
      {
        out[row] =
          m[0 * 4 + row] * in[0] + m[1 * 4 + row] * in[1] + m[2 * 4 + row] * in[2] + m[3 * 4 + row] * in[3];
      }
    }

    Matrix4d& CurrentMatrix()
    {
      fCachedMvpDirty = true;
      if (fMatrixMode == VK_PROJECTION) return fProjectionStack.back();
      if (fMatrixMode == VK_TEXTURE) return fTextureStack.back();
      return fModelViewStack.back();
    }

    std::vector<Matrix4d>& CurrentStack()
    {
      fCachedMvpDirty = true;
      if (fMatrixMode == VK_PROJECTION) return fProjectionStack;
      if (fMatrixMode == VK_TEXTURE) return fTextureStack;
      return fModelViewStack;
    }

    uint32_t FindMemoryType(uint32_t typeFilter, VkMemoryPropertyFlags properties) const
    {
      for (uint32_t i = 0; i < fMemProperties.memoryTypeCount; ++i)
      {
        if ((typeFilter & (1u << i)) &&
            (fMemProperties.memoryTypes[i].propertyFlags & properties) == properties)
        {
          return i;
        }
      }
      return UINT32_MAX;
    }

    void EnsureVulkanInitialised()
    {
      if (fVulkanInitAttempted) return;
      fVulkanInitAttempted = true;

      VkApplicationInfo appInfo{};
      appInfo.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO;
      appInfo.pApplicationName = "Geant4";
      appInfo.applicationVersion = VK_MAKE_VERSION(11, 5, 0);
      appInfo.pEngineName = "G4Vulkan";
      appInfo.engineVersion = VK_MAKE_VERSION(1, 0, 0);
      appInfo.apiVersion = VK_API_VERSION_1_0;

      VkInstanceCreateInfo createInfo{};
      createInfo.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
      createInfo.pApplicationInfo = &appInfo;

#ifdef VK_KHR_portability_enumeration
      const char* extNames[] = {VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME};
      createInfo.flags |= VK_INSTANCE_CREATE_ENUMERATE_PORTABILITY_BIT_KHR;
      createInfo.enabledExtensionCount = 1;
      createInfo.ppEnabledExtensionNames = extNames;
#endif

      if (vkCreateInstance(&createInfo, nullptr, &fInstance) != VK_SUCCESS)
      {
        createInfo.flags = 0;
        createInfo.enabledExtensionCount = 0;
        createInfo.ppEnabledExtensionNames = nullptr;
        if (vkCreateInstance(&createInfo, nullptr, &fInstance) != VK_SUCCESS)
        {
          fInstance = VK_NULL_HANDLE;
          return;
        }
      }

      uint32_t deviceCount = 0;
      if (vkEnumeratePhysicalDevices(fInstance, &deviceCount, nullptr) != VK_SUCCESS || deviceCount == 0)
      {
        return;
      }
      std::vector<VkPhysicalDevice> devices(deviceCount);
      vkEnumeratePhysicalDevices(fInstance, &deviceCount, devices.data());
      fPhysicalDevice = devices[0];
      vkGetPhysicalDeviceMemoryProperties(fPhysicalDevice, &fMemProperties);

      uint32_t queueFamilyCount = 0;
      vkGetPhysicalDeviceQueueFamilyProperties(fPhysicalDevice, &queueFamilyCount, nullptr);
      std::vector<VkQueueFamilyProperties> queueFamilies(queueFamilyCount);
      vkGetPhysicalDeviceQueueFamilyProperties(fPhysicalDevice, &queueFamilyCount, queueFamilies.data());

      uint32_t graphicsQueueFamily = UINT32_MAX;
      for (uint32_t i = 0; i < queueFamilyCount; ++i)
      {
        if (queueFamilies[i].queueFlags & VK_QUEUE_GRAPHICS_BIT)
        {
          graphicsQueueFamily = i;
          break;
        }
      }
      if (graphicsQueueFamily == UINT32_MAX) return;
      fGraphicsQueueFamily = graphicsQueueFamily;

      float queuePriority = 1.0f;
      VkDeviceQueueCreateInfo queueCreateInfo{};
      queueCreateInfo.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
      queueCreateInfo.queueFamilyIndex = fGraphicsQueueFamily;
      queueCreateInfo.queueCount = 1;
      queueCreateInfo.pQueuePriorities = &queuePriority;

      VkPhysicalDeviceFeatures deviceFeatures{};
      VkDeviceCreateInfo deviceCreateInfo{};
      deviceCreateInfo.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
      deviceCreateInfo.pQueueCreateInfos = &queueCreateInfo;
      deviceCreateInfo.queueCreateInfoCount = 1;
      deviceCreateInfo.pEnabledFeatures = &deviceFeatures;

      if (vkCreateDevice(fPhysicalDevice, &deviceCreateInfo, nullptr, &fDevice) != VK_SUCCESS)
      {
        fDevice = VK_NULL_HANDLE;
        return;
      }

      vkGetDeviceQueue(fDevice, fGraphicsQueueFamily, 0, &fGraphicsQueue);

      VkCommandPoolCreateInfo poolInfo{};
      poolInfo.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
      poolInfo.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
      poolInfo.queueFamilyIndex = fGraphicsQueueFamily;
      if (vkCreateCommandPool(fDevice, &poolInfo, nullptr, &fCommandPool) != VK_SUCCESS)
      {
        return;
      }

      VkCommandBufferAllocateInfo allocInfo{};
      allocInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
      allocInfo.commandPool = fCommandPool;
      allocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
      allocInfo.commandBufferCount = 1;
      vkAllocateCommandBuffers(fDevice, &allocInfo, &fCommandBuffer);

      allocInfo.commandBufferCount = kCmdRingSize;
      vkAllocateCommandBuffers(fDevice, &allocInfo, fCmdRing.data());

      VkFenceCreateInfo fenceInfo{};
      fenceInfo.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
      for (uint32_t i = 0; i < kCmdRingSize; ++i)
      {
        vkCreateFence(fDevice, &fenceInfo, nullptr, &fFenceRing[i]);
        fFenceSubmitted[i] = false;
      }
      vkCreateFence(fDevice, &fenceInfo, nullptr, &fTransferFence);

      CreateVulkanPipelines();
      CreateVulkanBuffers();
    }

    void CreateVulkanPipelines()
    {
      if (fDevice == VK_NULL_HANDLE) return;

      VkAttachmentDescription colorAttachment{};
      colorAttachment.format = VK_FORMAT_R8G8B8A8_UNORM;
      colorAttachment.samples = VK_SAMPLE_COUNT_1_BIT;
      colorAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
      colorAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
      colorAttachment.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
      colorAttachment.stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
      colorAttachment.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
      colorAttachment.finalLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;

      VkAttachmentDescription depthAttachment{};
      depthAttachment.format = VK_FORMAT_D32_SFLOAT;
      depthAttachment.samples = VK_SAMPLE_COUNT_1_BIT;
      depthAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
      depthAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
      depthAttachment.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
      depthAttachment.stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
      depthAttachment.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
      depthAttachment.finalLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;

      VkAttachmentReference colorRef{0, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL};
      VkAttachmentReference depthRef{1, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL};

      VkSubpassDescription subpass{};
      subpass.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
      subpass.colorAttachmentCount = 1;
      subpass.pColorAttachments = &colorRef;
      subpass.pDepthStencilAttachment = &depthRef;

      std::array<VkAttachmentDescription, 2> attachments = {colorAttachment, depthAttachment};
      VkRenderPassCreateInfo renderPassInfo{};
      renderPassInfo.sType = VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO;
      renderPassInfo.attachmentCount = static_cast<uint32_t>(attachments.size());
      renderPassInfo.pAttachments = attachments.data();
      renderPassInfo.subpassCount = 1;
      renderPassInfo.pSubpasses = &subpass;

      if (vkCreateRenderPass(fDevice, &renderPassInfo, nullptr, &fRenderPass) != VK_SUCCESS)
      {
        return;
      }

      attachments[0].loadOp = VK_ATTACHMENT_LOAD_OP_LOAD;
      attachments[0].initialLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;
      attachments[1].loadOp = VK_ATTACHMENT_LOAD_OP_LOAD;
      attachments[1].initialLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;
      vkCreateRenderPass(fDevice, &renderPassInfo, nullptr, &fRenderPassLoad);

      VkShaderModuleCreateInfo vertInfo{};
      vertInfo.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
      vertInfo.codeSize = sizeof(kG4VulkanVertSPV);
      vertInfo.pCode = kG4VulkanVertSPV;
      vkCreateShaderModule(fDevice, &vertInfo, nullptr, &fVertShader);

      VkShaderModuleCreateInfo fragInfo{};
      fragInfo.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
      fragInfo.codeSize = sizeof(kG4VulkanFragSPV);
      fragInfo.pCode = kG4VulkanFragSPV;
      vkCreateShaderModule(fDevice, &fragInfo, nullptr, &fFragShader);

      VkPushConstantRange pushRange{};
      pushRange.stageFlags = VK_SHADER_STAGE_VERTEX_BIT;
      pushRange.offset = 0;
      pushRange.size = sizeof(G4VulkanPushConstants);

      VkPipelineLayoutCreateInfo layoutInfo{};
      layoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
      layoutInfo.pushConstantRangeCount = 1;
      layoutInfo.pPushConstantRanges = &pushRange;
      vkCreatePipelineLayout(fDevice, &layoutInfo, nullptr, &fPipelineLayout);

      VkPipelineShaderStageCreateInfo shaderStages[2]{};
      shaderStages[0].sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
      shaderStages[0].stage = VK_SHADER_STAGE_VERTEX_BIT;
      shaderStages[0].module = fVertShader;
      shaderStages[0].pName = "main";
      shaderStages[1].sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
      shaderStages[1].stage = VK_SHADER_STAGE_FRAGMENT_BIT;
      shaderStages[1].module = fFragShader;
      shaderStages[1].pName = "main";

      VkVertexInputBindingDescription bindingDesc{};
      bindingDesc.binding = 0;
      bindingDesc.stride = sizeof(G4VulkanVertex);
      bindingDesc.inputRate = VK_VERTEX_INPUT_RATE_VERTEX;

      std::array<VkVertexInputAttributeDescription, 3> attrDescs{};
      attrDescs[0] = {0, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(G4VulkanVertex, pos)};
      attrDescs[1] = {1, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(G4VulkanVertex, normal)};
      attrDescs[2] = {2, 0, VK_FORMAT_R32G32B32A32_SFLOAT, offsetof(G4VulkanVertex, color)};

      VkPipelineVertexInputStateCreateInfo vertexInput{};
      vertexInput.sType = VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO;
      vertexInput.vertexBindingDescriptionCount = 1;
      vertexInput.pVertexBindingDescriptions = &bindingDesc;
      vertexInput.vertexAttributeDescriptionCount = static_cast<uint32_t>(attrDescs.size());
      vertexInput.pVertexAttributeDescriptions = attrDescs.data();

      VkPipelineViewportStateCreateInfo viewportState{};
      viewportState.sType = VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO;
      viewportState.viewportCount = 1;
      viewportState.scissorCount = 1;

      VkPipelineRasterizationStateCreateInfo rasterizer{};
      rasterizer.sType = VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO;
      rasterizer.polygonMode = VK_POLYGON_MODE_FILL;
      rasterizer.lineWidth = 1.0f;
      rasterizer.cullMode = VK_CULL_MODE_NONE;
      rasterizer.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;

      VkPipelineMultisampleStateCreateInfo multisampling{};
      multisampling.sType = VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO;
      multisampling.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;

      VkPipelineDepthStencilStateCreateInfo depthStencil{};
      depthStencil.sType = VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO;
      depthStencil.depthTestEnable = VK_TRUE;
      depthStencil.depthWriteEnable = VK_TRUE;
      depthStencil.depthCompareOp = VK_COMPARE_OP_LESS_OR_EQUAL;

      VkPipelineColorBlendAttachmentState colorBlendAttachment{};
      colorBlendAttachment.colorWriteMask = VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT |
                                            VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT;
      colorBlendAttachment.blendEnable = VK_TRUE;
      colorBlendAttachment.srcColorBlendFactor = VK_BLEND_FACTOR_SRC_ALPHA;
      colorBlendAttachment.dstColorBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
      colorBlendAttachment.colorBlendOp = VK_BLEND_OP_ADD;
      colorBlendAttachment.srcAlphaBlendFactor = VK_BLEND_FACTOR_ONE;
      colorBlendAttachment.dstAlphaBlendFactor = VK_BLEND_FACTOR_ZERO;
      colorBlendAttachment.alphaBlendOp = VK_BLEND_OP_ADD;

      VkPipelineColorBlendStateCreateInfo colorBlending{};
      colorBlending.sType = VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO;
      colorBlending.attachmentCount = 1;
      colorBlending.pAttachments = &colorBlendAttachment;

      VkDynamicState dynamicStates[] = {VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR};
      VkPipelineDynamicStateCreateInfo dynamicState{};
      dynamicState.sType = VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO;
      dynamicState.dynamicStateCount = 2;
      dynamicState.pDynamicStates = dynamicStates;

      auto createPipe = [&](VkPrimitiveTopology topo, VkPipeline& outPipe) {
        VkPipelineInputAssemblyStateCreateInfo inputAssembly{};
        inputAssembly.sType = VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO;
        inputAssembly.topology = topo;

        VkGraphicsPipelineCreateInfo pipelineInfo{};
        pipelineInfo.sType = VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO;
        pipelineInfo.stageCount = 2;
        pipelineInfo.pStages = shaderStages;
        pipelineInfo.pVertexInputState = &vertexInput;
        pipelineInfo.pInputAssemblyState = &inputAssembly;
        pipelineInfo.pViewportState = &viewportState;
        pipelineInfo.pRasterizationState = &rasterizer;
        pipelineInfo.pMultisampleState = &multisampling;
        pipelineInfo.pDepthStencilState = &depthStencil;
        pipelineInfo.pColorBlendState = &colorBlending;
        pipelineInfo.pDynamicState = &dynamicState;
        pipelineInfo.layout = fPipelineLayout;
        pipelineInfo.renderPass = fRenderPass;
        pipelineInfo.subpass = 0;
        vkCreateGraphicsPipelines(fDevice, VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, &outPipe);
      };

      createPipe(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST, fTrianglePipeline);
      createPipe(VK_PRIMITIVE_TOPOLOGY_LINE_LIST, fLinePipeline);
      createPipe(VK_PRIMITIVE_TOPOLOGY_POINT_LIST, fPointPipeline);
    }

    bool CreateBufferAndMemory(VkDeviceSize size, VkBufferUsageFlags usage,
                               VkMemoryPropertyFlags props, VkBuffer& buf, VkDeviceMemory& mem)
    {
      if (fDevice == VK_NULL_HANDLE) return false;
      VkBufferCreateInfo bufInfo{};
      bufInfo.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
      bufInfo.size = size;
      bufInfo.usage = usage;
      bufInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
      if (vkCreateBuffer(fDevice, &bufInfo, nullptr, &buf) != VK_SUCCESS) return false;

      VkMemoryRequirements memReq{};
      vkGetBufferMemoryRequirements(fDevice, buf, &memReq);
      uint32_t memType = FindMemoryType(memReq.memoryTypeBits, props);
      if (memType == UINT32_MAX)
      {
        vkDestroyBuffer(fDevice, buf, nullptr);
        buf = VK_NULL_HANDLE;
        return false;
      }

      VkMemoryAllocateInfo allocInfo{};
      allocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
      allocInfo.allocationSize = memReq.size;
      allocInfo.memoryTypeIndex = memType;
      if (vkAllocateMemory(fDevice, &allocInfo, nullptr, &mem) != VK_SUCCESS)
      {
        vkDestroyBuffer(fDevice, buf, nullptr);
        buf = VK_NULL_HANDLE;
        return false;
      }
      vkBindBufferMemory(fDevice, buf, mem, 0);
      return true;
    }

    void CreateVulkanBuffers()
    {
      if (fDevice == VK_NULL_HANDLE) return;
      const VkDeviceSize vbSize =
        static_cast<VkDeviceSize>(kMaxStaticVertices + kMaxDynamicVertices) * sizeof(G4VulkanVertex);
      if (CreateBufferAndMemory(vbSize, VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
                                VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
                                  VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                                fVertexBuffer, fVertexMemory))
      {
        void* ptr = nullptr;
        if (vkMapMemory(fDevice, fVertexMemory, 0, vbSize, 0, &ptr) == VK_SUCCESS)
        {
          fMappedVertices = static_cast<G4VulkanVertex*>(ptr);
        }
      }

      if (CreateBufferAndMemory(kStagingBufferSize, VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                                VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
                                  VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                                fStagingBuffer, fStagingMemory))
      {
        void* ptr = nullptr;
        if (vkMapMemory(fDevice, fStagingMemory, 0, kStagingBufferSize, 0, &ptr) == VK_SUCCESS)
        {
          fMappedStaging = static_cast<uint8_t*>(ptr);
        }
      }
    }

    void DestroyGpuFramebuffer()
    {
      if (fDevice == VK_NULL_HANDLE) return;
      WaitGpuIdle();
      if (fGpuFramebuffer != VK_NULL_HANDLE)
      {
        vkDestroyFramebuffer(fDevice, fGpuFramebuffer, nullptr);
        fGpuFramebuffer = VK_NULL_HANDLE;
      }
      if (fGpuColorView != VK_NULL_HANDLE)
      {
        vkDestroyImageView(fDevice, fGpuColorView, nullptr);
        fGpuColorView = VK_NULL_HANDLE;
      }
      if (fGpuColorImage != VK_NULL_HANDLE)
      {
        vkDestroyImage(fDevice, fGpuColorImage, nullptr);
        fGpuColorImage = VK_NULL_HANDLE;
      }
      if (fGpuColorMemory != VK_NULL_HANDLE)
      {
        vkFreeMemory(fDevice, fGpuColorMemory, nullptr);
        fGpuColorMemory = VK_NULL_HANDLE;
      }
      if (fGpuDepthView != VK_NULL_HANDLE)
      {
        vkDestroyImageView(fDevice, fGpuDepthView, nullptr);
        fGpuDepthView = VK_NULL_HANDLE;
      }
      if (fGpuDepthImage != VK_NULL_HANDLE)
      {
        vkDestroyImage(fDevice, fGpuDepthImage, nullptr);
        fGpuDepthImage = VK_NULL_HANDLE;
      }
      if (fGpuDepthMemory != VK_NULL_HANDLE)
      {
        vkFreeMemory(fDevice, fGpuDepthMemory, nullptr);
        fGpuDepthMemory = VK_NULL_HANDLE;
      }
      fGpuFbWidth = 0;
      fGpuFbHeight = 0;
      fGpuImageInitialised = false;
    }

    void EnsureGpuFramebuffer(VKsizei w, VKsizei h)
    {
      if (fDevice == VK_NULL_HANDLE || fRenderPass == VK_NULL_HANDLE) return;
      if (w <= 0) w = 1;
      if (h <= 0) h = 1;
      if (w == fGpuFbWidth && h == fGpuFbHeight && fGpuFramebuffer != VK_NULL_HANDLE) return;

      DestroyGpuFramebuffer();

      auto createImg = [&](VkFormat format, VkImageUsageFlags usage, VkImageAspectFlags aspect,
                           VkImage& img, VkDeviceMemory& mem, VkImageView& view) -> bool {
        VkImageCreateInfo imgInfo{};
        imgInfo.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
        imgInfo.imageType = VK_IMAGE_TYPE_2D;
        imgInfo.extent.width = static_cast<uint32_t>(w);
        imgInfo.extent.height = static_cast<uint32_t>(h);
        imgInfo.extent.depth = 1;
        imgInfo.mipLevels = 1;
        imgInfo.arrayLayers = 1;
        imgInfo.format = format;
        imgInfo.tiling = VK_IMAGE_TILING_OPTIMAL;
        imgInfo.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
        imgInfo.usage = usage;
        imgInfo.samples = VK_SAMPLE_COUNT_1_BIT;
        imgInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
        if (vkCreateImage(fDevice, &imgInfo, nullptr, &img) != VK_SUCCESS) return false;

        VkMemoryRequirements memReq{};
        vkGetImageMemoryRequirements(fDevice, img, &memReq);
        uint32_t memType = FindMemoryType(memReq.memoryTypeBits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
        if (memType == UINT32_MAX)
        {
          memType = FindMemoryType(memReq.memoryTypeBits, 0);
        }
        if (memType == UINT32_MAX) return false;

        VkMemoryAllocateInfo allocInfo{};
        allocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
        allocInfo.allocationSize = memReq.size;
        allocInfo.memoryTypeIndex = memType;
        if (vkAllocateMemory(fDevice, &allocInfo, nullptr, &mem) != VK_SUCCESS) return false;
        vkBindImageMemory(fDevice, img, mem, 0);

        VkImageViewCreateInfo viewInfo{};
        viewInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
        viewInfo.image = img;
        viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
        viewInfo.format = format;
        viewInfo.subresourceRange.aspectMask = aspect;
        viewInfo.subresourceRange.baseMipLevel = 0;
        viewInfo.subresourceRange.levelCount = 1;
        viewInfo.subresourceRange.baseArrayLayer = 0;
        viewInfo.subresourceRange.layerCount = 1;
        return vkCreateImageView(fDevice, &viewInfo, nullptr, &view) == VK_SUCCESS;
      };

      if (!createImg(VK_FORMAT_R8G8B8A8_UNORM,
                     VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT,
                     VK_IMAGE_ASPECT_COLOR_BIT, fGpuColorImage, fGpuColorMemory, fGpuColorView))
      {
        DestroyGpuFramebuffer();
        return;
      }

      if (!createImg(VK_FORMAT_D32_SFLOAT, VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT,
                     VK_IMAGE_ASPECT_DEPTH_BIT, fGpuDepthImage, fGpuDepthMemory, fGpuDepthView))
      {
        DestroyGpuFramebuffer();
        return;
      }

      VkImageView attachments[2] = {fGpuColorView, fGpuDepthView};
      VkFramebufferCreateInfo fbInfo{};
      fbInfo.sType = VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO;
      fbInfo.renderPass = fRenderPass;
      fbInfo.attachmentCount = 2;
      fbInfo.pAttachments = attachments;
      fbInfo.width = static_cast<uint32_t>(w);
      fbInfo.height = static_cast<uint32_t>(h);
      fbInfo.layers = 1;
      if (vkCreateFramebuffer(fDevice, &fbInfo, nullptr, &fGpuFramebuffer) != VK_SUCCESS)
      {
        DestroyGpuFramebuffer();
        return;
      }

      fGpuFbWidth = w;
      fGpuFbHeight = h;
      fGpuImageInitialised = false;
    }

    void WaitGpuIdle()
    {
      if (fDevice == VK_NULL_HANDLE) return;
      for (uint32_t i = 0; i < kCmdRingSize; ++i)
      {
        if (fFenceSubmitted[i])
        {
          vkWaitForFences(fDevice, 1, &fFenceRing[i], VK_TRUE, UINT64_MAX);
          vkResetFences(fDevice, 1, &fFenceRing[i]);
          fFenceSubmitted[i] = false;
        }
      }
      fDynamicVertexOffset = 0;
    }

    void DestroyVulkan()
    {
      if (fDevice != VK_NULL_HANDLE)
      {
        vkDeviceWaitIdle(fDevice);
        DestroyGpuFramebuffer();
        if (fMappedVertices)
        {
          vkUnmapMemory(fDevice, fVertexMemory);
          fMappedVertices = nullptr;
        }
        if (fVertexBuffer != VK_NULL_HANDLE) vkDestroyBuffer(fDevice, fVertexBuffer, nullptr);
        if (fVertexMemory != VK_NULL_HANDLE) vkFreeMemory(fDevice, fVertexMemory, nullptr);
        if (fMappedStaging)
        {
          vkUnmapMemory(fDevice, fStagingMemory);
          fMappedStaging = nullptr;
        }
        if (fStagingBuffer != VK_NULL_HANDLE) vkDestroyBuffer(fDevice, fStagingBuffer, nullptr);
        if (fStagingMemory != VK_NULL_HANDLE) vkFreeMemory(fDevice, fStagingMemory, nullptr);
        for (uint32_t i = 0; i < kCmdRingSize; ++i)
        {
          if (fFenceRing[i] != VK_NULL_HANDLE) vkDestroyFence(fDevice, fFenceRing[i], nullptr);
        }
        if (fTransferFence != VK_NULL_HANDLE) vkDestroyFence(fDevice, fTransferFence, nullptr);
        if (fTrianglePipeline != VK_NULL_HANDLE) vkDestroyPipeline(fDevice, fTrianglePipeline, nullptr);
        if (fLinePipeline != VK_NULL_HANDLE) vkDestroyPipeline(fDevice, fLinePipeline, nullptr);
        if (fPointPipeline != VK_NULL_HANDLE) vkDestroyPipeline(fDevice, fPointPipeline, nullptr);
        if (fPipelineLayout != VK_NULL_HANDLE) vkDestroyPipelineLayout(fDevice, fPipelineLayout, nullptr);
        if (fVertShader != VK_NULL_HANDLE) vkDestroyShaderModule(fDevice, fVertShader, nullptr);
        if (fFragShader != VK_NULL_HANDLE) vkDestroyShaderModule(fDevice, fFragShader, nullptr);
        if (fRenderPass != VK_NULL_HANDLE) vkDestroyRenderPass(fDevice, fRenderPass, nullptr);
        if (fRenderPassLoad != VK_NULL_HANDLE) vkDestroyRenderPass(fDevice, fRenderPassLoad, nullptr);
        if (fCommandPool != VK_NULL_HANDLE) vkDestroyCommandPool(fDevice, fCommandPool, nullptr);
        vkDestroyDevice(fDevice, nullptr);
        fDevice = VK_NULL_HANDLE;
      }
      if (fInstance != VK_NULL_HANDLE)
      {
        vkDestroyInstance(fInstance, nullptr);
        fInstance = VK_NULL_HANDLE;
      }
    }

    G4VulkanPushConstants BuildPushConstants()
    {
      if (fCachedMvpDirty)
      {
        const Matrix4d m = MultiplyMatrix(fProjectionStack.back(), fModelViewStack.back());
        for (int col = 0; col < 4; ++col)
        {
          fCachedMvp[col * 4 + 0] = static_cast<float>(m[col * 4 + 0]);
          fCachedMvp[col * 4 + 1] = static_cast<float>(m[col * 4 + 1]);
          // Map OpenGL clip-space Z [-w, +w] to Vulkan clip-space Z [0, +w]
          fCachedMvp[col * 4 + 2] = static_cast<float>(0.5 * (m[col * 4 + 2] + m[col * 4 + 3]));
          fCachedMvp[col * 4 + 3] = static_cast<float>(m[col * 4 + 3]);
        }
        const auto& mv = fModelViewStack.back();
        const double lx = fLight0Pos[0], ly = fLight0Pos[1], lz = fLight0Pos[2];
        fCachedLightDir[0] = static_cast<float>(mv[0] * lx + mv[1] * ly + mv[2] * lz);
        fCachedLightDir[1] = static_cast<float>(mv[4] * lx + mv[5] * ly + mv[6] * lz);
        fCachedLightDir[2] = static_cast<float>(mv[8] * lx + mv[9] * ly + mv[10] * lz);
        fCachedLightDir[3] = 0.0f;
        fCachedMvpDirty = false;
      }

      G4VulkanPushConstants pc{};
      std::memcpy(pc.mvp, fCachedMvp, sizeof(fCachedMvp));
      std::memcpy(pc.lightDir, fCachedLightDir, sizeof(fCachedLightDir));
      pc.pointSize = fPointSize;
      pc.lighting = fEnabled.count(VK_LIGHTING) ? 1 : 0;
      return pc;
    }

    void EnqueueDrawCmd(VkPipeline pipe, uint32_t firstVertex, uint32_t vertexCount,
                        const G4VulkanPushConstants& pc)
    {
      if (pipe == VK_NULL_HANDLE || vertexCount == 0) return;
      if (!fPendingDrawCmds.empty())
      {
        auto& back = fPendingDrawCmds.back();
        if (back.pipeline == pipe && back.firstVertex + back.vertexCount == firstVertex &&
            std::memcmp(&back.pc, &pc, sizeof(G4VulkanPushConstants)) == 0)
        {
          back.vertexCount += vertexCount;
          fFrameHasWork = true;
          return;
        }
      }
      fPendingDrawCmds.push_back({pipe, firstVertex, vertexCount, pc});
      fFrameHasWork = true;
    }

    void UploadDisplayListToGpu(G4VulkanDisplayList& dl)
    {
      if (!fMappedVertices) return;
      const uint32_t needed = static_cast<uint32_t>(dl.triangles.size() + dl.lines.size() + dl.points.size());
      if (needed == 0) return;
      if (fStaticVertexCount + needed > kMaxStaticVertices)
      {
        WaitGpuIdle();
        fStaticVertexCount = 0;
      }
      if (!dl.triangles.empty())
      {
        dl.gpuTriFirst = fStaticVertexCount;
        dl.gpuTriCount = static_cast<uint32_t>(dl.triangles.size());
        std::memcpy(fMappedVertices + fStaticVertexCount, dl.triangles.data(),
                    dl.triangles.size() * sizeof(G4VulkanVertex));
        fStaticVertexCount += dl.gpuTriCount;
      }
      if (!dl.lines.empty())
      {
        dl.gpuLineFirst = fStaticVertexCount;
        dl.gpuLineCount = static_cast<uint32_t>(dl.lines.size());
        std::memcpy(fMappedVertices + fStaticVertexCount, dl.lines.data(),
                    dl.lines.size() * sizeof(G4VulkanVertex));
        fStaticVertexCount += dl.gpuLineCount;
      }
      if (!dl.points.empty())
      {
        dl.gpuPointFirst = fStaticVertexCount;
        dl.gpuPointCount = static_cast<uint32_t>(dl.points.size());
        std::memcpy(fMappedVertices + fStaticVertexCount, dl.points.data(),
                    dl.points.size() * sizeof(G4VulkanVertex));
        fStaticVertexCount += dl.gpuPointCount;
      }
    }

    void FlushGpuCommands()
    {
      if (fDevice == VK_NULL_HANDLE || !fFrameHasWork) return;
      EnsureGpuFramebuffer(fFbWidth, fFbHeight);
      if (fGpuFramebuffer == VK_NULL_HANDLE) return;

      const uint32_t slot = fCmdRingIndex;
      fCmdRingIndex = (fCmdRingIndex + 1) % kCmdRingSize;

      if (fFenceSubmitted[slot])
      {
        vkWaitForFences(fDevice, 1, &fFenceRing[slot], VK_TRUE, UINT64_MAX);
        vkResetFences(fDevice, 1, &fFenceRing[slot]);
        fFenceSubmitted[slot] = false;
      }

      VkCommandBuffer cmd = fCmdRing[slot];
      vkResetCommandBuffer(cmd, 0);

      VkCommandBufferBeginInfo beginInfo{};
      beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
      beginInfo.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
      vkBeginCommandBuffer(cmd, &beginInfo);

      const bool doClear = fFrameNeedsClear || !fGpuImageInitialised;
      VkClearValue clearValues[2]{};
      clearValues[0].color = {{fClearColor[0], fClearColor[1], fClearColor[2], fClearColor[3]}};
      clearValues[1].depthStencil = {static_cast<float>(fClearDepth), 0};

      VkRenderPassBeginInfo rpBegin{};
      rpBegin.sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO;
      rpBegin.renderPass = doClear ? fRenderPass : fRenderPassLoad;
      rpBegin.framebuffer = fGpuFramebuffer;
      rpBegin.renderArea.offset = {0, 0};
      rpBegin.renderArea.extent = {static_cast<uint32_t>(fGpuFbWidth),
                                   static_cast<uint32_t>(fGpuFbHeight)};
      rpBegin.clearValueCount = 2;
      rpBegin.pClearValues = clearValues;

      vkCmdBeginRenderPass(cmd, &rpBegin, VK_SUBPASS_CONTENTS_INLINE);

      VkViewport viewport{};
      viewport.x = static_cast<float>(fViewport[0]);
      viewport.y = static_cast<float>(fViewport[1]);
      viewport.width = static_cast<float>(std::max(1, fViewport[2]));
      viewport.height = static_cast<float>(std::max(1, fViewport[3]));
      viewport.minDepth = 0.0f;
      viewport.maxDepth = 1.0f;
      vkCmdSetViewport(cmd, 0, 1, &viewport);

      VkRect2D scissor{};
      scissor.offset = {0, 0};
      scissor.extent = {static_cast<uint32_t>(fGpuFbWidth), static_cast<uint32_t>(fGpuFbHeight)};
      vkCmdSetScissor(cmd, 0, 1, &scissor);

      if (!fPendingDrawCmds.empty() && fVertexBuffer != VK_NULL_HANDLE)
      {
        VkDeviceSize offset = 0;
        vkCmdBindVertexBuffers(cmd, 0, 1, &fVertexBuffer, &offset);
        VkPipeline lastPipe = VK_NULL_HANDLE;
        for (const auto& dc : fPendingDrawCmds)
        {
          if (dc.pipeline == VK_NULL_HANDLE || dc.vertexCount == 0) continue;
          if (dc.pipeline != lastPipe)
          {
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, dc.pipeline);
            lastPipe = dc.pipeline;
          }
          vkCmdPushConstants(cmd, fPipelineLayout, VK_SHADER_STAGE_VERTEX_BIT, 0,
                             sizeof(G4VulkanPushConstants), &dc.pc);
          vkCmdDraw(cmd, dc.vertexCount, 1, dc.firstVertex, 0);
        }
      }

      vkCmdEndRenderPass(cmd);
      vkEndCommandBuffer(cmd);

      VkSubmitInfo submitInfo{};
      submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
      submitInfo.commandBufferCount = 1;
      submitInfo.pCommandBuffers = &cmd;

      if (vkQueueSubmit(fGraphicsQueue, 1, &submitInfo, fFenceRing[slot]) == VK_SUCCESS)
      {
        fFenceSubmitted[slot] = true;
        fGpuImageInitialised = true;
        fGpuColorDirtyForCpu = true;
      }

      fPendingDrawCmds.clear();
      fFrameNeedsClear = false;
      fFrameHasWork = false;
    }

    bool ReadGpuPixels(VKint x, VKint y, VKsizei width, VKsizei height, uint8_t* dstRgba)
    {
      if (fDevice == VK_NULL_HANDLE || !fMappedStaging) return false;
      FlushGpuCommands();
      WaitGpuIdle();
      if (!fGpuImageInitialised || fGpuColorImage == VK_NULL_HANDLE) return false;

      const int cx = std::clamp(x, 0, std::max(0, fGpuFbWidth - 1));
      const int cy = std::clamp(y, 0, std::max(0, fGpuFbHeight - 1));
      const int cw = std::clamp(width, 1, std::max(1, fGpuFbWidth - cx));
      const int ch = std::clamp(height, 1, std::max(1, fGpuFbHeight - cy));
      if (static_cast<VkDeviceSize>(cw) * static_cast<VkDeviceSize>(ch) * 4 > kStagingBufferSize)
      {
        return false;
      }

      vkResetCommandBuffer(fCommandBuffer, 0);
      VkCommandBufferBeginInfo beginInfo{};
      beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
      beginInfo.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
      vkBeginCommandBuffer(fCommandBuffer, &beginInfo);

      VkBufferImageCopy region{};
      region.bufferOffset = 0;
      region.bufferRowLength = 0;
      region.bufferImageHeight = 0;
      region.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
      region.imageSubresource.mipLevel = 0;
      region.imageSubresource.baseArrayLayer = 0;
      region.imageSubresource.layerCount = 1;
      region.imageOffset = {cx, cy, 0};
      region.imageExtent = {static_cast<uint32_t>(cw), static_cast<uint32_t>(ch), 1};

      vkCmdCopyImageToBuffer(fCommandBuffer, fGpuColorImage, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                             fStagingBuffer, 1, &region);
      vkEndCommandBuffer(fCommandBuffer);

      vkResetFences(fDevice, 1, &fTransferFence);
      VkSubmitInfo submitInfo{};
      submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
      submitInfo.commandBufferCount = 1;
      submitInfo.pCommandBuffers = &fCommandBuffer;
      if (vkQueueSubmit(fGraphicsQueue, 1, &submitInfo, fTransferFence) != VK_SUCCESS)
      {
        return false;
      }
      vkWaitForFences(fDevice, 1, &fTransferFence, VK_TRUE, UINT64_MAX);

      std::memcpy(dstRgba, fMappedStaging, static_cast<std::size_t>(cw) * static_cast<std::size_t>(ch) * 4);
      return true;
    }

    void SyncGpuFramebufferToCpu()
    {
      if (fDevice == VK_NULL_HANDLE || !fGpuColorDirtyForCpu) return;
      const std::size_t n = static_cast<std::size_t>(fFbWidth) * static_cast<std::size_t>(fFbHeight);
      if (fColorBuffer.size() < n * 4) fColorBuffer.resize(n * 4, 0);
      if (ReadGpuPixels(0, 0, fFbWidth, fFbHeight, fColorBuffer.data()))
      {
        fGpuColorDirtyForCpu = false;
      }
    }

    void RecordIfCompiling(std::function<void()> op)
    {
      if (fCurrentCompileList != 0)
      {
        fDisplayLists[fCurrentCompileList].stateOps.push_back(op);
      }
      op();
    }

    void ResizeFramebuffer(VKsizei w, VKsizei h)
    {
      if (w <= 0) w = 1;
      if (h <= 0) h = 1;
      if (w == fFbWidth && h == fFbHeight && (fDevice != VK_NULL_HANDLE || !fColorBuffer.empty()))
      {
        return;
      }
      fFbWidth = w;
      fFbHeight = h;
      if (fDevice == VK_NULL_HANDLE)
      {
        const std::size_t n = static_cast<std::size_t>(w) * static_cast<std::size_t>(h);
        fColorBuffer.assign(n * 4, 0);
        fDepthBuffer.assign(n, 1.0f);
      }
    }

    void ClearBuffers(VKbitfield mask)
    {
      EnsureVulkanInitialised();
      ResizeFramebuffer(fViewport[2], fViewport[3]);
      if (fDevice != VK_NULL_HANDLE)
      {
        if (mask & (VK_COLOR_BUFFER_BIT | VK_DEPTH_BUFFER_BIT))
        {
          fFrameNeedsClear = true;
          fFrameHasWork = true;
          fPendingDrawCmds.clear();
        }
        return;
      }
      const std::size_t n = static_cast<std::size_t>(fFbWidth) * static_cast<std::size_t>(fFbHeight);
      if (fColorBuffer.size() < n * 4) fColorBuffer.resize(n * 4, 0);
      if (fDepthBuffer.size() < n) fDepthBuffer.resize(n, 1.0f);
      if (mask & VK_COLOR_BUFFER_BIT)
      {
        const uint8_t r = static_cast<uint8_t>(std::clamp(fClearColor[0], 0.0f, 1.0f) * 255.0f);
        const uint8_t g = static_cast<uint8_t>(std::clamp(fClearColor[1], 0.0f, 1.0f) * 255.0f);
        const uint8_t b = static_cast<uint8_t>(std::clamp(fClearColor[2], 0.0f, 1.0f) * 255.0f);
        const uint8_t a = static_cast<uint8_t>(std::clamp(fClearColor[3], 0.0f, 1.0f) * 255.0f);
        for (std::size_t i = 0; i < n; ++i)
        {
          fColorBuffer[i * 4 + 0] = r;
          fColorBuffer[i * 4 + 1] = g;
          fColorBuffer[i * 4 + 2] = b;
          fColorBuffer[i * 4 + 3] = a;
        }
      }
      if (mask & VK_DEPTH_BUFFER_BIT)
      {
        std::fill(fDepthBuffer.begin(), fDepthBuffer.end(), static_cast<float>(fClearDepth));
      }
    }

    bool IsClippedWorld(const double eyePos[4]) const
    {
      for (int i = 0; i < 6; ++i)
      {
        if (fEnabled.count(VK_CLIP_PLANE0 + i))
        {
          const auto& eq = fClipPlanes[i];
          double d = eq[0] * eyePos[0] + eq[1] * eyePos[1] + eq[2] * eyePos[2] + eq[3] * eyePos[3];
          if (d < 0.0) return true;
        }
      }
      return false;
    }

    struct ProjectedVertex
    {
      double winX{0.0};
      double winY{0.0};
      double winZ{0.0};
      double clipX{0.0};
      double clipY{0.0};
      double clipZ{0.0};
      double clipW{1.0};
      float color[4]{1.0f, 1.0f, 1.0f, 1.0f};
      bool clipped{false};
    };

    ProjectedVertex TransformAndProject(const G4VulkanVertex& v) const
    {
      ProjectedVertex pv{};
      const double in[4] = {v.pos[0], v.pos[1], v.pos[2], 1.0};
      double eye[4];
      TransformVec4(fModelViewStack.back(), in, eye);
      if (IsClippedWorld(eye))
      {
        pv.clipped = true;
        return pv;
      }
      double clip[4];
      TransformVec4(fProjectionStack.back(), eye, clip);
      pv.clipX = clip[0];
      pv.clipY = clip[1];
      pv.clipZ = clip[2];
      pv.clipW = clip[3];

      for (int c = 0; c < 4; ++c) pv.color[c] = v.color[c];
      if (fEnabled.count(VK_LIGHTING))
      {
        const auto& mv = fModelViewStack.back();
        double nx = mv[0] * v.normal[0] + mv[4] * v.normal[1] + mv[8] * v.normal[2];
        double ny = mv[1] * v.normal[0] + mv[5] * v.normal[1] + mv[9] * v.normal[2];
        double nz = mv[2] * v.normal[0] + mv[6] * v.normal[1] + mv[10] * v.normal[2];
        double nlen = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (nlen > 1e-12)
        {
          nx /= nlen;
          ny /= nlen;
          nz /= nlen;
        }
        double lx = fLight0Pos[0], ly = fLight0Pos[1], lz = fLight0Pos[2];
        double llen = std::sqrt(lx * lx + ly * ly + lz * lz);
        if (llen > 1e-12)
        {
          lx /= llen;
          ly /= llen;
          lz /= llen;
        }
        else
        {
          lz = 1.0;
        }
        double ndotl = std::abs(nx * lx + ny * ly + nz * lz);
        float factor = static_cast<float>(std::clamp(0.25 + 0.75 * ndotl, 0.0, 1.0));
        pv.color[0] *= factor;
        pv.color[1] *= factor;
        pv.color[2] *= factor;
      }

      if (std::abs(clip[3]) < 1e-12)
      {
        pv.clipped = true;
        return pv;
      }
      double ndcX = clip[0] / clip[3];
      double ndcY = clip[1] / clip[3];
      double ndcZ = clip[2] / clip[3];
      if (ndcZ < -1.0 || ndcZ > 1.0)
      {
        pv.clipped = true;
      }
      pv.winX = fViewport[0] + (ndcX + 1.0) * 0.5 * fViewport[2];
      pv.winY = fViewport[1] + (ndcY + 1.0) * 0.5 * fViewport[3];
      pv.winZ = (ndcZ + 1.0) * 0.5;
      return pv;
    }

    void WritePixelFast(int x, int y, float z, const float color[4], bool depthTest, bool blend)
    {
      if (x < 0 || x >= fFbWidth || y < 0 || y >= fFbHeight) return;
      const std::size_t idx = static_cast<std::size_t>(y) * static_cast<std::size_t>(fFbWidth) +
                              static_cast<std::size_t>(x);
      if (depthTest)
      {
        const float curZ = fDepthBuffer[idx];
        bool pass = true;
        switch (fDepthFunc)
        {
          case VK_NEVER: pass = false; break;
          case VK_LESS: pass = (z < curZ); break;
          case VK_EQUAL: pass = (std::abs(z - curZ) < 1e-6f); break;
          case VK_LEQUAL: pass = (z <= curZ); break;
          case VK_GREATER: pass = (z > curZ); break;
          case VK_NOTEQUAL: pass = (std::abs(z - curZ) >= 1e-6f); break;
          case VK_GEQUAL: pass = (z >= curZ); break;
          case VK_ALWAYS:
          default: pass = true; break;
        }
        if (!pass) return;
        if (fDepthMask) fDepthBuffer[idx] = z;
      }
      if (!fColorMask[0] && !fColorMask[1] && !fColorMask[2] && !fColorMask[3]) return;

      float r = std::clamp(color[0], 0.0f, 1.0f);
      float g = std::clamp(color[1], 0.0f, 1.0f);
      float b = std::clamp(color[2], 0.0f, 1.0f);
      float a = std::clamp(color[3], 0.0f, 1.0f);

      if (blend)
      {
        const float dstR = fColorBuffer[idx * 4 + 0] / 255.0f;
        const float dstG = fColorBuffer[idx * 4 + 1] / 255.0f;
        const float dstB = fColorBuffer[idx * 4 + 2] / 255.0f;
        const float dstA = fColorBuffer[idx * 4 + 3] / 255.0f;
        r = r * a + dstR * (1.0f - a);
        g = g * a + dstG * (1.0f - a);
        b = b * a + dstB * (1.0f - a);
        a = std::max(a, dstA);
      }

      if (fColorMask[0]) fColorBuffer[idx * 4 + 0] = static_cast<uint8_t>(r * 255.0f);
      if (fColorMask[1]) fColorBuffer[idx * 4 + 1] = static_cast<uint8_t>(g * 255.0f);
      if (fColorMask[2]) fColorBuffer[idx * 4 + 2] = static_cast<uint8_t>(b * 255.0f);
      if (fColorMask[3]) fColorBuffer[idx * 4 + 3] = static_cast<uint8_t>(a * 255.0f);
    }

    void WritePixel(int x, int y, float z, const float color[4])
    {
      WritePixelFast(x, y, z, color, fEnabled.count(VK_DEPTH_TEST) != 0,
                     fEnabled.count(VK_BLEND) != 0);
    }

    void RasterizePoint(const ProjectedVertex& p)
    {
      if (p.clipped) return;
      const bool depthTest = (fEnabled.count(VK_DEPTH_TEST) != 0);
      const bool blend = (fEnabled.count(VK_BLEND) != 0);
      const int rad = std::max(0, static_cast<int>(fPointSize * 0.5f));
      const int cx = static_cast<int>(std::round(p.winX));
      const int cy = static_cast<int>(std::round(p.winY));
      for (int dy = -rad; dy <= rad; ++dy)
      {
        for (int dx = -rad; dx <= rad; ++dx)
        {
          WritePixelFast(cx + dx, cy + dy, static_cast<float>(p.winZ), p.color, depthTest, blend);
        }
      }
    }

    void RasterizeLine(const ProjectedVertex& p0, const ProjectedVertex& p1)
    {
      if (p0.clipped && p1.clipped) return;
      const bool depthTest = (fEnabled.count(VK_DEPTH_TEST) != 0);
      const bool blend = (fEnabled.count(VK_BLEND) != 0);
      double x0 = p0.winX, y0 = p0.winY, z0 = p0.winZ;
      double x1 = p1.winX, y1 = p1.winY, z1 = p1.winZ;
      const double dx = x1 - x0;
      const double dy = y1 - y0;
      const int steps = std::max({1, static_cast<int>(std::ceil(std::max(std::abs(dx), std::abs(dy))))});
      if (steps > 8192) return;
      const int halfW = std::max(0, static_cast<int>((fLineWidth - 1.0f) * 0.5f));
      for (int i = 0; i <= steps; ++i)
      {
        const double t = static_cast<double>(i) / static_cast<double>(steps);
        const int x = static_cast<int>(std::round(x0 + t * dx));
        const int y = static_cast<int>(std::round(y0 + t * dy));
        const float z = static_cast<float>(z0 + t * (z1 - z0));
        float col[4];
        for (int c = 0; c < 4; ++c)
          col[c] = static_cast<float>(p0.color[c] + t * (p1.color[c] - p0.color[c]));
        for (int wy = -halfW; wy <= halfW; ++wy)
        {
          for (int wx = -halfW; wx <= halfW; ++wx)
          {
            WritePixelFast(x + wx, y + wy, z, col, depthTest, blend);
          }
        }
      }
    }

    void RasterizeTriangle(const ProjectedVertex& v0, const ProjectedVertex& v1,
                           const ProjectedVertex& v2)
    {
      if (v0.clipped && v1.clipped && v2.clipped) return;
      if (fPolygonMode == VK_LINE)
      {
        RasterizeLine(v0, v1);
        RasterizeLine(v1, v2);
        RasterizeLine(v2, v0);
        return;
      }
      if (fPolygonMode == VK_POINT)
      {
        RasterizePoint(v0);
        RasterizePoint(v1);
        RasterizePoint(v2);
        return;
      }

      const int minX = std::max(0, static_cast<int>(std::floor(std::min({v0.winX, v1.winX, v2.winX}))));
      const int maxX = std::min(fFbWidth - 1, static_cast<int>(std::ceil(std::max({v0.winX, v1.winX, v2.winX}))));
      const int minY = std::max(0, static_cast<int>(std::floor(std::min({v0.winY, v1.winY, v2.winY}))));
      const int maxY = std::min(fFbHeight - 1, static_cast<int>(std::ceil(std::max({v0.winY, v1.winY, v2.winY}))));
      if (minX > maxX || minY > maxY) return;

      const double area = (v1.winX - v0.winX) * (v2.winY - v0.winY) -
                          (v1.winY - v0.winY) * (v2.winX - v0.winX);
      if (std::abs(area) < 1e-8) return;
      const double invArea = 1.0 / area;

      const bool depthTest = (fEnabled.count(VK_DEPTH_TEST) != 0);
      const bool blend = (fEnabled.count(VK_BLEND) != 0);
      const float polyOffset =
        fEnabled.count(VK_POLYGON_OFFSET_FILL) ? (fPolyOffsetUnits * 1e-5f) : 0.0f;

      for (int y = minY; y <= maxY; ++y)
      {
        const double py = y + 0.5;
        for (int x = minX; x <= maxX; ++x)
        {
          const double px = x + 0.5;
          const double w0 = ((v1.winX - px) * (v2.winY - py) - (v1.winY - py) * (v2.winX - px)) * invArea;
          const double w1 = ((v2.winX - px) * (v0.winY - py) - (v2.winY - py) * (v0.winX - px)) * invArea;
          const double w2 = 1.0 - w0 - w1;
          if (w0 >= -1e-5 && w1 >= -1e-5 && w2 >= -1e-5)
          {
            const float z = static_cast<float>(w0 * v0.winZ + w1 * v1.winZ + w2 * v2.winZ) + polyOffset;
            float col[4] = {
              static_cast<float>(w0 * v0.color[0] + w1 * v1.color[0] + w2 * v2.color[0]),
              static_cast<float>(w0 * v0.color[1] + w1 * v1.color[1] + w2 * v2.color[1]),
              static_cast<float>(w0 * v0.color[2] + w1 * v1.color[2] + w2 * v2.color[2]),
              static_cast<float>(w0 * v0.color[3] + w1 * v1.color[3] + w2 * v2.color[3])};
            WritePixelFast(x, y, z, col, depthTest, blend);
          }
        }
      }
    }

    void AppendFeedbackVertex(const ProjectedVertex& pv)
    {
      if (!fFeedbackBuffer) return;
      auto push = [&](VKfloat val) {
        if (fFeedbackOffset < fFeedbackSize) fFeedbackBuffer[fFeedbackOffset++] = val;
      };
      push(static_cast<VKfloat>(pv.winX));
      push(static_cast<VKfloat>(pv.winY));
      if (fFeedbackType != VK_2D) push(static_cast<VKfloat>(pv.winZ));
      if (fFeedbackType == VK_3D_COLOR || fFeedbackType == VK_3D_COLOR_TEXTURE ||
          fFeedbackType == VK_4D_COLOR_TEXTURE)
      {
        for (int c = 0; c < 4; ++c) push(pv.color[c]);
      }
    }

    void RecordSelectHit(double zMin, double zMax)
    {
      if (!fSelectBuffer || fSelectBufferSize <= 0) return;
      const uint32_t z1 = static_cast<uint32_t>(std::clamp(zMin, 0.0, 1.0) * 4294967295.0);
      const uint32_t z2 = static_cast<uint32_t>(std::clamp(zMax, 0.0, 1.0) * 4294967295.0);

      if (fLastHitRecordOffset >= 0 &&
          fSelectOffset == fLastHitRecordOffset + 3 + static_cast<VKint>(fNameStack.size()))
      {
        bool sameStack = true;
        for (std::size_t i = 0; i < fNameStack.size(); ++i)
        {
          if (fSelectBuffer[fLastHitRecordOffset + 3 + i] != fNameStack[i])
          {
            sameStack = false;
            break;
          }
        }
        if (sameStack)
        {
          fSelectBuffer[fLastHitRecordOffset + 1] = std::min(fSelectBuffer[fLastHitRecordOffset + 1], z1);
          fSelectBuffer[fLastHitRecordOffset + 2] = std::max(fSelectBuffer[fLastHitRecordOffset + 2], z2);
          return;
        }
      }

      if (fSelectOffset + 3 + static_cast<VKint>(fNameStack.size()) <= fSelectBufferSize)
      {
        fLastHitRecordOffset = fSelectOffset;
        fSelectBuffer[fSelectOffset++] = static_cast<VKuint>(fNameStack.size());
        fSelectBuffer[fSelectOffset++] = z1;
        fSelectBuffer[fSelectOffset++] = z2;
        for (VKuint n : fNameStack)
        {
          fSelectBuffer[fSelectOffset++] = n;
        }
        ++fSelectHits;
      }
    }

    static bool InPickVolume(const ProjectedVertex& pv)
    {
      if (pv.clipped || std::abs(pv.clipW) < 1e-12) return false;
      const double nx = pv.clipX / pv.clipW;
      const double ny = pv.clipY / pv.clipW;
      const double nz = pv.clipZ / pv.clipW;
      return (nx >= -1.0 && nx <= 1.0 && ny >= -1.0 && ny <= 1.0 && nz >= -1.0 && nz <= 1.0);
    }

    static bool SegmentIntersectsPickBox(const ProjectedVertex& a, const ProjectedVertex& b)
    {
      if (a.clipped && b.clipped) return false;
      if (InPickVolume(a) || InPickVolume(b)) return true;
      if (std::abs(a.clipW) < 1e-12 || std::abs(b.clipW) < 1e-12) return false;
      double x0 = a.clipX / a.clipW, y0 = a.clipY / a.clipW;
      double x1 = b.clipX / b.clipW, y1 = b.clipY / b.clipW;
      double t0 = 0.0, t1 = 1.0;
      double dx = x1 - x0, dy = y1 - y0;
      auto clipTest = [&](double p, double q) {
        if (std::abs(p) < 1e-12) return q >= 0.0;
        double r = q / p;
        if (p < 0.0)
        {
          if (r > t1) return false;
          if (r > t0) t0 = r;
        }
        else
        {
          if (r < t0) return false;
          if (r < t1) t1 = r;
        }
        return true;
      };
      return clipTest(-dx, x0 + 1.0) && clipTest(dx, 1.0 - x0) && clipTest(-dy, y0 + 1.0) &&
             clipTest(dy, 1.0 - y0);
    }

    static bool TriangleIntersectsPickBox(const ProjectedVertex& a, const ProjectedVertex& b,
                                          const ProjectedVertex& c)
    {
      if (SegmentIntersectsPickBox(a, b) || SegmentIntersectsPickBox(b, c) ||
          SegmentIntersectsPickBox(c, a))
        return true;
      if (std::abs(a.clipW) < 1e-12 || std::abs(b.clipW) < 1e-12 || std::abs(c.clipW) < 1e-12)
        return false;
      double x0 = a.clipX / a.clipW, y0 = a.clipY / a.clipW;
      double x1 = b.clipX / b.clipW, y1 = b.clipY / b.clipW;
      double x2 = c.clipX / c.clipW, y2 = c.clipY / c.clipW;
      double area = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
      if (std::abs(area) < 1e-12) return false;
      double w0 = (x1 * y2 - y1 * x2) / area;
      double w1 = (x2 * y0 - y2 * x0) / area;
      double w2 = 1.0 - w0 - w1;
      return (w0 >= 0.0 && w1 >= 0.0 && w2 >= 0.0);
    }

    static bool SamePosition(const G4VulkanVertex& a, const G4VulkanVertex& b)
    {
      return a.pos[0] == b.pos[0] && a.pos[1] == b.pos[1] && a.pos[2] == b.pos[2];
    }

    template <typename PointFn, typename LineFn, typename TriFn>
    void DecomposeCurrentPrimitive(PointFn&& onPoint, LineFn&& onLine, TriFn&& onTri) const
    {
      const std::size_t n = fCurrentVertices.size();
      auto emitTri = [&](const G4VulkanVertex& a, const G4VulkanVertex& b, const G4VulkanVertex& c) {
        if (SamePosition(a, c) || SamePosition(a, b) || SamePosition(b, c)) return;
        if (fPolygonMode == VK_LINE)
        {
          onLine(a, b);
          onLine(b, c);
          onLine(c, a);
        }
        else if (fPolygonMode == VK_POINT)
        {
          onPoint(a);
          onPoint(b);
          onPoint(c);
        }
        else
        {
          onTri(a, b, c);
        }
      };

      switch (fCurrentTopology)
      {
        case VK_POINTS:
          for (std::size_t i = 0; i < n; ++i) onPoint(fCurrentVertices[i]);
          break;
        case VK_LINES:
          for (std::size_t i = 0; i + 1 < n; i += 2) onLine(fCurrentVertices[i], fCurrentVertices[i + 1]);
          break;
        case VK_LINE_STRIP:
          for (std::size_t i = 0; i + 1 < n; ++i) onLine(fCurrentVertices[i], fCurrentVertices[i + 1]);
          break;
        case VK_LINE_LOOP:
          for (std::size_t i = 0; i + 1 < n; ++i) onLine(fCurrentVertices[i], fCurrentVertices[i + 1]);
          if (n >= 2) onLine(fCurrentVertices[n - 1], fCurrentVertices[0]);
          break;
        case VK_TRIANGLES:
          for (std::size_t i = 0; i + 2 < n; i += 3)
            emitTri(fCurrentVertices[i], fCurrentVertices[i + 1], fCurrentVertices[i + 2]);
          break;
        case VK_TRIANGLE_STRIP:
          for (std::size_t i = 0; i + 2 < n; ++i)
          {
            if (i % 2 == 0)
              emitTri(fCurrentVertices[i], fCurrentVertices[i + 1], fCurrentVertices[i + 2]);
            else
              emitTri(fCurrentVertices[i + 1], fCurrentVertices[i], fCurrentVertices[i + 2]);
          }
          break;
        case VK_TRIANGLE_FAN:
        case VK_POLYGON:
          for (std::size_t i = 1; i + 1 < n; ++i)
            emitTri(fCurrentVertices[0], fCurrentVertices[i], fCurrentVertices[i + 1]);
          break;
        case VK_QUADS:
          for (std::size_t i = 0; i + 3 < n; i += 4)
          {
            emitTri(fCurrentVertices[i], fCurrentVertices[i + 1], fCurrentVertices[i + 2]);
            if (!SamePosition(fCurrentVertices[i], fCurrentVertices[i + 3]))
            {
              emitTri(fCurrentVertices[i], fCurrentVertices[i + 2], fCurrentVertices[i + 3]);
            }
          }
          break;
        case VK_QUAD_STRIP:
          for (std::size_t i = 0; i + 3 < n; i += 2)
          {
            emitTri(fCurrentVertices[i], fCurrentVertices[i + 1], fCurrentVertices[i + 2]);
            emitTri(fCurrentVertices[i + 1], fCurrentVertices[i + 3], fCurrentVertices[i + 2]);
          }
          break;
        default:
          break;
      }
    }

    void AppendPrimitiveToDisplayList(G4VulkanDisplayList& dl)
    {
      for (const auto& v : fCurrentVertices) dl.AddBounds(v);
      DecomposeCurrentPrimitive(
        [&](const G4VulkanVertex& p) { dl.points.push_back(p); },
        [&](const G4VulkanVertex& a, const G4VulkanVertex& b) {
          dl.lines.push_back(a);
          dl.lines.push_back(b);
        },
        [&](const G4VulkanVertex& a, const G4VulkanVertex& b, const G4VulkanVertex& c) {
          dl.triangles.push_back(a);
          dl.triangles.push_back(b);
          dl.triangles.push_back(c);
        });
    }

    void ExecuteDisplayList(const G4VulkanDisplayList& dl)
    {
      for (const auto& op : dl.stateOps) op();
      if (dl.hasLightingState)
      {
        if (dl.lighting)
          fEnabled.insert(VK_LIGHTING);
        else
          fEnabled.erase(VK_LIGHTING);
      }

      if (fRenderMode == VK_RENDER)
      {
        if (fDevice != VK_NULL_HANDLE)
        {
          const G4VulkanPushConstants pc = BuildPushConstants();
          if (dl.gpuTriCount > 0) EnqueueDrawCmd(fTrianglePipeline, dl.gpuTriFirst, dl.gpuTriCount, pc);
          if (dl.gpuLineCount > 0) EnqueueDrawCmd(fLinePipeline, dl.gpuLineFirst, dl.gpuLineCount, pc);
          if (dl.gpuPointCount > 0) EnqueueDrawCmd(fPointPipeline, dl.gpuPointFirst, dl.gpuPointCount, pc);
          return;
        }

        ResizeFramebuffer(fViewport[2], fViewport[3]);
        for (std::size_t i = 0; i + 2 < dl.triangles.size(); i += 3)
        {
          RasterizeTriangle(TransformAndProject(dl.triangles[i]),
                            TransformAndProject(dl.triangles[i + 1]),
                            TransformAndProject(dl.triangles[i + 2]));
        }
        for (std::size_t i = 0; i + 1 < dl.lines.size(); i += 2)
        {
          RasterizeLine(TransformAndProject(dl.lines[i]), TransformAndProject(dl.lines[i + 1]));
        }
        for (const auto& pt : dl.points)
        {
          RasterizePoint(TransformAndProject(pt));
        }
      }
      else if (fRenderMode == VK_SELECT)
      {
        if (!dl.hasBounds) return;
        double minNx = 1e30, maxNx = -1e30;
        double minNy = 1e30, maxNy = -1e30;
        double minNz = 1e30, maxNz = -1e30;
        double minWinZ = 1.0, maxWinZ = 0.0;
        bool anyCornerIn = false;
        int validCorners = 0;
        for (int c = 0; c < 8; ++c)
        {
          G4VulkanVertex cv{};
          cv.pos[0] = (c & 1) ? dl.maxPos[0] : dl.minPos[0];
          cv.pos[1] = (c & 2) ? dl.maxPos[1] : dl.minPos[1];
          cv.pos[2] = (c & 4) ? dl.maxPos[2] : dl.minPos[2];
          ProjectedVertex pv = TransformAndProject(cv);
          if (std::abs(pv.clipW) > 1e-12)
          {
            ++validCorners;
            const double nx = pv.clipX / pv.clipW;
            const double ny = pv.clipY / pv.clipW;
            const double nz = pv.clipZ / pv.clipW;
            minNx = std::min(minNx, nx);
            maxNx = std::max(maxNx, nx);
            minNy = std::min(minNy, ny);
            maxNy = std::max(maxNy, ny);
            minNz = std::min(minNz, nz);
            maxNz = std::max(maxNz, nz);
            minWinZ = std::min(minWinZ, pv.winZ);
            maxWinZ = std::max(maxWinZ, pv.winZ);
            if (InPickVolume(pv)) anyCornerIn = true;
          }
        }
        if (validCorners == 0) return;
        if (maxNx < -1.0 || minNx > 1.0 || maxNy < -1.0 || minNy > 1.0 || maxNz < -1.0 ||
            minNz > 1.0)
        {
          return;
        }
        if (anyCornerIn || (minNx <= 0.0 && maxNx >= 0.0 && minNy <= 0.0 && maxNy >= 0.0))
        {
          RecordSelectHit(std::clamp(minWinZ, 0.0, 1.0), std::clamp(maxWinZ, 0.0, 1.0));
          return;
        }
        for (std::size_t i = 0; i + 2 < dl.triangles.size(); i += 3)
        {
          auto p0 = TransformAndProject(dl.triangles[i]);
          auto p1 = TransformAndProject(dl.triangles[i + 1]);
          auto p2 = TransformAndProject(dl.triangles[i + 2]);
          if (TriangleIntersectsPickBox(p0, p1, p2))
          {
            RecordSelectHit(std::min({p0.winZ, p1.winZ, p2.winZ}),
                            std::max({p0.winZ, p1.winZ, p2.winZ}));
            return;
          }
        }
        for (std::size_t i = 0; i + 1 < dl.lines.size(); i += 2)
        {
          auto p0 = TransformAndProject(dl.lines[i]);
          auto p1 = TransformAndProject(dl.lines[i + 1]);
          if (SegmentIntersectsPickBox(p0, p1))
          {
            RecordSelectHit(std::min(p0.winZ, p1.winZ), std::max(p0.winZ, p1.winZ));
            return;
          }
        }
        for (const auto& pt : dl.points)
        {
          auto p0 = TransformAndProject(pt);
          if (InPickVolume(p0))
          {
            RecordSelectHit(p0.winZ, p0.winZ);
            return;
          }
        }
      }
      else if (fRenderMode == VK_FEEDBACK && fFeedbackBuffer)
      {
        for (std::size_t i = 0; i + 2 < dl.triangles.size(); i += 3)
        {
          auto p0 = TransformAndProject(dl.triangles[i]);
          auto p1 = TransformAndProject(dl.triangles[i + 1]);
          auto p2 = TransformAndProject(dl.triangles[i + 2]);
          if (p0.clipped && p1.clipped && p2.clipped) continue;
          if (fFeedbackOffset < fFeedbackSize) fFeedbackBuffer[fFeedbackOffset++] = VK_POLYGON_TOKEN;
          if (fFeedbackOffset < fFeedbackSize) fFeedbackBuffer[fFeedbackOffset++] = 3.0f;
          AppendFeedbackVertex(p0);
          AppendFeedbackVertex(p1);
          AppendFeedbackVertex(p2);
        }
        for (std::size_t i = 0; i + 1 < dl.lines.size(); i += 2)
        {
          auto p0 = TransformAndProject(dl.lines[i]);
          auto p1 = TransformAndProject(dl.lines[i + 1]);
          if (p0.clipped && p1.clipped) continue;
          if (fFeedbackOffset < fFeedbackSize) fFeedbackBuffer[fFeedbackOffset++] = VK_LINE_TOKEN;
          AppendFeedbackVertex(p0);
          AppendFeedbackVertex(p1);
        }
        for (const auto& pt : dl.points)
        {
          auto p0 = TransformAndProject(pt);
          if (p0.clipped) continue;
          if (fFeedbackOffset < fFeedbackSize) fFeedbackBuffer[fFeedbackOffset++] = VK_POINT_TOKEN;
          AppendFeedbackVertex(p0);
        }
      }
    }

    void EndPrimitive()
    {
      if (!fInBegin) return;
      fInBegin = false;
      if (fCurrentVertices.empty()) return;

      if (fCurrentCompileList != 0)
      {
        AppendPrimitiveToDisplayList(fDisplayLists[fCurrentCompileList]);
        if (fCompileMode != VK_COMPILE_AND_EXECUTE)
        {
          fCurrentVertices.clear();
          return;
        }
      }

      if (fRenderMode == VK_RENDER && fDevice != VK_NULL_HANDLE && fMappedVertices)
      {
        const std::size_t maxOut = fCurrentVertices.size() * 3;
        if (fDynamicVertexOffset + maxOut > kMaxDynamicVertices)
        {
          FlushGpuCommands();
          WaitGpuIdle();
        }
        const uint32_t baseVert = kMaxStaticVertices + fDynamicVertexOffset;
        uint32_t triCount = 0;
        uint32_t lineCount = 0;
        uint32_t pointCount = 0;
        G4VulkanVertex* outPtr = fMappedVertices + baseVert;

        DecomposeCurrentPrimitive(
          [&](const G4VulkanVertex& p) { outPtr[pointCount++] = p; },
          [&](const G4VulkanVertex& a, const G4VulkanVertex& b) {
            outPtr[lineCount++] = a;
            outPtr[lineCount++] = b;
          },
          [&](const G4VulkanVertex& a, const G4VulkanVertex& b, const G4VulkanVertex& c) {
            outPtr[triCount++] = a;
            outPtr[triCount++] = b;
            outPtr[triCount++] = c;
          });

        const G4VulkanPushConstants pc = BuildPushConstants();
        if (triCount > 0)
        {
          fDynamicVertexOffset += triCount;
          EnqueueDrawCmd(fTrianglePipeline, baseVert, triCount, pc);
        }
        else if (lineCount > 0)
        {
          fDynamicVertexOffset += lineCount;
          EnqueueDrawCmd(fLinePipeline, baseVert, lineCount, pc);
        }
        else if (pointCount > 0)
        {
          fDynamicVertexOffset += pointCount;
          EnqueueDrawCmd(fPointPipeline, baseVert, pointCount, pc);
        }
        fCurrentVertices.clear();
        return;
      }

      ResizeFramebuffer(fViewport[2], fViewport[3]);

      std::vector<ProjectedVertex> proj;
      proj.reserve(fCurrentVertices.size());
      for (const auto& v : fCurrentVertices)
      {
        proj.push_back(TransformAndProject(v));
      }

      auto emitPoint = [&](const ProjectedVertex& p) {
        if (fRenderMode == VK_RENDER)
        {
          RasterizePoint(p);
        }
        else if (fRenderMode == VK_SELECT)
        {
          if (InPickVolume(p)) RecordSelectHit(p.winZ, p.winZ);
        }
        else if (fRenderMode == VK_FEEDBACK && fFeedbackBuffer && !p.clipped)
        {
          if (fFeedbackOffset < fFeedbackSize) fFeedbackBuffer[fFeedbackOffset++] = VK_POINT_TOKEN;
          AppendFeedbackVertex(p);
        }
      };

      auto emitLine = [&](const ProjectedVertex& p0, const ProjectedVertex& p1) {
        if (fRenderMode == VK_RENDER)
        {
          RasterizeLine(p0, p1);
        }
        else if (fRenderMode == VK_SELECT)
        {
          if (SegmentIntersectsPickBox(p0, p1))
          {
            RecordSelectHit(std::min(p0.winZ, p1.winZ), std::max(p0.winZ, p1.winZ));
          }
        }
        else if (fRenderMode == VK_FEEDBACK && fFeedbackBuffer && !(p0.clipped && p1.clipped))
        {
          if (fFeedbackOffset < fFeedbackSize) fFeedbackBuffer[fFeedbackOffset++] = VK_LINE_TOKEN;
          AppendFeedbackVertex(p0);
          AppendFeedbackVertex(p1);
        }
      };

      auto emitTriangle = [&](const ProjectedVertex& p0, const ProjectedVertex& p1,
                              const ProjectedVertex& p2) {
        if (fRenderMode == VK_RENDER)
        {
          RasterizeTriangle(p0, p1, p2);
        }
        else if (fRenderMode == VK_SELECT)
        {
          if (TriangleIntersectsPickBox(p0, p1, p2))
          {
            RecordSelectHit(std::min({p0.winZ, p1.winZ, p2.winZ}),
                            std::max({p0.winZ, p1.winZ, p2.winZ}));
          }
        }
        else if (fRenderMode == VK_FEEDBACK && fFeedbackBuffer &&
                 !(p0.clipped && p1.clipped && p2.clipped))
        {
          if (fPolygonMode == VK_LINE)
          {
            emitLine(p0, p1);
            emitLine(p1, p2);
            emitLine(p2, p0);
          }
          else
          {
            if (fFeedbackOffset < fFeedbackSize) fFeedbackBuffer[fFeedbackOffset++] = VK_POLYGON_TOKEN;
            if (fFeedbackOffset < fFeedbackSize) fFeedbackBuffer[fFeedbackOffset++] = 3.0f;
            AppendFeedbackVertex(p0);
            AppendFeedbackVertex(p1);
            AppendFeedbackVertex(p2);
          }
        }
      };

      const std::size_t n = proj.size();
      switch (fCurrentTopology)
      {
        case VK_POINTS:
          for (std::size_t i = 0; i < n; ++i) emitPoint(proj[i]);
          break;
        case VK_LINES:
          for (std::size_t i = 0; i + 1 < n; i += 2) emitLine(proj[i], proj[i + 1]);
          break;
        case VK_LINE_STRIP:
          for (std::size_t i = 0; i + 1 < n; ++i) emitLine(proj[i], proj[i + 1]);
          break;
        case VK_LINE_LOOP:
          for (std::size_t i = 0; i + 1 < n; ++i) emitLine(proj[i], proj[i + 1]);
          if (n >= 2) emitLine(proj[n - 1], proj[0]);
          break;
        case VK_TRIANGLES:
          for (std::size_t i = 0; i + 2 < n; i += 3) emitTriangle(proj[i], proj[i + 1], proj[i + 2]);
          break;
        case VK_TRIANGLE_STRIP:
          for (std::size_t i = 0; i + 2 < n; ++i)
          {
            if (i % 2 == 0)
              emitTriangle(proj[i], proj[i + 1], proj[i + 2]);
            else
              emitTriangle(proj[i + 1], proj[i], proj[i + 2]);
          }
          break;
        case VK_TRIANGLE_FAN:
        case VK_POLYGON:
          for (std::size_t i = 1; i + 1 < n; ++i) emitTriangle(proj[0], proj[i], proj[i + 1]);
          break;
        case VK_QUADS:
          for (std::size_t i = 0; i + 3 < n; i += 4)
          {
            emitTriangle(proj[i], proj[i + 1], proj[i + 2]);
            if (!SamePosition(fCurrentVertices[i], fCurrentVertices[i + 3]))
            {
              emitTriangle(proj[i], proj[i + 2], proj[i + 3]);
            }
          }
          break;
        case VK_QUAD_STRIP:
          for (std::size_t i = 0; i + 3 < n; i += 2)
          {
            emitTriangle(proj[i], proj[i + 1], proj[i + 2]);
            emitTriangle(proj[i + 1], proj[i + 3], proj[i + 2]);
          }
          break;
        default:
          break;
      }

      fCurrentVertices.clear();
    }

    // Vulkan hardware handles
    bool fVulkanInitAttempted{false};
    VkInstance fInstance{VK_NULL_HANDLE};
    VkPhysicalDevice fPhysicalDevice{VK_NULL_HANDLE};
    VkPhysicalDeviceMemoryProperties fMemProperties{};
    uint32_t fGraphicsQueueFamily{0};
    VkDevice fDevice{VK_NULL_HANDLE};
    VkQueue fGraphicsQueue{VK_NULL_HANDLE};
    VkCommandPool fCommandPool{VK_NULL_HANDLE};
    VkCommandBuffer fCommandBuffer{VK_NULL_HANDLE};
    std::array<VkCommandBuffer, kCmdRingSize> fCmdRing{};
    std::array<VkFence, kCmdRingSize> fFenceRing{};
    std::array<bool, kCmdRingSize> fFenceSubmitted{};
    uint32_t fCmdRingIndex{0};
    VkFence fTransferFence{VK_NULL_HANDLE};
    VkRenderPass fRenderPass{VK_NULL_HANDLE};
    VkRenderPass fRenderPassLoad{VK_NULL_HANDLE};
    VkShaderModule fVertShader{VK_NULL_HANDLE};
    VkShaderModule fFragShader{VK_NULL_HANDLE};
    VkPipelineLayout fPipelineLayout{VK_NULL_HANDLE};
    VkPipeline fTrianglePipeline{VK_NULL_HANDLE};
    VkPipeline fLinePipeline{VK_NULL_HANDLE};
    VkPipeline fPointPipeline{VK_NULL_HANDLE};

    // GPU offscreen framebuffer & persistent mapped buffers
    VkImage fGpuColorImage{VK_NULL_HANDLE};
    VkDeviceMemory fGpuColorMemory{VK_NULL_HANDLE};
    VkImageView fGpuColorView{VK_NULL_HANDLE};
    VkImage fGpuDepthImage{VK_NULL_HANDLE};
    VkDeviceMemory fGpuDepthMemory{VK_NULL_HANDLE};
    VkImageView fGpuDepthView{VK_NULL_HANDLE};
    VkFramebuffer fGpuFramebuffer{VK_NULL_HANDLE};
    VKsizei fGpuFbWidth{0};
    VKsizei fGpuFbHeight{0};
    bool fGpuImageInitialised{false};
    bool fGpuColorDirtyForCpu{false};
    bool fFrameNeedsClear{true};
    bool fFrameHasWork{false};

    VkBuffer fVertexBuffer{VK_NULL_HANDLE};
    VkDeviceMemory fVertexMemory{VK_NULL_HANDLE};
    G4VulkanVertex* fMappedVertices{nullptr};
    uint32_t fStaticVertexCount{0};
    uint32_t fDynamicVertexOffset{0};

    VkBuffer fStagingBuffer{VK_NULL_HANDLE};
    VkDeviceMemory fStagingMemory{VK_NULL_HANDLE};
    uint8_t* fMappedStaging{nullptr};

    std::vector<G4VulkanDrawCmd> fPendingDrawCmds;
    bool fCachedMvpDirty{true};
    float fCachedMvp[16]{};
    float fCachedLightDir[4]{0.0f, 0.0f, 1.0f, 0.0f};

    // State machine members
    VKenum fMatrixMode{VK_MODELVIEW};
    std::vector<Matrix4d> fModelViewStack;
    std::vector<Matrix4d> fProjectionStack;
    std::vector<Matrix4d> fTextureStack;
    VKint fViewport[4]{0, 0, 600, 600};

    float fClearColor[4]{0.0f, 0.0f, 0.0f, 0.0f};
    double fClearDepth{1.0};
    float fCurrentColor[4]{1.0f, 1.0f, 1.0f, 1.0f};
    float fCurrentNormal[3]{0.0f, 0.0f, 1.0f};
    double fRasterPos[3]{0.0, 0.0, 0.0};

    float fPointSize{1.0f};
    float fLineWidth{1.0f};
    VKushort fLineStipplePattern{0xFFFF};
    VKint fLineStippleRepeat{1};
    VKenum fPolygonMode{VK_FILL};
    VKenum fShadeModel{VK_SMOOTH};
    VKenum fFrontFace{VK_CCW};
    VKenum fCullFaceMode{VK_BACK};
    VKenum fDepthFunc{VK_LEQUAL};
    VKboolean fDepthMask{VK_TRUE};
    VKboolean fColorMask[4]{VK_TRUE, VK_TRUE, VK_TRUE, VK_TRUE};
    VKenum fBlendSrc{VK_SRC_ALPHA};
    VKenum fBlendDst{VK_ONE_MINUS_SRC_ALPHA};
    float fPolyOffsetFactor{0.0f};
    float fPolyOffsetUnits{0.0f};
    float fLight0Pos[4]{0.0f, 0.0f, 1.0f, 0.0f};
    std::array<std::array<double, 4>, 6> fClipPlanes{};
    std::unordered_set<VKenum> fEnabled;

    // Primitive assembly
    bool fInBegin{false};
    VKenum fCurrentTopology{VK_TRIANGLES};
    std::vector<G4VulkanVertex> fCurrentVertices;

    // Display list recording store
    VKuint fNextListId{1};
    VKuint fCurrentCompileList{0};
    VKenum fCompileMode{VK_COMPILE};
    VKuint fListBase{0};
    std::unordered_map<VKuint, G4VulkanDisplayList> fDisplayLists;

    // Selection / Picking
    VKenum fRenderMode{VK_RENDER};
    VKuint* fSelectBuffer{nullptr};
    VKsizei fSelectBufferSize{0};
    VKint fSelectOffset{0};
    VKint fSelectHits{0};
    VKint fLastHitRecordOffset{-1};
    std::vector<VKuint> fNameStack;

    // Feedback (gl2ps vector export)
    VKfloat* fFeedbackBuffer{nullptr};
    VKsizei fFeedbackSize{0};
    VKint fFeedbackOffset{0};
    VKenum fFeedbackType{VK_3D_COLOR};

    // Framebuffer
    VKsizei fFbWidth{600};
    VKsizei fFbHeight{600};
    std::vector<uint8_t> fColorBuffer;
    std::vector<float> fDepthBuffer;
};

// State & rendering entrypoints used by G4Vulkan* classes and G4vk2ps
inline void vkEnable(VKenum cap)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fEnabled.insert(cap);
  if (ctx.fCurrentCompileList != 0 && cap == VK_LIGHTING)
  {
    auto& dl = ctx.fDisplayLists[ctx.fCurrentCompileList];
    dl.hasLightingState = true;
    dl.lighting = true;
  }
}

inline void vkDisable(VKenum cap)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fEnabled.erase(cap);
  if (ctx.fCurrentCompileList != 0 && cap == VK_LIGHTING)
  {
    auto& dl = ctx.fDisplayLists[ctx.fCurrentCompileList];
    dl.hasLightingState = true;
    dl.lighting = false;
  }
}

inline VKboolean vkIsEnabled(VKenum cap)
{
  return G4VulkanContext::Instance().fEnabled.count(cap) ? VK_TRUE : VK_FALSE;
}

inline void vkClearColor(VKclampf r, VKclampf g, VKclampf b, VKclampf a)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fClearColor[0] = r;
  ctx.fClearColor[1] = g;
  ctx.fClearColor[2] = b;
  ctx.fClearColor[3] = a;
}

inline void vkClearDepth(VKclampd d) { G4VulkanContext::Instance().fClearDepth = d; }

inline void vkClear(VKbitfield mask) { G4VulkanContext::Instance().ClearBuffers(mask); }

inline void vkViewport(VKint x, VKint y, VKsizei w, VKsizei h)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fViewport[0] = x;
  ctx.fViewport[1] = y;
  ctx.fViewport[2] = std::max(1, w);
  ctx.fViewport[3] = std::max(1, h);
  ctx.ResizeFramebuffer(ctx.fViewport[2], ctx.fViewport[3]);
}

inline void vkScissor(VKint, VKint, VKsizei, VKsizei) {}
inline void vkDrawBuffer(VKenum) {}

inline void vkMatrixMode(VKenum mode)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fMatrixMode = mode;
}

inline void vkLoadIdentity()
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.CurrentMatrix() = G4VulkanContext::IdentityMatrix();
}

inline void vkPushMatrix()
{
  auto& stack = G4VulkanContext::Instance().CurrentStack();
  stack.push_back(stack.back());
}

inline void vkPopMatrix()
{
  auto& stack = G4VulkanContext::Instance().CurrentStack();
  if (stack.size() > 1) stack.pop_back();
}

inline void vkLoadMatrixd(const VKdouble* m)
{
  auto& cur = G4VulkanContext::Instance().CurrentMatrix();
  std::copy(m, m + 16, cur.begin());
}

inline void vkMultMatrixd(const VKdouble* m)
{
  G4VulkanContext::Matrix4d mat{};
  std::copy(m, m + 16, mat.begin());
  auto& cur = G4VulkanContext::Instance().CurrentMatrix();
  cur = G4VulkanContext::MultiplyMatrix(cur, mat);
}

inline void vkMultMatrixf(const VKfloat* m)
{
  G4VulkanContext::Matrix4d mat{};
  for (int i = 0; i < 16; ++i) mat[i] = static_cast<double>(m[i]);
  auto& cur = G4VulkanContext::Instance().CurrentMatrix();
  cur = G4VulkanContext::MultiplyMatrix(cur, mat);
}

inline void vkTranslated(VKdouble x, VKdouble y, VKdouble z)
{
  G4VulkanContext::Matrix4d t = G4VulkanContext::IdentityMatrix();
  t[12] = x;
  t[13] = y;
  t[14] = z;
  vkMultMatrixd(t.data());
}

inline void vkTranslatef(VKfloat x, VKfloat y, VKfloat z) { vkTranslated(x, y, z); }

inline void vkScaled(VKdouble x, VKdouble y, VKdouble z)
{
  G4VulkanContext::Matrix4d s = G4VulkanContext::IdentityMatrix();
  s[0] = x;
  s[5] = y;
  s[10] = z;
  vkMultMatrixd(s.data());
}

inline void vkScalef(VKfloat x, VKfloat y, VKfloat z) { vkScaled(x, y, z); }

inline void vkRotated(VKdouble angleDeg, VKdouble x, VKdouble y, VKdouble z)
{
  const double len = std::sqrt(x * x + y * y + z * z);
  if (len < 1e-12) return;
  x /= len;
  y /= len;
  z /= len;
  const double rad = angleDeg * 3.14159265358979323846 / 180.0;
  const double c = std::cos(rad);
  const double s = std::sin(rad);
  const double omc = 1.0 - c;
  G4VulkanContext::Matrix4d r = G4VulkanContext::IdentityMatrix();
  r[0] = x * x * omc + c;
  r[1] = y * x * omc + z * s;
  r[2] = x * z * omc - y * s;
  r[4] = x * y * omc - z * s;
  r[5] = y * y * omc + c;
  r[6] = y * z * omc + x * s;
  r[8] = x * z * omc + y * s;
  r[9] = y * z * omc - x * s;
  r[10] = z * z * omc + c;
  vkMultMatrixd(r.data());
}

inline void vkRotatef(VKfloat angleDeg, VKfloat x, VKfloat y, VKfloat z)
{
  vkRotated(angleDeg, x, y, z);
}

inline void vkOrtho(VKdouble left, VKdouble right, VKdouble bottom, VKdouble top, VKdouble zNear,
                    VKdouble zFar)
{
  if (std::abs(right - left) < 1e-12 || std::abs(top - bottom) < 1e-12 ||
      std::abs(zFar - zNear) < 1e-12)
    return;
  G4VulkanContext::Matrix4d m = G4VulkanContext::IdentityMatrix();
  m[0] = 2.0 / (right - left);
  m[5] = 2.0 / (top - bottom);
  m[10] = -2.0 / (zFar - zNear);
  m[12] = -(right + left) / (right - left);
  m[13] = -(top + bottom) / (top - bottom);
  m[14] = -(zFar + zNear) / (zFar - zNear);
  vkMultMatrixd(m.data());
}

inline void vkFrustum(VKdouble left, VKdouble right, VKdouble bottom, VKdouble top, VKdouble zNear,
                      VKdouble zFar)
{
  if (zNear <= 0.0 || zFar <= 0.0 || std::abs(right - left) < 1e-12 ||
      std::abs(top - bottom) < 1e-12 || std::abs(zFar - zNear) < 1e-12)
    return;
  G4VulkanContext::Matrix4d m{};
  m[0] = (2.0 * zNear) / (right - left);
  m[5] = (2.0 * zNear) / (top - bottom);
  m[8] = (right + left) / (right - left);
  m[9] = (top + bottom) / (top - bottom);
  m[10] = -(zFar + zNear) / (zFar - zNear);
  m[11] = -1.0;
  m[14] = -(2.0 * zFar * zNear) / (zFar - zNear);
  vkMultMatrixd(m.data());
}

inline void vkBegin(VKenum mode)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fInBegin = true;
  ctx.fCurrentTopology = mode;
  ctx.fCurrentVertices.clear();
}

inline void vkEnd() { G4VulkanContext::Instance().EndPrimitive(); }

inline void vkVertex3d(VKdouble x, VKdouble y, VKdouble z)
{
  auto& ctx = G4VulkanContext::Instance();
  G4VulkanVertex v{};
  v.pos[0] = static_cast<float>(x);
  v.pos[1] = static_cast<float>(y);
  v.pos[2] = static_cast<float>(z);
  v.normal[0] = ctx.fCurrentNormal[0];
  v.normal[1] = ctx.fCurrentNormal[1];
  v.normal[2] = ctx.fCurrentNormal[2];
  v.color[0] = ctx.fCurrentColor[0];
  v.color[1] = ctx.fCurrentColor[1];
  v.color[2] = ctx.fCurrentColor[2];
  v.color[3] = ctx.fCurrentColor[3];
  ctx.fCurrentVertices.push_back(v);
}

inline void vkVertex3f(VKfloat x, VKfloat y, VKfloat z) { vkVertex3d(x, y, z); }
inline void vkVertex2d(VKdouble x, VKdouble y) { vkVertex3d(x, y, 0.0); }
inline void vkVertex2f(VKfloat x, VKfloat y) { vkVertex3d(x, y, 0.0); }
inline void vkVertex3dv(const VKdouble* v) { vkVertex3d(v[0], v[1], v[2]); }
inline void vkVertex3fv(const VKfloat* v) { vkVertex3d(v[0], v[1], v[2]); }

inline void vkNormal3d(VKdouble nx, VKdouble ny, VKdouble nz)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fCurrentNormal[0] = static_cast<float>(nx);
  ctx.fCurrentNormal[1] = static_cast<float>(ny);
  ctx.fCurrentNormal[2] = static_cast<float>(nz);
}

inline void vkNormal3f(VKfloat nx, VKfloat ny, VKfloat nz) { vkNormal3d(nx, ny, nz); }
inline void vkNormal3dv(const VKdouble* v) { vkNormal3d(v[0], v[1], v[2]); }
inline void vkNormal3fv(const VKfloat* v) { vkNormal3d(v[0], v[1], v[2]); }

inline void vkColor4d(VKdouble r, VKdouble g, VKdouble b, VKdouble a)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fCurrentColor[0] = static_cast<float>(r);
  ctx.fCurrentColor[1] = static_cast<float>(g);
  ctx.fCurrentColor[2] = static_cast<float>(b);
  ctx.fCurrentColor[3] = static_cast<float>(a);
}

inline void vkColor4f(VKfloat r, VKfloat g, VKfloat b, VKfloat a) { vkColor4d(r, g, b, a); }
inline void vkColor3d(VKdouble r, VKdouble g, VKdouble b) { vkColor4d(r, g, b, 1.0); }
inline void vkColor3f(VKfloat r, VKfloat g, VKfloat b) { vkColor4d(r, g, b, 1.0); }
inline void vkColor4dv(const VKdouble* v) { vkColor4d(v[0], v[1], v[2], v[3]); }
inline void vkColor4fv(const VKfloat* v) { vkColor4d(v[0], v[1], v[2], v[3]); }
inline void vkColor3dv(const VKdouble* v) { vkColor4d(v[0], v[1], v[2], 1.0); }
inline void vkColor3fv(const VKfloat* v) { vkColor4d(v[0], v[1], v[2], 1.0); }

inline void vkPointSize(VKfloat size)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fPointSize = std::max(1.0f, size);
  if (ctx.fRenderMode == VK_FEEDBACK && ctx.fFeedbackBuffer)
  {
    if (ctx.fFeedbackOffset + 2 <= ctx.fFeedbackSize)
    {
      ctx.fFeedbackBuffer[ctx.fFeedbackOffset++] = VK_PASS_THROUGH_TOKEN;
      ctx.fFeedbackBuffer[ctx.fFeedbackOffset++] = 700.0f + size;
    }
  }
}

inline void vkLineWidth(VKfloat width)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fLineWidth = std::max(1.0f, width);
}

inline void vkLineStipple(VKint factor, VKushort pattern)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fLineStippleRepeat = factor;
  ctx.fLineStipplePattern = pattern;
}

inline void vkPolygonStipple(const VKubyte*) {}
inline void vkPolygonMode(VKenum, VKenum mode)
{
  G4VulkanContext::Instance().fPolygonMode = mode;
}

inline void vkPolygonOffset(VKfloat factor, VKfloat units)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fPolyOffsetFactor = factor;
  ctx.fPolyOffsetUnits = units;
}

inline void vkShadeModel(VKenum mode) { G4VulkanContext::Instance().fShadeModel = mode; }
inline void vkFrontFace(VKenum mode) { G4VulkanContext::Instance().fFrontFace = mode; }
inline void vkCullFace(VKenum mode) { G4VulkanContext::Instance().fCullFaceMode = mode; }
inline void vkEdgeFlag(VKboolean) {}

inline void vkDepthFunc(VKenum func) { G4VulkanContext::Instance().fDepthFunc = func; }

inline void vkDepthMask(VKboolean flag) { G4VulkanContext::Instance().fDepthMask = flag; }

inline void vkDepthRange(VKclampd, VKclampd) {}

inline void vkColorMask(VKboolean r, VKboolean g, VKboolean b, VKboolean a)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fColorMask[0] = r;
  ctx.fColorMask[1] = g;
  ctx.fColorMask[2] = b;
  ctx.fColorMask[3] = a;
}

inline void vkBlendFunc(VKenum sfactor, VKenum dfactor)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fBlendSrc = sfactor;
  ctx.fBlendDst = dfactor;
}

inline void vkLightfv(VKenum light, VKenum pname, const VKfloat* params)
{
  if (light == VK_LIGHT0 && pname == VK_POSITION && params)
  {
    auto& ctx = G4VulkanContext::Instance();
    for (int i = 0; i < 4; ++i) ctx.fLight0Pos[i] = params[i];
    ctx.fCachedMvpDirty = true;
  }
}

inline void vkLightf(VKenum, VKenum, VKfloat) {}
inline void vkLightModelfv(VKenum, const VKfloat*) {}
inline void vkLightModeli(VKenum, VKint) {}
inline void vkMaterialfv(VKenum, VKenum, const VKfloat*) {}
inline void vkMaterialf(VKenum, VKenum, VKfloat) {}
inline void vkColorMaterial(VKenum, VKenum) {}
inline void vkHint(VKenum, VKenum) {}

inline void vkClipPlane(VKenum plane, const VKdouble* equation)
{
  if (!equation) return;
  const int idx = static_cast<int>(plane) - static_cast<int>(VK_CLIP_PLANE0);
  if (idx < 0 || idx >= 6) return;
  auto& ctx = G4VulkanContext::Instance();
  ctx.fClipPlanes[idx] = {equation[0], equation[1], equation[2], equation[3]};
}

// Display Lists
inline VKuint vkGenLists(VKsizei range)
{
  if (range <= 0) return 0;
  auto& ctx = G4VulkanContext::Instance();
  VKuint base = ctx.fNextListId;
  ctx.fNextListId += static_cast<VKuint>(range);
  for (VKsizei i = 0; i < range; ++i)
  {
    ctx.fDisplayLists[base + static_cast<VKuint>(i)] = {};
  }
  return base;
}

inline void vkNewList(VKuint list, VKenum mode)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.EnsureVulkanInitialised();
  ctx.fCurrentCompileList = list;
  ctx.fCompileMode = mode;
  ctx.fDisplayLists[list].Clear();
}

inline void vkEndList()
{
  auto& ctx = G4VulkanContext::Instance();
  if (ctx.fCurrentCompileList != 0)
  {
    auto it = ctx.fDisplayLists.find(ctx.fCurrentCompileList);
    if (it != ctx.fDisplayLists.end())
    {
      ctx.UploadDisplayListToGpu(it->second);
    }
    ctx.fCurrentCompileList = 0;
  }
}

inline void vkCallList(VKuint list)
{
  auto& ctx = G4VulkanContext::Instance();
  auto it = ctx.fDisplayLists.find(list);
  if (it != ctx.fDisplayLists.end())
  {
    ctx.ExecuteDisplayList(it->second);
  }
}

inline void vkListBase(VKuint base) { G4VulkanContext::Instance().fListBase = base; }

inline void vkCallLists(VKsizei n, VKenum type, const VKvoid* lists)
{
  if (!lists || n <= 0) return;
  const VKuint base = G4VulkanContext::Instance().fListBase;
  if (type == VK_UNSIGNED_BYTE)
  {
    const auto* p = static_cast<const VKubyte*>(lists);
    for (VKsizei i = 0; i < n; ++i) vkCallList(base + p[i]);
  }
}

inline void vkDeleteLists(VKuint list, VKsizei range)
{
  auto& ctx = G4VulkanContext::Instance();
  for (VKsizei i = 0; i < range; ++i)
  {
    ctx.fDisplayLists.erase(list + static_cast<VKuint>(i));
  }
}

inline VKboolean vkIsList(VKuint list)
{
  return G4VulkanContext::Instance().fDisplayLists.count(list) ? VK_TRUE : VK_FALSE;
}

inline void vkPushAttrib(VKbitfield) {}
inline void vkPopAttrib() {}

inline void vkRasterPos3d(VKdouble x, VKdouble y, VKdouble z)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fRasterPos[0] = x;
  ctx.fRasterPos[1] = y;
  ctx.fRasterPos[2] = z;
  if (ctx.fRenderMode == VK_FEEDBACK && ctx.fFeedbackBuffer)
  {
    G4VulkanVertex v{};
    v.pos[0] = static_cast<float>(x);
    v.pos[1] = static_cast<float>(y);
    v.pos[2] = static_cast<float>(z);
    for (int i = 0; i < 4; ++i) v.color[i] = ctx.fCurrentColor[i];
    auto pv = ctx.TransformAndProject(v);
    if (!pv.clipped && ctx.fFeedbackOffset < ctx.fFeedbackSize)
    {
      ctx.fFeedbackBuffer[ctx.fFeedbackOffset++] = VK_BITMAP_TOKEN;
      ctx.AppendFeedbackVertex(pv);
    }
  }
}

#ifndef VK_LUMINANCE
#  define VK_LUMINANCE 0x1909
#  define VK_UNPACK_SWAP_BYTES 0x0CF0
#  define VK_UNPACK_LSB_FIRST 0x0CF1
#  define VK_UNPACK_ROW_LENGTH 0x0CF2
#  define VK_UNPACK_SKIP_ROWS 0x0CF3
#  define VK_UNPACK_SKIP_PIXELS 0x0CF4
#  define VK_PACK_SWAP_BYTES 0x0D00
#  define VK_PACK_LSB_FIRST 0x0D01
#  define VK_PACK_ROW_LENGTH 0x0D02
#  define VK_PACK_SKIP_ROWS 0x0D03
#  define VK_PACK_SKIP_PIXELS 0x0D04
#endif

inline void vkBitmap(VKsizei, VKsizei, VKfloat, VKfloat, VKfloat, VKfloat, const VKubyte*) {}
inline void vkPixelStorei(VKenum, VKint) {}
inline void vkReadBuffer(VKenum) {}

inline void vkReadPixels(VKint x, VKint y, VKsizei width, VKsizei height, VKenum format,
                         VKenum type, VKvoid* pixels)
{
  if (!pixels || width <= 0 || height <= 0) return;
  auto& ctx = G4VulkanContext::Instance();
  ctx.ResizeFramebuffer(ctx.fViewport[2], ctx.fViewport[3]);

  if (type == VK_UNSIGNED_BYTE &&
      (format == VK_RGB || format == VK_RGBA || format == VK_LUMINANCE))
  {
    const int channels = (format == VK_RGBA) ? 4 : (format == VK_RGB ? 3 : 1);
    auto* out = static_cast<uint8_t*>(pixels);

    if (ctx.fDevice != VK_NULL_HANDLE)
    {
      std::vector<uint8_t> gpuRgba(static_cast<std::size_t>(width) *
                                   static_cast<std::size_t>(height) * 4,
                                   0);
      if (ctx.ReadGpuPixels(x, y, width, height, gpuRgba.data()))
      {
        if (channels == 4)
        {
          std::memcpy(out, gpuRgba.data(), gpuRgba.size());
          return;
        }
        const std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
        for (std::size_t i = 0; i < count; ++i)
        {
          const uint8_t r = gpuRgba[i * 4 + 0];
          const uint8_t g = gpuRgba[i * 4 + 1];
          const uint8_t b = gpuRgba[i * 4 + 2];
          if (channels == 1)
          {
            out[i] = static_cast<uint8_t>(0.299 * r + 0.587 * g + 0.114 * b);
          }
          else
          {
            out[i * 3 + 0] = r;
            out[i * 3 + 1] = g;
            out[i * 3 + 2] = b;
          }
        }
        return;
      }
    }

    for (VKsizei row = 0; row < height; ++row)
    {
      const int sy = y + row;
      for (VKsizei col = 0; col < width; ++col)
      {
        const int sx = x + col;
        uint8_t r = 0, g = 0, b = 0, a = 255;
        if (sx >= 0 && sx < ctx.fFbWidth && sy >= 0 && sy < ctx.fFbHeight)
        {
          const std::size_t idx =
            (static_cast<std::size_t>(sy) * static_cast<std::size_t>(ctx.fFbWidth) +
             static_cast<std::size_t>(sx)) *
            4;
          r = ctx.fColorBuffer[idx + 0];
          g = ctx.fColorBuffer[idx + 1];
          b = ctx.fColorBuffer[idx + 2];
          a = ctx.fColorBuffer[idx + 3];
        }
        const std::size_t dstIdx =
          (static_cast<std::size_t>(row) * static_cast<std::size_t>(width) +
           static_cast<std::size_t>(col)) *
          static_cast<std::size_t>(channels);
        if (channels == 1)
        {
          out[dstIdx] = static_cast<uint8_t>(0.299 * r + 0.587 * g + 0.114 * b);
        }
        else
        {
          out[dstIdx + 0] = r;
          out[dstIdx + 1] = g;
          out[dstIdx + 2] = b;
          if (channels == 4) out[dstIdx + 3] = a;
        }
      }
    }
  }
}

// Selection & Feedback
inline void vkSelectBuffer(VKsizei size, VKuint* buffer)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fSelectBuffer = buffer;
  ctx.fSelectBufferSize = size;
  ctx.fSelectOffset = 0;
  ctx.fSelectHits = 0;
  ctx.fLastHitRecordOffset = -1;
}

inline void vkFeedbackBuffer(VKsizei size, VKenum type, VKfloat* buffer)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fFeedbackBuffer = buffer;
  ctx.fFeedbackSize = size;
  ctx.fFeedbackType = type;
  ctx.fFeedbackOffset = 0;
}

inline VKint vkRenderMode(VKenum mode)
{
  auto& ctx = G4VulkanContext::Instance();
  const VKenum prev = ctx.fRenderMode;
  ctx.fRenderMode = mode;
  if (mode == VK_SELECT)
  {
    ctx.fSelectOffset = 0;
    ctx.fSelectHits = 0;
    ctx.fLastHitRecordOffset = -1;
    ctx.fNameStack.clear();
    return 0;
  }
  if (mode == VK_FEEDBACK)
  {
    ctx.fFeedbackOffset = 0;
    return 0;
  }
  if (prev == VK_SELECT) return ctx.fSelectHits;
  if (prev == VK_FEEDBACK) return ctx.fFeedbackOffset;
  return 0;
}

inline void vkInitNames()
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fNameStack.clear();
  ctx.fLastHitRecordOffset = -1;
}

inline void vkPushName(VKuint name)
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.fNameStack.push_back(name);
  ctx.fLastHitRecordOffset = -1;
}

inline void vkPopName()
{
  auto& ctx = G4VulkanContext::Instance();
  if (!ctx.fNameStack.empty()) ctx.fNameStack.pop_back();
  ctx.fLastHitRecordOffset = -1;
}

inline void vkLoadName(VKuint name)
{
  auto& ctx = G4VulkanContext::Instance();
  if (ctx.fNameStack.empty())
    ctx.fNameStack.push_back(name);
  else
    ctx.fNameStack.back() = name;
  ctx.fLastHitRecordOffset = -1;
}

inline void vkPassThrough(VKfloat token)
{
  auto& ctx = G4VulkanContext::Instance();
  if (ctx.fRenderMode == VK_FEEDBACK && ctx.fFeedbackBuffer)
  {
    if (ctx.fFeedbackOffset + 2 <= ctx.fFeedbackSize)
    {
      ctx.fFeedbackBuffer[ctx.fFeedbackOffset++] = VK_PASS_THROUGH_TOKEN;
      ctx.fFeedbackBuffer[ctx.fFeedbackOffset++] = token;
    }
  }
}

// State queries
inline void vkGetDoublev(VKenum pname, VKdouble* params)
{
  if (!params) return;
  auto& ctx = G4VulkanContext::Instance();
  switch (pname)
  {
    case VK_MODELVIEW_MATRIX:
      std::copy(ctx.fModelViewStack.back().begin(), ctx.fModelViewStack.back().end(), params);
      break;
    case VK_PROJECTION_MATRIX:
      std::copy(ctx.fProjectionStack.back().begin(), ctx.fProjectionStack.back().end(), params);
      break;
    case VK_CURRENT_COLOR:
      for (int i = 0; i < 4; ++i) params[i] = ctx.fCurrentColor[i];
      break;
    case VK_COLOR_CLEAR_VALUE:
      for (int i = 0; i < 4; ++i) params[i] = ctx.fClearColor[i];
      break;
    case VK_VIEWPORT:
      for (int i = 0; i < 4; ++i) params[i] = ctx.fViewport[i];
      break;
    default:
      params[0] = 0.0;
      break;
  }
}

inline void vkGetFloatv(VKenum pname, VKfloat* params)
{
  if (!params) return;
  auto& ctx = G4VulkanContext::Instance();
  switch (pname)
  {
    case VK_MODELVIEW_MATRIX:
      for (int i = 0; i < 16; ++i) params[i] = static_cast<float>(ctx.fModelViewStack.back()[i]);
      break;
    case VK_PROJECTION_MATRIX:
      for (int i = 0; i < 16; ++i) params[i] = static_cast<float>(ctx.fProjectionStack.back()[i]);
      break;
    case VK_CURRENT_COLOR:
      for (int i = 0; i < 4; ++i) params[i] = ctx.fCurrentColor[i];
      break;
    case VK_COLOR_CLEAR_VALUE:
      for (int i = 0; i < 4; ++i) params[i] = ctx.fClearColor[i];
      break;
    case VK_POINT_SIZE:
      params[0] = ctx.fPointSize;
      break;
    case VK_LINE_WIDTH:
      params[0] = ctx.fLineWidth;
      break;
    case VK_VIEWPORT:
      for (int i = 0; i < 4; ++i) params[i] = static_cast<float>(ctx.fViewport[i]);
      break;
    default:
      params[0] = 0.0f;
      break;
  }
}

inline void vkGetIntegerv(VKenum pname, VKint* params)
{
  if (!params) return;
  auto& ctx = G4VulkanContext::Instance();
  switch (pname)
  {
    case VK_VIEWPORT:
      for (int i = 0; i < 4; ++i) params[i] = ctx.fViewport[i];
      break;
    case VK_RENDER_MODE:
      params[0] = static_cast<VKint>(ctx.fRenderMode);
      break;
    case VK_MAX_CLIP_PLANES:
      params[0] = 6;
      break;
    case VK_MAX_LIGHTS:
      params[0] = 8;
      break;
    case VK_MAX_VIEWPORT_DIMS:
      params[0] = 8192;
      params[1] = 8192;
      break;
    case VK_BLEND_SRC:
      params[0] = static_cast<VKint>(ctx.fBlendSrc);
      break;
    case VK_BLEND_DST:
      params[0] = static_cast<VKint>(ctx.fBlendDst);
      break;
    case VK_LINE_STIPPLE_PATTERN:
      params[0] = static_cast<VKint>(ctx.fLineStipplePattern);
      break;
    case VK_LINE_STIPPLE_REPEAT:
      params[0] = ctx.fLineStippleRepeat;
      break;
    default:
      params[0] = 0;
      break;
  }
}

inline void vkGetBooleanv(VKenum pname, VKboolean* params)
{
  if (!params) return;
  params[0] = vkIsEnabled(pname);
}

inline VKenum vkGetError() { return VK_NO_ERROR; }

inline void vkFlush()
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.EnsureVulkanInitialised();
  if (ctx.fDevice != VK_NULL_HANDLE)
  {
    ctx.FlushGpuCommands();
  }
}

inline void vkFinish()
{
  auto& ctx = G4VulkanContext::Instance();
  ctx.EnsureVulkanInitialised();
  if (ctx.fDevice != VK_NULL_HANDLE)
  {
    ctx.FlushGpuCommands();
    ctx.WaitGpuIdle();
  }
}

#ifndef VK_KEEP
#  define VK_KEEP 0x1E00
#  define VK_REPLACE 0x1E01
#  define VK_INCR 0x1E02
#  define VK_DECR 0x1E03
#  define VK_INVERT 0x150A
#endif

inline void vkStencilFunc(VKenum, VKint, VKuint) {}
inline void vkStencilOp(VKenum, VKenum, VKenum) {}
inline void vkStencilMask(VKuint) {}
inline void vkClearStencil(VKint) {}

#define G4VULKAN_FLT_BIG 1.e20

#endif
