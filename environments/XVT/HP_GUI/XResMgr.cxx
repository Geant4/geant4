// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: XResMgr.cxx,v 1.2 1999-12-15 14:48:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/********************************************************************

   Copyright (c) 1995-96 XVT Software, Inc.

   File Name:           ResMgr.cxx

   Usage:

     If you plan to use Motif-specific icons or cursors in your application,
     you must edit this file and add the appropriate definitions.  Using
     XVT-Architect's Drafting Board module, you can lay out CIcon objects
     in a project's windows.  However, the resources corresponding to CIcon
     objects cannot be specified via XVT's Universal Resource Language
     (URL).  You must instead define them in this file with a small amount
     of C code.

     For more information on the purpose and contents of this file, see
     the "Cursors and Drawn Icons" section in the XVT Platform-Specific
     Book for Motif.

 ********************************************************************/

#include "XVTPwr.h"
#include "xvt_xres.h"
#include Global_i		// for NULLicon

#define pwr_width 32
#define pwr_height 30
static unsigned char pwr_bits[] =
{
	0xff, 0xff, 0xff, 0xff, 0x01, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80,
	0xf9, 0x03, 0x07, 0x87, 0x09, 0x04, 0x05, 0x85, 0x09, 0x08, 0x05, 0x85,
	0x09, 0xe8, 0xfd, 0xbd, 0x09, 0x28, 0x20, 0xa0, 0x09, 0xe8, 0xfd, 0xbd,
	0x09, 0x08, 0x05, 0x85, 0x09, 0x08, 0x05, 0x85, 0x09, 0x08, 0x07, 0x87,
	0x09, 0x04, 0x00, 0x80, 0xf9, 0x03, 0x00, 0x80, 0x09, 0x00, 0x00, 0x80,
	0x09, 0x4e, 0xf4, 0x9e, 0x09, 0x51, 0x14, 0xa2, 0x09, 0x51, 0x14, 0xa2,
	0x09, 0x51, 0x14, 0xa2, 0x09, 0x51, 0x14, 0xa2, 0x09, 0x51, 0x74, 0x9e,
	0x09, 0x51, 0x15, 0x8a, 0x09, 0x51, 0x15, 0x92, 0x09, 0x51, 0x15, 0x92,
	0x09, 0x91, 0x12, 0xa2, 0x09, 0x8e, 0xf2, 0xa2, 0x01, 0x00, 0x00, 0x80,
	0x01, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80, 0xff, 0xff, 0xff, 0xff
};

// Define an icon or cursor resource type variable:
static ICON_RESOURCE PwrIcon;

// Place resouce into XVT's resource information table:
RESOURCE_INFO rtable[] =
{
	{ "ICON", NULLicon, (char *)&PwrIcon },
	{ 0 }
};

// Build the resource:
RESOURCE_INFO* xvt_xres_create_table()
{
	xvt_xres_build_icon( &PwrIcon, NULLheight, NULLwidth, (DATA_PTR)pwr_bits );
	return rtable;
}
