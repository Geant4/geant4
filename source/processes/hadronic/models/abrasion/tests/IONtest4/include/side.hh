#ifndef side_h
#define side_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              side.hh
//
// Version:		1.0
// Date:		09/03/00
// Author:		P R Truscott
// Organisation:	DERA UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		12115/96/NL/JG Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------
//
// This module defines the side type which can have values LEFT or RIGHT.
// It is used to define whether the datapoints in a histogram that fall
// exactly on a bin-edge belong to the bin to the left of the edge or to the
// right.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
//
// (None)
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 30 June 1999, P R Truscott, DERA UK
// Version number update 0.b.2 -> 0.b.3, but no functional change.
//
// 28 August 1999, F Lei, DERA UK
// Version number update 0.b.3 -> 0.b.4, but no functional change.
//
// 17 September 1999, P R Truscott, DERA UK
// Version number update 0.b.4 -> 0.b.5, but no functional change.
//
// 09 March 2000, P R Truscott, DERA UK
// Update 0.b.3 -> 1.0, for compliance with ISO ANSI C++ (no functional change).
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
enum side{LEFT,RIGHT};
////////////////////////////////////////////////////////////////////////////////
#endif




