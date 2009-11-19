#ifndef Histograms_h
#define Histograms_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              Histograms.hh
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
// This module defines the HistSpecialBin type which can have values
// underflow_bin (==0), overflow_bin (==1), or inrange(==2).
//
// C-preprocessor commands are also used to replace expressions BIN_UNDERFLOW
// and BIN_OVERFLOW with numerical values (-1 and-2 respectively).
// BIN_UNDERFLOW and BIN_OVERFLOW are used to denote underflow and overflow bin
// conditions when the return type required is integer rather than a
// HistSpecialBin type.
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

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
enum HistSpecialBin {underflow_bin=0, overflow_bin=1, inrange=2};
#define BIN_UNDERFLOW -1
#define BIN_OVERFLOW -2
////////////////////////////////////////////////////////////////////////////////
#endif


