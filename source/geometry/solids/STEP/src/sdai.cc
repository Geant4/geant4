
// Use limits.h instead of values.h for portability. Therefore,
// use FLT_MIN instead of MINFLOAT, LONG_MAX instead of MAXLONG.
// Also commented inclusion of values.h in sdai.h - GC
#include <limits.h>
#include <sdai.h>
#include "globals.hh"
//  The values used here to represent NULL values for numeric types come from 
//  the ANSI C standard

const char *SCLversion = "STEP Class Library version 3.1";

const SCLP23(Integer) SCLP23(INT_NULL) = LONG_MAX;
const SCLP23(Real)    SCLP23(REAL_NULL) = FLT_MIN;
const SCLP23(Real)    SCLP23(NUMBER_NULL) = FLT_MIN;
