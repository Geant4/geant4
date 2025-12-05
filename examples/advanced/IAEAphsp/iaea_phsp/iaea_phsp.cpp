/*
 * Copyright (C) 2006 International Atomic Energy Agency
 * -----------------------------------------------------------------------------
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is furnished
 * to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 *-----------------------------------------------------------------------------
 *
 *   AUTHORS:
 *
 *   Roberto Capote Noy, PhD
 *   e-mail: R.CapoteNoy@iaea.org (rcapotenoy@yahoo.com)
 *   International Atomic Energy Agency
 *   Nuclear Data Section, P.O.Box 100
 *   Wagramerstrasse 5, Vienna A-1400, AUSTRIA
 *   Phone: +431-260021713; Fax: +431-26007
 *
 *   Iwan Kawrakow, PhD
 *   e-mail iwan@irs.phy.nrc.ca
 *   Ionizing Radiation Standards
 *   Institute for National Measurement Standards
 *   National Research Council of Canada Ottawa, ON, K1A 0R6 Canada
 *   Phone: +1-613-993 2197, ext.241; Fax: +1-613-952 9865
 *
 **********************************************************************************
 * For documentation
 * see http://www-nds.iaea.org/reports-new/indc-reports/indc-nds/indc-nds-0484.pdf
 **********************************************************************************/
//
#define MAIN

#if (defined WIN32) || (defined WIN64)
#include <iostream>  // so that namespace std becomes defined
#endif
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cctype>
#include <algorithm>  // for max() and min()

#include<sys/types.h>
#include<sys/stat.h>

#if !(defined WIN32) && !(defined WIN64)
using namespace std;
#endif

#include "utilities.h"
#include "iaea_record.h"
#include "iaea_header.h"
#include "iaea_phsp.h"

//#define false 0  //(MACG) to silent -Wkeyword-macro build warning
//#define true  1  //(MACG) to silent -Wkeyword-macro build warning

// These variables are defined globally. They contain pointers
// to header and record structures defined by calling iaea_new_source()
// routine to maintain a list of already initialized IAEA sources.

static iaea_header_type *p_iaea_header[MAX_NUM_SOURCES];
static iaea_record_type *p_iaea_record[MAX_NUM_SOURCES];

/************************************************************************
* Initialization
*
* Given a file name header_file of length hf_length, initialize a
* new IAEA particle source, assign an unique Id to it and return this
* Id in result. Dont assume header_file is null-terminated as the
* function may be called from a Fortran program.
* The need for an Id arises from the fact that some applications may
* want to use several IAEA sources at once. The implementation must therefore
* maintain a list of already initialized sources.
* If an error occures (e.g. header file does not exist, there are errors
* in the header file, etc.), assign a negative number to result.
* (one may want to specify a list of error codes so that the application
* knows what went wrong). This function *must* be called before using
* any of the following functions for a given source id.
*
* access = 1 => opening read-only file
* access = 2 => opening file for writing
* access = 3 => opening file for appending/updating
*
***********************************************************************/

static int __iaea_source_used[MAX_NUM_SOURCES];
static int __iaea_n_source = 0;

IAEA_EXTERN_C IAEA_EXPORT
void iaea_new_source(IAEA_I32 *source_ID, char *header_file,
                     const IAEA_I32 *access, IAEA_I32 *result,
                     int hf_length) {

   if( !header_file ) {
       *result = 105; *source_ID = -1; return;
   } // null header file name
   if(*access != 1 && *access != 2 && *access != 3) {
       *result = -99 ; *source_ID = -1; return;
   } // Wrong access requested

   if( hf_length >= MAX_STR_LEN) {
       *result = -100 ; *source_ID = -1; return;
   } // Too long string
   if( hf_length < 1)            {
       *result = -101 ; *source_ID = -1; return;
   } // String length < 1

   if( !__iaea_n_source ) { // called for the first time
       for(int j=0; j<MAX_NUM_SOURCES; j++) __iaea_source_used[j] = false;
   }
   int sid=-1;
   // do we have a spare spot in the arrays ?
   // (e.g. because a source was destroyed)
   for(int j=0; j<__iaea_n_source; j++) {
       if( !__iaea_source_used[j] ) { sid = j; break; }
   }
   if( sid < 0 ) {
       // so, we don't => increase source count and check if
       // space left in arrays.
       if( ++__iaea_n_source >= MAX_NUM_SOURCES ) {
           *result = -98; *source_ID = -1; return;
       }
       sid = __iaea_n_source-1; *source_ID = sid;
   }
   __iaea_source_used[sid] = true;

   //int ilen = strlen(header_file);
   // the above requires a null-terminated string, which may not be
   // satisfied, if we are called from a fortran progrmam.

   int ilen = hf_length;
   if( header_file[hf_length-1] != 0 ) { // not null-terminated
       while(ilen > 0 && isspace(header_file[--ilen]) );
       if( ilen < hf_length-1 ) header_file[ilen+1] = '\0';
   }

   // Creating IAEA phsp header and allocating memory for it
   p_iaea_header[*source_ID] = (iaea_header_type *) calloc(1, sizeof(iaea_header_type));
   // Opening header file
   if(*access == 1) p_iaea_header[*source_ID]->fheader =
         open_file(header_file,".IAEAheader","rb");
   if(*access == 2) p_iaea_header[*source_ID]->fheader =
         open_file(header_file,".IAEAheader","wb");
   if(*access == 3) p_iaea_header[*source_ID]->fheader =
         open_file(header_file,".IAEAheader","r+b");

   if(p_iaea_header[*source_ID]->fheader == NULL) { *result = -96; return;} // phsp failed to open

   // Creating IAEA record and allocating memory for it
   p_iaea_record[*source_ID] = (iaea_record_type *) calloc(1, sizeof(iaea_record_type));

   p_iaea_header[*source_ID]->initialize_counters();

   switch( *access )
   {
         case 2: // writing a new phsp

             strcpy(p_iaea_header[*source_ID]->title,"PHASESPACE in IAEA format");
             // Default IAEA index
             *result = p_iaea_header[*source_ID]->iaea_index = 1000;

             p_iaea_record[*source_ID]->p_file =
                 open_file(header_file, ".IAEAphsp", "wb");

             if(p_iaea_record[*source_ID]->p_file == NULL) { *result = -94 ; return; }

             // Setting default i/o flags
             if(p_iaea_record[*source_ID]->initialize() != OK)
                 {*result = -1; return;}

             if( p_iaea_header[*source_ID]->set_record_contents(p_iaea_record[*source_ID])
                 == FAIL ) { *result = -95; return;}

             return;

         case 3 : // appending to the existing phsp

             if( p_iaea_header[*source_ID]->read_header() != OK)
                 { *result = -93; return;}

             int i;
             // Setting up Average Kinetic Energy counters to usable values
             for(i=0;i<MAX_NUM_PARTICLES;i++)
                 p_iaea_header[*source_ID]->averageKineticEnergy[i] *=
                 p_iaea_header[*source_ID]->sumParticleWeight[i];

             // Opening phsp file to append
             p_iaea_record[*source_ID]->p_file =
                 open_file(header_file, ".IAEAphsp", "a+b");

             if(p_iaea_record[*source_ID]->p_file == NULL) { *result = -94 ; return; }

             if(p_iaea_record[*source_ID]->initialize() != OK) {*result = -1; return;}

             // Get read/write logical block from the header
             if( p_iaea_header[*source_ID]->get_record_contents(p_iaea_record[*source_ID])
                 == FAIL) { *result = -91; return;}

             *result = p_iaea_header[*source_ID]->iaea_index; // returning IAEA index

             break;

         case 1 : // reading existing phsp

             if( p_iaea_header[*source_ID]->read_header() != OK) { *result = -93; return;}

             // Opening phsp file to read
             p_iaea_record[*source_ID]->p_file =
                 open_file(header_file, ".IAEAphsp", "rb");

             if(p_iaea_record[*source_ID]->p_file == NULL)
                 { *result = -94 ; return; }

             if(p_iaea_record[*source_ID]->initialize() != OK) {*result = -1; return;}

             // Get read/write logical block from the header
             if( p_iaea_header[*source_ID]->get_record_contents(p_iaea_record[*source_ID])
                 == FAIL) { *result = -91; return;}

             *result = p_iaea_header[*source_ID]->iaea_index; // returning IAEA index

             break;
   }

   return;
}

IAEA_EXTERN_C IAEA_EXPORT
void iaea_new_source_(IAEA_I32 *source_ID, char *header_file,
                      const IAEA_I32 *access, IAEA_I32 *result,
                      int hf_length) {
    iaea_new_source(source_ID,header_file,access,result,hf_length);
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_new_source__(IAEA_I32 *source_ID, char *header_file,
                      const IAEA_I32 *access, IAEA_I32 *result,
                      IAEA_I32 hf_length) {
    iaea_new_source(source_ID,header_file,access,result,hf_length);
}
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_NEW_SOURCE(IAEA_I32 *source_ID, char *header_file,
                      const IAEA_I32 *access, IAEA_I32 *result,
                      IAEA_I32 hf_length) {
    iaea_new_source(source_ID,header_file,access,result,hf_length);
}
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_NEW_SOURCE_(IAEA_I32 *source_ID, char *header_file,
                      const IAEA_I32 *access, IAEA_I32 *result,
                      IAEA_I32 hf_length) {
    iaea_new_source(source_ID,header_file,access,result,hf_length);
}
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_NEW_SOURCE__(IAEA_I32 *source_ID, char *header_file,
                      const IAEA_I32 *access, IAEA_I32 *result,
                      IAEA_I32 hf_length) {
    iaea_new_source(source_ID,header_file,access,result,hf_length);
}

/************************************************************************
* Maximum number of particles
*
* Set n_particle to the maximum number of particle of type type the
* source with Id id can return. If type<0, set n_particle to the
* total number of all particles. For event generators, this function
* should set n_particle to the maximum integer that can be stored in a
* signed 64 bit integer. If the source with Id id does not exist,
* set n_particle to a negative number.
*************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_max_particles(const IAEA_I32 *id, const IAEA_I32 *type,
                                  IAEA_I64 *n_particle)
{
      // No header found
      if(p_iaea_header[*id]->fheader == NULL) {*n_particle = -1; return;}

      int file_type = p_iaea_header[*id]->file_type;

      if(file_type == 1) // Event generator
      {
	    if( sizeof(IAEA_I64) >= 8 ) {
                #if (defined WIN32) || (defined WIN64)
                *n_particle = 9223372036854775807;   // 2^63 - 1
                #else
                *n_particle = 9223372036854775807LL; // 2^63 - 1
                #endif
	    }
	    else *n_particle = 2147483647;           // 2^31 - 1
            return;
      }
      // phsp file
      if (*type < 0) {*n_particle = p_iaea_header[*id]->nParticles; return;}
      if ( (*type >= MAX_NUM_PARTICLES) || (*type ==0) ) {*n_particle = 0; return;}

      *n_particle = p_iaea_header[*id]->particle_number[*type-1];

      return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_max_particles_(const IAEA_I32 *id,
                                           const IAEA_I32 *type,
                                                 IAEA_I64 *n_particle)
{ iaea_get_max_particles(id, type, n_particle); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_max_particles__(const IAEA_I32 *id,
                                           const IAEA_I32 *type,
                                                 IAEA_I64 *n_particle)
{ iaea_get_max_particles(id, type, n_particle); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_MAX_PARTICLES(const IAEA_I32 *id,
                                           const IAEA_I32 *type,
                                                 IAEA_I64 *n_particle)
{ iaea_get_max_particles(id, type, n_particle); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_MAX_PARTICLES_(const IAEA_I32 *id,
                                           const IAEA_I32 *type,
                                                 IAEA_I64 *n_particle)
{ iaea_get_max_particles(id, type, n_particle); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_MAX_PARTICLES__(const IAEA_I32 *id,
                                           const IAEA_I32 *type,
                                                 IAEA_I64 *n_particle)
{ iaea_get_max_particles(id, type, n_particle); }

/************************************************************************
* Maximum energy
*
* Return the maximum energy of an initialized IAEA source with Id id
* in Emax. Set Emax to negative if a source with that Id does not exist.
************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_maximum_energy(const IAEA_I32 *id, IAEA_Float *Emax)
{
      // No header found
      if(p_iaea_header[*id]->fheader == NULL) {*Emax = -1.f; return;}

      int file_type = p_iaea_header[*id]->file_type;

      if(file_type == 1) {*Emax = -1.f; return;}// Event generator

      // phsp file
      *Emax = 0.f;
      for(int i=0;i<MAX_NUM_PARTICLES;i++)
      {
        // *Emax = (IAEA_Float) max(*Emax,p_iaea_header[*id]->maximumKineticEnergy[i]);
	//(MACG) using std::max from algorithm libraries
        *Emax = (IAEA_Float) max(*Emax,
				 (IAEA_Float)p_iaea_header[*id]->maximumKineticEnergy[i]);
      }
      return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_maximum_energy_(const IAEA_I32 *id, IAEA_Float *Emax)
{ iaea_get_maximum_energy(id, Emax); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_maximum_energy__(const IAEA_I32 *id, IAEA_Float *Emax)
{ iaea_get_maximum_energy(id, Emax); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_MAXIMUM_ENERGY(const IAEA_I32 *id, IAEA_Float *Emax)
{ iaea_get_maximum_energy(id, Emax); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_MAXIMUM_ENERGY_(const IAEA_I32 *id, IAEA_Float *Emax)
{ iaea_get_maximum_energy(id, Emax); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_MAXIMUM_ENERGY__(const IAEA_I32 *id, IAEA_Float *Emax)
{ iaea_get_maximum_energy(id, Emax); }

/*************************************************************************
* Number of additional floats and integers returned by the source
*
* Return the number of additional floats in n_extra_float and the number
* of additional integers in n_extra_integer for the source with Id id.
* Set n_extra_integer and/or n_extra_float to be negative if such a
* source does not exist.
*************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_extra_numbers(const IAEA_I32 *id, IAEA_I32 *n_extra_float,
                                          IAEA_I32 *n_extra_int)
{
      // No header found
      if(p_iaea_header[*id]->fheader == NULL)
            {*n_extra_float = *n_extra_int = -1; return;}

      *n_extra_float = p_iaea_header[*id]->record_contents[7];
      *n_extra_int   = p_iaea_header[*id]->record_contents[8];
      return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_extra_numbers_(const IAEA_I32 *id,
                                                 IAEA_I32 *n_extra_float,
                                                 IAEA_I32 *n_extra_int)
{ iaea_get_extra_numbers(id, n_extra_float,n_extra_int); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_extra_numbers__(const IAEA_I32 *id,
                                                 IAEA_I32 *n_extra_float,
                                                 IAEA_I32 *n_extra_int)
{ iaea_get_extra_numbers(id, n_extra_float,n_extra_int); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_EXTRA_NUMBERS(const IAEA_I32 *id,
                                                 IAEA_I32 *n_extra_float,
                                                 IAEA_I32 *n_extra_int)
{ iaea_get_extra_numbers(id, n_extra_float,n_extra_int); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_EXTRA_NUMBERS_(const IAEA_I32 *id,
                                                 IAEA_I32 *n_extra_float,
                                                 IAEA_I32 *n_extra_int)
{ iaea_get_extra_numbers(id, n_extra_float,n_extra_int); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_EXTRA_NUMBERS__(const IAEA_I32 *id,
                                                 IAEA_I32 *n_extra_float,
                                                 IAEA_I32 *n_extra_int)
{ iaea_get_extra_numbers(id, n_extra_float,n_extra_int); }

/*************************************************************************
* Number of additional floats and integers to be stored
*
* Set the number of additional floats in n_extra_float and the number
* of additional integers in n_extra_integer for the source with Id id
* to be stored in the corresponding file.
*************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_extra_numbers(const IAEA_I32 *id, IAEA_I32 *n_extra_float,
                                  IAEA_I32 *n_extra_int)
{
      // No header found
      if(p_iaea_header[*id]->fheader == NULL) return;

      p_iaea_header[*id]->record_contents[7] = *n_extra_float;
      p_iaea_header[*id]->record_contents[8] = *n_extra_int;

    // Store read/write logical block in the PHSP header
    if( p_iaea_header[*id]->get_record_contents(p_iaea_record[*id])
         == FAIL) return;

      return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_extra_numbers_(const IAEA_I32 *id,
                                                 IAEA_I32 *n_extra_float,
                                                 IAEA_I32 *n_extra_int)
{ iaea_set_extra_numbers(id, n_extra_float, n_extra_int); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_extra_numbers__(const IAEA_I32 *id,
                                                 IAEA_I32 *n_extra_float,
                                                 IAEA_I32 *n_extra_int)
{ iaea_set_extra_numbers(id, n_extra_float, n_extra_int); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_EXTRA_NUMBERS(const IAEA_I32 *id,
                                                 IAEA_I32 *n_extra_float,
                                                 IAEA_I32 *n_extra_int)
{ iaea_set_extra_numbers(id, n_extra_float, n_extra_int); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_EXTRA_NUMBERS_(const IAEA_I32 *id,
                                                 IAEA_I32 *n_extra_float,
                                                 IAEA_I32 *n_extra_int)
{ iaea_set_extra_numbers(id, n_extra_float, n_extra_int); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_EXTRA_NUMBERS__(const IAEA_I32 *id,
                                                 IAEA_I32 *n_extra_float,
                                                 IAEA_I32 *n_extra_int)
{ iaea_set_extra_numbers(id, n_extra_float, n_extra_int); }


/*******************************************************************************
* Set a type type of the extra long variable corresponding to the "index" number
* for a corresponding header of the phsp "id". Index is running from zero.
*
* The current list of types for extra long variables is:
*   0: User defined generic type
*   1: Incremental history number (EGS,PENELOPE)
*      = 0 indicates a nonprimary particle event
*      > 0 indicates a primary particle. The value is equal to the number of
*          primaries particles employed to get to this history after the last
*          primary event was recorded.
*   2: LATCH (EGS)
*   3: ILB5 (PENELOPE)
*   4: ILB4 (PENELOPE)
*   5: ILB3 (PENELOPE)
*   6: ILB2 (PENELOPE)
*   7: ILB1 (PENELOPE)
*   more to be defined
*
* Usually called before writing phsp header to set the type of extra long
* variables to be stored. It must be called once for every extralong variable.
*
* type = -1 means the source's header file does not exist
*           or source was not properly initialized (call iaea_new_...)
* type = -2 means the index is out of range ( 0 <= index < NUM_EXTRA_LONG )
* type = -3 means the type to be set is out of range
*           ( 1 <= type < MAX_NUMB_EXTRALONG_TYPES )
*******************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_type_extralong_variable(const IAEA_I32 *id,
                                      const IAEA_I32 *index,
                                            IAEA_I32 *type)
{
   // No header found
   if(p_iaea_header[*id]->fheader == NULL) {*type = -1; return;}

   if((*index < 0) || (*index >= NUM_EXTRA_LONG) ) {*type = -2; return;}

   if((*type < 0) || (*type > MAX_NUMB_EXTRALONG_TYPES) ) {*type = -3; return;}

   p_iaea_header[*id]->extralong_contents[*index] = *type;

   return;
}

IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_type_extralong_variable_(const IAEA_I32 *id,
                                                     const IAEA_I32 *index,
                                                     IAEA_I32 *type)
{ iaea_set_type_extralong_variable(id, index, type); }

IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_type_extralong_variable__(const IAEA_I32 *id,
                                                     const IAEA_I32 *index,
                                                     IAEA_I32 *type)
{ iaea_set_type_extralong_variable(id, index, type); }

IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_TYPE_EXTRALONG_VARIABLE(const IAEA_I32 *id,
                                                     const IAEA_I32 *index,
                                                     IAEA_I32 *type)
{ iaea_set_type_extralong_variable(id, index, type); }

IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_TYPE_EXTRALONG_VARIABLE_(const IAEA_I32 *id,
                                                     const IAEA_I32 *index,
                                                     IAEA_I32 *type)
{ iaea_set_type_extralong_variable(id, index, type); }

IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_TYPE_EXTRALONG_VARIABLE__(const IAEA_I32 *id,
                                                     const IAEA_I32 *index,
                                                     IAEA_I32 *type)
{ iaea_set_type_extralong_variable(id, index, type); }


/********************************************************************************
* Set a type type of the extra float variable corresponding to the "index" number
* for a corresponding header of the phsp "id". Index is running from zero.
*
* The current list of types for extra float variables is:
*   0: User defined generic type
*   1: XLAST (x coord. of the last interaction)
*   2: YLAST (y coord. of the last interaction)
*   3: ZLAST (z coord. of the last interaction)
*   more to be defined
*
* Usually called before writing phsp header to set the type of extra float
* variables to be stored. It must be called once for every extra float variable.
*
* type = -1 means the source's header file does not exist
*           or source was not properly initialized (call iaea_new_...)
* type = -2 means the index is out of range ( 0 <= index < NUM_EXTRA_FLOAT )
* type = -3 means the type to be set is out of range
*           ( 1 <= type < MAX_NUMB_EXTRAFLOAT_TYPES )
*******************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_type_extrafloat_variable(const IAEA_I32 *id,
                                       const IAEA_I32 *index,
                                             IAEA_I32 *type)
{
   // No header found
   if(p_iaea_header[*id]->fheader == NULL) {*type = -1; return;}

   if((*index < 0) || (*index >= NUM_EXTRA_FLOAT) ) {*type = -2; return;}

   if((*type < 0) || (*type > MAX_NUMB_EXTRAFLOAT_TYPES) ) {*type = -3; return;}

   p_iaea_header[*id]->extrafloat_contents[*index] = *type;

   return;
}

IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_type_extrafloat_variable_(const IAEA_I32 *id,
                                                     const IAEA_I32 *index,
                                                     IAEA_I32 *type)
{ iaea_set_type_extrafloat_variable(id, index, type); }

IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_type_extrafloat_variable__(const IAEA_I32 *id,
                                                     const IAEA_I32 *index,
                                                     IAEA_I32 *type)
{ iaea_set_type_extrafloat_variable(id, index, type); }

IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_TYPE_EXTRAFLOAT_VARIABLE(const IAEA_I32 *id,
                                                     const IAEA_I32 *index,
                                                     IAEA_I32 *type)
{ iaea_set_type_extrafloat_variable(id, index, type); }

IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_TYPE_EXTRAFLOAT_VARIABLE_(const IAEA_I32 *id,
                                                     const IAEA_I32 *index,
                                                     IAEA_I32 *type)
{ iaea_set_type_extrafloat_variable(id, index, type); }

IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_TYPE_EXTRAFLOAT_VARIABLE__(const IAEA_I32 *id,
                                                     const IAEA_I32 *index,
                                                     IAEA_I32 *type)
{ iaea_set_type_extrafloat_variable(id, index, type); }


/****************************************************************************
* Get a type type of all extra variables from a header of the phsp "id".
*
* extralong_types[] AND extrafloat_types[] must have a dimension bigger than
* MAX_NUMB_EXTRALONG_TYPES and MAX_NUMB_EXTRAFLOAT_TYPES correspondingly
*
* The current list of types for extra long variables is:
*   0: User defined generic type
*   1: Incremental history number (EGS,PENELOPE)
*      = 0 indicates a nonprimary particle event
*      > 0 indicates a primary particle. The value is equal to the number of
*          primaries particles employed to get to this history after the last
*          primary event was recorded.
*   2: LATCH (EGS)
*   3: ILB5 (PENELOPE)
*   4: ILB4 (PENELOPE)
*   5: ILB3 (PENELOPE)
*   6: ILB2 (PENELOPE)
*   7: ILB1 (PENELOPE)
*   more to be defined
*
* The current list of types for extra float variables is:
*   0: User defined generic type
*   1: XLAST (x coord. of the last interaction)
*   2: YLAST (y coord. of the last interaction)
*   3: ZLAST (z coord. of the last interaction)
*   more to be defined
*
* Usually called before reading phsp header to know the type of extra long
* variables to be read. It must be called once for every extra float variable.
*
* result = -1 means the source's header file does not exist
*             or source was not properly initialized (call iaea_new_...)
*******************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_type_extra_variables(const IAEA_I32 *id, IAEA_I32 *result,
      IAEA_I32 extralong_types[], IAEA_I32 extrafloat_types[])
{
      // No header found
      if(p_iaea_header[*id]->fheader == NULL) {*result = -1; return;}

      for (int i=0;i<p_iaea_header[*id]->record_contents[8];i++ )
          extralong_types[i] = p_iaea_header[*id]->extralong_contents[i];

      for (int j=0;j<p_iaea_header[*id]->record_contents[7];j++ )
          extrafloat_types[j] = p_iaea_header[*id]->extrafloat_contents[j];

      *result = +1;
      return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_type_extra_variables_(const IAEA_I32 *id, IAEA_I32 *result,
      IAEA_I32 extralong_types[], IAEA_I32 extrafloat_types[])
{ iaea_get_type_extra_variables(id, result, extralong_types, extrafloat_types); }

IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_type_extra_variables__(const IAEA_I32 *id, IAEA_I32 *result,
      IAEA_I32 extralong_types[], IAEA_I32 extrafloat_types[])
{ iaea_get_type_extra_variables(id, result, extralong_types, extrafloat_types); }

IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_TYPE_EXTRA_VARIABLES(const IAEA_I32 *id, IAEA_I32 *result,
      IAEA_I32 extralong_types[], IAEA_I32 extrafloat_types[])
{ iaea_get_type_extra_variables(id, result, extralong_types, extrafloat_types); }

IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_TYPE_EXTRA_VARIABLES_(const IAEA_I32 *id, IAEA_I32 *result,
      IAEA_I32 extralong_types[], IAEA_I32 extrafloat_types[])
{ iaea_get_type_extra_variables(id, result, extralong_types, extrafloat_types); }

IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_TYPE_EXTRA_VARIABLES__(const IAEA_I32 *id, IAEA_I32 *result,
      IAEA_I32 extralong_types[], IAEA_I32 extrafloat_types[])
{ iaea_get_type_extra_variables(id, result, extralong_types, extrafloat_types); }

/*************************************************************************
* Set variable corresponding to the "index" number to a "constant" value
* for a corresponding header of the phsp "id". Index is running from zero.
*
*       (Usually called as needed before MC loop started)
*
*                index  =  0 1 2 3 4 5 6
*          corresponds to  x,y,z,u,v,w,wt
*
* Usually called before writing phsp files to set those variables which
* are not going to be stored. It must be called once for every variable
*
* constant = -1 means the source's header file does not exist
*               or source was not properly initialized (call iaea_new_...)
* constant = -2 means the index is out of range ( 0 <= index < 7 )
*************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_constant_variable(const IAEA_I32 *id, const IAEA_I32 *index,
                                               IAEA_Float *constant)
{
      // No header found
      if(p_iaea_header[*id]->fheader == NULL) {*constant = -1.f; return;}

      if((*index < 0) || (*index > 6) ) {*constant = -2.f; return;}

      p_iaea_header[*id]->record_contents[*index] = 0; // variable is constant
      p_iaea_header[*id]->record_constant[*index] = *constant;

      // Store read/write logical block changes in the PHSP header
      p_iaea_header[*id]->get_record_contents(p_iaea_record[*id]);

      return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_constant_variable_(const IAEA_I32 *id,
                                               const IAEA_I32 *index,
                                                     IAEA_Float *constant)
{iaea_set_constant_variable(id, index, constant); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_constant_variable__(const IAEA_I32 *id,
                                               const IAEA_I32 *index,
                                                     IAEA_Float *constant)
{iaea_set_constant_variable(id, index, constant); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_CONSTANT_VARIABLE(const IAEA_I32 *id,
                                               const IAEA_I32 *index,
                                                     IAEA_Float *constant)
{iaea_set_constant_variable(id, index, constant); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_CONSTANT_VARIABLE_(const IAEA_I32 *id,
                                               const IAEA_I32 *index,
                                                     IAEA_Float *constant)
{iaea_set_constant_variable(id, index, constant); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_CONSTANT_VARIABLE__(const IAEA_I32 *id,
                                               const IAEA_I32 *index,
                                                     IAEA_Float *constant)
{iaea_set_constant_variable(id, index, constant); }

/*************************************************************************
* Get value of constant corresponding to the "index" number
* for a corresponding header of the phsp "id". Index is running from zero.
*
*                index  =  0 1 2 3 4 5 6
*          corresponds to  x,y,z,u,v,w,wt
*
* Usually called when reading phsp header info.
* It must be called once for every variable
*
*  result = -1 means the source's header file does not exist
*               or source was not properly initialized (call iaea_new_...)
*  result = -2 means the index is out of range ( 0 <= index < 7 )
*  result = -3 means that the parameter indicated by index is not a constant
*************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_constant_variable(const IAEA_I32 *id, const IAEA_I32 *index,
                                 IAEA_Float *constant, IAEA_I32 *result)
{
      *result=0;
      // No header found
      if(p_iaea_header[*id]->fheader == NULL) {*result = -1; return;}

      if((*index < 0) || (*index > 6) ) {*result = -2; return;}

      if( p_iaea_header[*id]->record_contents[*index] == 0) {
          *constant=p_iaea_header[*id]->record_constant[*index];
      }
      else { *result=-3;}

      return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_constant_variable_(const IAEA_I32 *id,
                                        const IAEA_I32 *index,
                                     IAEA_Float *constant, IAEA_I32 *result)
{iaea_get_constant_variable(id, index, constant, result); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_constant_variable__(const IAEA_I32 *id,
                                        const IAEA_I32 *index,
                                     IAEA_Float *constant, IAEA_I32 *result)
{iaea_get_constant_variable(id, index, constant, result); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_CONSTANT_VARIABLE(const IAEA_I32 *id,
                                        const IAEA_I32 *index,
                                     IAEA_Float *constant, IAEA_I32 *result)
{iaea_get_constant_variable(id, index, constant, result); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_CONSTANT_VARIABLE_(const IAEA_I32 *id,
                                        const IAEA_I32 *index,
                                     IAEA_Float *constant, IAEA_I32 *result)
{iaea_get_constant_variable(id, index, constant, result); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_CONSTANT_VARIABLE__(const IAEA_I32 *id,
                                        const IAEA_I32 *index,
                                     IAEA_Float *constant, IAEA_I32 *result)
{iaea_get_constant_variable(id, index, constant, result); }

/*****************************************************************************
* Get n_indep_particles number of statistically independent particles read
* so far from the Source with Id id.
*
* Set n_indep_particles to negative if such source does not exist.
******************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_used_original_particles(const IAEA_I32 *id, IAEA_I64 *n_indep_particles)
{
      // No header found
      if(p_iaea_header[*id]->fheader == NULL) {*n_indep_particles = -1; return;}

      // (Number of electron histories for linacs)
      *n_indep_particles = p_iaea_header[*id]->read_indep_histories;
      return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_used_original_particles_(const IAEA_I32 *id, IAEA_I64 *n_indep_particles)
{ iaea_get_used_original_particles(id, n_indep_particles); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_used_original_particles__(const IAEA_I32 *id, IAEA_I64 *n_indep_particles)
{ iaea_get_used_original_particles(id, n_indep_particles); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_USED_ORIGINAL_PARTICLES(const IAEA_I32 *id, IAEA_I64 *n_indep_particles)
{ iaea_get_used_original_particles(id, n_indep_particles); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_USED_ORIGINAL_PARTICLES_(const IAEA_I32 *id, IAEA_I64 *n_indep_particles)
{ iaea_get_used_original_particles(id, n_indep_particles); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_USED_ORIGINAL_PARTICLES__(const IAEA_I32 *id, IAEA_I64 *n_indep_particles)
{ iaea_get_used_original_particles(id, n_indep_particles); }

/*****************************************************************************
* Get Total Number of Original Particles from the Source with Id id.
*
* For a typical linac it should be equal to the total number of electrons
* incident on the primary target.
*
* Set number_of_original_particles to negative if such source does not exist.
******************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_total_original_particles(const IAEA_I32 *id,
                                       IAEA_I64 *number_of_original_particles)
{
      // No header found
      if(p_iaea_header[*id]->fheader == NULL)
          {*number_of_original_particles = -1; return;}

      *number_of_original_particles = p_iaea_header[*id]->orig_histories;

      return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_total_original_particles_(const IAEA_I32 *id,
                                               IAEA_I64 *number_of_original_particles)
{ iaea_get_total_original_particles(id, number_of_original_particles); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_total_original_particles__(const IAEA_I32 *id,
                                                IAEA_I64 *number_of_original_particles)
{ iaea_get_total_original_particles(id, number_of_original_particles); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_TOTAL_ORIGINAL_PARTICLES(const IAEA_I32 *id,
                                              IAEA_I64 *number_of_original_particles)
{ iaea_get_total_original_particles(id, number_of_original_particles); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_TOTAL_ORIGINAL_PARTICLES_(const IAEA_I32 *id,
                                               IAEA_I64 *number_of_original_particles)
{ iaea_get_total_original_particles(id, number_of_original_particles); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_TOTAL_ORIGINAL_PARTICLES__(const IAEA_I32 *id,
                                                IAEA_I64 *number_of_original_particles)
{ iaea_get_total_original_particles(id, number_of_original_particles); }

/*****************************************************************************
* Set Total Number of Original Particles for the Source with Id id.
*
* For a typical linac it should be equal to the total number of electrons
* incident on the primary target.
*
* Set number_of_original_particles to negative if such source does not exist.
******************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_total_original_particles(const IAEA_I32 *id,
                                       IAEA_I64 *number_of_original_particles)
{
      // No header found
      if(p_iaea_header[*id]->fheader == NULL)
            {*number_of_original_particles = -1; return;}

      // Bug corrected, RCN, dec. 2006
      p_iaea_header[*id]->orig_histories = *number_of_original_particles;
      return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_total_original_particles_(const IAEA_I32 *id,
                                               IAEA_I64 *number_of_original_particles)
{ iaea_set_total_original_particles(id, number_of_original_particles); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_total_original_particles__(const IAEA_I32 *id,
                                                IAEA_I64 *number_of_original_particles)
{ iaea_set_total_original_particles(id, number_of_original_particles); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_TOTAL_ORIGINAL_PARTICLES(const IAEA_I32 *id,
                                              IAEA_I64 *number_of_original_particles)
{ iaea_set_total_original_particles(id, number_of_original_particles); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_TOTAL_ORIGINAL_PARTICLES_(const IAEA_I32 *id,
                                               IAEA_I64 *number_of_original_particles)
{ iaea_set_total_original_particles(id, number_of_original_particles); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_TOTAL_ORIGINAL_PARTICLES__(const IAEA_I32 *id,
                                                IAEA_I64 *number_of_original_particles)
{ iaea_set_total_original_particles(id, number_of_original_particles); }

/**************************************************************************
* Partitioning for parallel runs
*
* i_parallel is the job number, i_chunk the calculation chunk,
* n_chunk the total number of calculation chunks. This function
* should divide the available phase space of source with Id id
* into n_chunk equal portions and from now on deliver particles
* from the i_chunk-th portion. (i_chunk must be between 1 and n_chunk)
* The extra parameter i_parallel is needed
* for the cases where the source is an event generator and should
* be used to adjust the random number sequence.
* The variable result should be set to 0 if everything went smoothly,
* or to some error code if it didnt.
**************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
//(MACG) void iaea_set_parallel(const IAEA_I32 *id, const IAEA_I32 *i_parallel,
void iaea_set_parallel(const IAEA_I32 *id, const IAEA_I32* /*i_parallel*/,
                       const IAEA_I32 *i_chunk, const IAEA_I32 *n_chunk,
                                       IAEA_I32 *result)
{
   if(p_iaea_header[*id]->fheader == NULL) {*result = -1; return;}
   if(*n_chunk <= 0) {*result = -2; return;}
   if( (*i_chunk < 1) || (*i_chunk > *n_chunk) ) {*result = -3; return;}

   if(p_iaea_header[*id]->file_type == 1)
   {
         // set i_parallel for event generators
         *result = 0;
         return;
   }

   IAEA_I64 nrecords =  p_iaea_header[*id]->nParticles;
   IAEA_I32 record_length =  p_iaea_header[*id]->record_length;
   // IAEA_I32 number_record_per_chunk = (IAEA_I32)nrecords/(*n_chunk); // changed, May 2011
   IAEA_I64 number_record_per_chunk = nrecords/(*n_chunk);


   // IAEA_I32 offset = ((*i_chunk)-1)*record_length * number_record_per_chunk;	// changed, May 2011
   IAEA_I64 offset = ((*i_chunk)-1)*record_length * number_record_per_chunk;
   /*
   SEEK_CUR   Current position of file pointer
   SEEK_END   End of file
   SEEK_SET   Beginning of file
   stream   Pointer to FILE structure
   offset   Number of bytes from origin
   origin   Initial position
   */
   if( fseek(p_iaea_record[*id]->p_file, offset ,SEEK_SET) == 0)
   {
         // IAEA_I32 pos = ftell(p_iaea_record[*id]->p_file);  // changed, May 2011
	     IAEA_I64 pos = (IAEA_I64)ftell(p_iaea_record[*id]->p_file);
	 (void) pos; //(MACG) -Wunused-variable
         *result = 0;
         return;
   }

   *result = -1;
   return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_parallel_(const IAEA_I32 *id,
                                      const IAEA_I32 *i_parallel,
                                      const IAEA_I32 *i_chunk,
                                      const IAEA_I32 *n_chunk,
                                            IAEA_I32 *is_ok)
{ iaea_set_parallel(id, i_parallel, i_chunk, n_chunk, is_ok); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_parallel__(const IAEA_I32 *id,
                                      const IAEA_I32 *i_parallel,
                                      const IAEA_I32 *i_chunk,
                                      const IAEA_I32 *n_chunk,
                                            IAEA_I32 *is_ok)
{ iaea_set_parallel(id, i_parallel, i_chunk, n_chunk, is_ok); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_PARALLEL(const IAEA_I32 *id,
                                      const IAEA_I32 *i_parallel,
                                      const IAEA_I32 *i_chunk,
                                      const IAEA_I32 *n_chunk,
                                            IAEA_I32 *is_ok)
{ iaea_set_parallel(id, i_parallel, i_chunk, n_chunk, is_ok); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_PARALLEL_(const IAEA_I32 *id,
                                      const IAEA_I32 *i_parallel,
                                      const IAEA_I32 *i_chunk,
                                      const IAEA_I32 *n_chunk,
                                            IAEA_I32 *is_ok)
{ iaea_set_parallel(id, i_parallel, i_chunk, n_chunk, is_ok); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_PARALLEL__(const IAEA_I32 *id,
                                      const IAEA_I32 *i_parallel,
                                      const IAEA_I32 *i_chunk,
                                      const IAEA_I32 *n_chunk,
                                            IAEA_I32 *is_ok)
{ iaea_set_parallel(id, i_parallel, i_chunk, n_chunk, is_ok); }

/**************************************************************************
* check that the file size equals the value of checksum in the header
* and that the byte order of the machine being run on matches that of
* the file
*
* id is the phase space file identifier.  Returns 0 if the size of the file=
* checksum and the byte order = the byte order of the machine being run on.
* Returns -1 if the header does not exist; -2 if the function fseek fails
* for some reason; -3 if there is a file size mismatch; -4 if there is a
* byte order mismatch; -5 if there is a mismatch in both
**************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_check_file_size_byte_order(const IAEA_I32 *id,
                                       IAEA_I32 *result)
{
   if(p_iaea_header[*id]->fheader == NULL) {*result = -1; return;}

   int machine_byte_order = check_byte_order();

   #if (defined WIN32) || (defined WIN64)
     // IAEA_I64 size = _filelengthi64(fileno(p_iaea_record[*id]->p_file));
     struct _stati64 fileStatus;
     _fstati64(fileno(p_iaea_record[*id]->p_file),&fileStatus);
   #else
     struct stat fileStatus;
     fstat(fileno(p_iaea_record[*id]->p_file),&fileStatus);
   #endif

   IAEA_I64 size = fileStatus.st_size;
   printf(" phsp size = %llu\n",size);

   // bug found, changed to filelength use. May 2011
   // avoiding ftell and fseek use for I64 compatibility
   if ( size == p_iaea_header[*id]->checksum )
   {
       if (machine_byte_order==p_iaea_header[*id]->byte_order)
       {
          *result=0;
       }
       else {
          *result=-4;
       }
   }
   else {
       if (machine_byte_order==p_iaea_header[*id]->byte_order)
       {
          *result=-3;
       }
       else {
          *result=-5;
       }
   }

   return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_check_file_size_byte_order_(const IAEA_I32 *id,
                                       IAEA_I32 *result)
{ iaea_check_file_size_byte_order(id, result);}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_check_file_size_byte_order__(const IAEA_I32 *id,
                                       IAEA_I32 *result)
{ iaea_check_file_size_byte_order( id, result);}
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_CHECK_FILE_SIZE_BYTE_ORDER(const IAEA_I32 *id,
                                       IAEA_I32 *result)
{ iaea_check_file_size_byte_order( id, result);}
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_CHECK_FILE_SIZE_BYTE_ORDER_(const IAEA_I32 *id,
                                       IAEA_I32 *result)
{ iaea_check_file_size_byte_order( id, result);}
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_CHECK_FILE_SIZE_BYTE_ORDER__(const IAEA_I32 *id,
                                       IAEA_I32 *result)
{ iaea_check_file_size_byte_order( id, result);}


/**************************************************************************
* setting the pointer to a user-specified record no. in the file
*
* record_num is the user-specified record number passed to the function.
* id is the phase space file identifier.  If record_num = the number of
* particles in the file + 1, then the position is set to the end of the
* file.
* The variable result should be set to 0 if everything went smoothly,
* or to some error code if it didnt.
**************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_record(const IAEA_I32 *id, const IAEA_I64 *record_num,
                                       IAEA_I32 *result)
{
   if(p_iaea_header[*id]->fheader == NULL) {*result = -1; return;}
   if(*record_num <= 0) {*result = -2; return;}
   if(*record_num > p_iaea_header[*id]->nParticles+1) {*result = -3; return;}

   IAEA_I32 record_length =  p_iaea_header[*id]->record_length;

   IAEA_I64 offset = (*record_num-1) * record_length;
   /*
   SEEK_CUR   Current position of file pointer
   SEEK_END   End of file
   SEEK_SET   Beginning of file
   stream   Pointer to FILE structure
   offset   Number of bytes from origin
   origin   Initial position
   */

   if( fseek(p_iaea_record[*id]->p_file, offset ,SEEK_SET) == 0)
   {
     IAEA_I32 pos = ftell(p_iaea_record[*id]->p_file);
     (void) pos; //(MACG)-Wunused-variable
     *result = 0;
     return;
   }

   *result = -1;
   return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_record_(const IAEA_I32 *id,
                                      const IAEA_I64 *record_num,
                                            IAEA_I32 *is_ok)
{ iaea_set_record(id, record_num, is_ok); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_record__(const IAEA_I32 *id,
                                      const IAEA_I64 *record_num,
                                            IAEA_I32 *is_ok)
{ iaea_set_record(id, record_num, is_ok); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_RECORD(const IAEA_I32 *id,
                                      const IAEA_I64 *record_num,
                                            IAEA_I32 *is_ok)
{ iaea_set_record(id, record_num, is_ok); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_RECORD_(const IAEA_I32 *id,
                                      const IAEA_I64 *record_num,
                                            IAEA_I32 *is_ok)
{ iaea_set_record(id, record_num, is_ok); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_SET_RECORD__(const IAEA_I32 *id,
                                      const IAEA_I64 *record_num,
                                            IAEA_I32 *is_ok)
{ iaea_set_record(id, record_num, is_ok); }

/**************************************************************************
* Get a particle
*
* Return the next particle from the sequence of particles from source
* with Id id. Set n_stat to the number of statistically independent
* events since the last call to this function (i.e. n_stat = 0, if
* the particle resulted from the same incident electron, n_stat = 377
* if there were 377 statistically independent events sinc the last particle
* returned, etc.). If this information is not available,
* simply set n_stat to 1 if the particle belongs to a new statistically
* independent event. Set n_stat to -1, if a source with Id id does not
* exist. Set n_stat to -2, if end of file of the phase space source reached

**************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_particle(const IAEA_I32 *id, IAEA_I32 *n_stat,
IAEA_I32 *type, /* particle type */
IAEA_Float *E,  /* kinetic energy in MeV */
IAEA_Float *wt, /* statistical weight */
IAEA_Float *x,
IAEA_Float *y,
IAEA_Float *z,  /* position in cartesian coordinates*/
IAEA_Float *u,
IAEA_Float *v,
IAEA_Float *w,  /* direction in cartesian coordinates*/
IAEA_Float *extra_floats,
IAEA_I32 *extra_ints)
{
      if(feof(p_iaea_record[*id]->p_file)) {
         *n_stat = -2;
         rewind (p_iaea_record[*id]->p_file);
         return;
      }

      if( p_iaea_record[*id]->read_particle() == FAIL ) { *n_stat = -1; return;}

      iaea_record_type *p = p_iaea_record[*id];

      // Corrected on Dec. 2006. Before n_stat was not assigned if
      // (p->iextralong > 0)  and (p_iaea_header[*id]->extralong_contents[j] != 1)

      if( p->IsNewHistory > 0) *n_stat=1;
      else                     *n_stat=0;

      if(p->iextralong > 0) {
          for(int j=0;j<p->iextralong ;j++) {
              // Looking for incremental number of histories
              // (Type 1 of the extralong stored variable)
              if(p_iaea_header[*id]->extralong_contents[j] == 1) {
                  *n_stat = p->extralong[j];
                  p->IsNewHistory = *n_stat;
              }
          }
      }

      *type  = p->particle; /* particle type */
      *E     = p->energy;   /* kinetic energy in MeV */

      if(p->ix > 0) *x = p->x;
      else          *x = p_iaea_header[*id]->record_constant[0];

      if(p->iy > 0) *y = p->y;
      else          *y = p_iaea_header[*id]->record_constant[1];

      if(p->iz > 0) *z = p->z; /* position in cartesian coordinates*/
      else          *z = p_iaea_header[*id]->record_constant[2];

      if(p->iu > 0) *u = p->u;
      else          *u = p_iaea_header[*id]->record_constant[3];

      if(p->iv > 0) *v = p->v;
      else          *v = p_iaea_header[*id]->record_constant[4];

      if(p->iw > 0) *w = p->w; /* direction in cartesian coordinates*/
      else          *w = p_iaea_header[*id]->record_constant[5];

      if(p->iweight > 0) *wt = p->weight;   /* statistical weight */
      else               *wt = p_iaea_header[*id]->record_constant[6];

      for(int k=0;k<p->iextrafloat;k++) extra_floats[k] = p->extrafloat[k];
      for(int j=0;j<p->iextralong ;j++) extra_ints[j] = p->extralong[j];

      /*
        Updating counters including:
        Min and Max Weight per particle,
        Total Weight per particle,
        Min and Max Kinetic energy per particle,
        Min and Max X,Y,Z coordinates,
        Total number of particles,
        Total number of each particle type
        Number of statistically independent histories
      */
      p_iaea_header[*id]->update_counters(p_iaea_record[*id]);

      return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_particle_(const IAEA_I32 *id, IAEA_I32 *n_stat,
IAEA_I32 *type, /* particle type */
IAEA_Float *E,  /* kinetic energy in MeV */
IAEA_Float *wt, /* statistical weight */
IAEA_Float *x,
IAEA_Float *y,
IAEA_Float *z,  /* position in cartesian coordinates*/
IAEA_Float *u,
IAEA_Float *v,
IAEA_Float *w,  /* direction in cartesian coordinates*/
IAEA_Float *extra_floats,
IAEA_I32 *extra_ints)
{ iaea_get_particle(id, n_stat, type,
                                E, wt, x, y, z, u, v, w, extra_floats, extra_ints); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_particle__(const IAEA_I32 *id, IAEA_I32 *n_stat,
IAEA_I32 *type, /* particle type */
IAEA_Float *E,  /* kinetic energy in MeV */
IAEA_Float *wt, /* statistical weight */
IAEA_Float *x,
IAEA_Float *y,
IAEA_Float *z,  /* position in cartesian coordinates*/
IAEA_Float *u,
IAEA_Float *v,
IAEA_Float *w,  /* direction in cartesian coordinates*/
IAEA_Float *extra_floats,
IAEA_I32 *extra_ints)
{ iaea_get_particle(id, n_stat, type,
                                E, wt, x, y, z, u, v, w, extra_floats, extra_ints); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_PARTICLE(const IAEA_I32 *id, IAEA_I32 *n_stat,
IAEA_I32 *type, /* particle type */
IAEA_Float *E,  /* kinetic energy in MeV */
IAEA_Float *wt, /* statistical weight */
IAEA_Float *x,
IAEA_Float *y,
IAEA_Float *z,  /* position in cartesian coordinates*/
IAEA_Float *u,
IAEA_Float *v,
IAEA_Float *w,  /* direction in cartesian coordinates*/
IAEA_Float *extra_floats,
IAEA_I32 *extra_ints)
{ iaea_get_particle(id, n_stat, type,
                                E, wt, x, y, z, u, v, w, extra_floats, extra_ints); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_PARTICLE_(const IAEA_I32 *id, IAEA_I32 *n_stat,
IAEA_I32 *type, /* particle type */
IAEA_Float *E,  /* kinetic energy in MeV */
IAEA_Float *wt, /* statistical weight */
IAEA_Float *x,
IAEA_Float *y,
IAEA_Float *z,  /* position in cartesian coordinates*/
IAEA_Float *u,
IAEA_Float *v,
IAEA_Float *w,  /* direction in cartesian coordinates*/
IAEA_Float *extra_floats,
IAEA_I32 *extra_ints)
{ iaea_get_particle(id, n_stat, type,
                                E, wt, x, y, z, u, v, w, extra_floats, extra_ints); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_GET_PARTICLE__(const IAEA_I32 *id, IAEA_I32 *n_stat,
IAEA_I32 *type, /* particle type */
IAEA_Float *E,  /* kinetic energy in MeV */
IAEA_Float *wt, /* statistical weight */
IAEA_Float *x,
IAEA_Float *y,
IAEA_Float *z,  /* position in cartesian coordinates*/
IAEA_Float *u,
IAEA_Float *v,
IAEA_Float *w,  /* direction in cartesian coordinates*/
IAEA_Float *extra_floats,
IAEA_I32 *extra_ints)
{ iaea_get_particle(id, n_stat, type,
                                E, wt, x, y, z, u, v, w, extra_floats, extra_ints); }

/**************************************************************************
* Write a particle
* n_stat = 0 for a secondary particle
* n_stat > 0 for an independent particle
*
* Write a particle to the source with Id id.
* Set n_stat to -1, if ERROR (source with Id id does not exist).
**************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_write_particle(const IAEA_I32 *id, IAEA_I32 *n_stat,
const IAEA_I32 *type, /* particle type */
const IAEA_Float *E,  /* kinetic energy in MeV */
const IAEA_Float *wt, /* statistical weight */
const IAEA_Float *x,
const IAEA_Float *y,
const IAEA_Float *z,  /* position in cartesian coordinates*/
const IAEA_Float *u,
const IAEA_Float *v,
const IAEA_Float *w,  /* direction in cartesian coordinates*/
const IAEA_Float *extra_floats,
const IAEA_I32 *extra_ints)
{

      iaea_record_type *p = p_iaea_record[*id];

      if( *n_stat > 0 ) p->IsNewHistory = *n_stat;
      else              p->IsNewHistory = 0;

      p->particle = (short)*type; /* particle type */
      p->energy   = *E;    /* kinetic energy in MeV */
      if(p->iweight > 0) p->weight = *wt;   /* statistical weight */
      if(p->ix > 0) p->x = *x; /* position in cartesian coordinates*/
      if(p->iy > 0) p->y = *y;
      if(p->iz > 0) p->z = *z;
      if(p->iu > 0) p->u = *u; /* direction in cartesian coordinates*/
      if(p->iv > 0) p->v = *v;
      if(p->iw > 0) p->w = *w;

      for(int k=0;k<p->iextrafloat;k++) p->extrafloat[k] = extra_floats[k];
      for(int j=0;j<p->iextralong ;j++)  p->extralong[j] = extra_ints[j];

      if( p->write_particle() == FAIL ) {*n_stat = -1; return;}

      /*
        Updating counters including:
        Min and Max Weight per particle,
        Total Weight per particle,
        Min and Max Kinetic energy per particle,
        Min and Max X,Y,Z coordinates,
        Total number of particles,
        Total number of each particle type
        Number of statistically independent histories so far
      */
      p_iaea_header[*id]->update_counters(p_iaea_record[*id]);

      return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_write_particle_(const IAEA_I32 *id, IAEA_I32 *n_stat,
const IAEA_I32 *type, /* particle type */
const IAEA_Float *E,  /* kinetic energy in MeV */
const IAEA_Float *wt, /* statistical weight */
const IAEA_Float *x,
const IAEA_Float *y,
const IAEA_Float *z,  /* position in cartesian coordinates*/
const IAEA_Float *u,
const IAEA_Float *v,
const IAEA_Float *w,  /* direction in cartesian coordinates*/
const IAEA_Float *extra_floats,
const IAEA_I32 *extra_ints)
{ iaea_write_particle(id, n_stat, type,
                                E, wt, x, y, z, u, v, w, extra_floats, extra_ints); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_write_particle__(const IAEA_I32 *id, IAEA_I32 *n_stat,
const IAEA_I32 *type, /* particle type */
const IAEA_Float *E,  /* kinetic energy in MeV */
const IAEA_Float *wt, /* statistical weight */
const IAEA_Float *x,
const IAEA_Float *y,
const IAEA_Float *z,  /* position in cartesian coordinates*/
const IAEA_Float *u,
const IAEA_Float *v,
const IAEA_Float *w,  /* direction in cartesian coordinates*/
const IAEA_Float *extra_floats,
const IAEA_I32 *extra_ints)
{ iaea_write_particle(id, n_stat, type,
                                E, wt, x, y, z, u, v, w, extra_floats, extra_ints); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_WRITE_PARTICLE(const IAEA_I32 *id, IAEA_I32 *n_stat,
const IAEA_I32 *type, /* particle type */
const IAEA_Float *E,  /* kinetic energy in MeV */
const IAEA_Float *wt, /* statistical weight */
const IAEA_Float *x,
const IAEA_Float *y,
const IAEA_Float *z,  /* position in cartesian coordinates*/
const IAEA_Float *u,
const IAEA_Float *v,
const IAEA_Float *w,  /* direction in cartesian coordinates*/
const IAEA_Float *extra_floats,
const IAEA_I32 *extra_ints)
{ iaea_write_particle(id, n_stat, type,
                                E, wt, x, y, z, u, v, w, extra_floats, extra_ints); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_WRITE_PARTICLE_(const IAEA_I32 *id, IAEA_I32 *n_stat,
const IAEA_I32 *type, /* particle type */
const IAEA_Float *E,  /* kinetic energy in MeV */
const IAEA_Float *wt, /* statistical weight */
const IAEA_Float *x,
const IAEA_Float *y,
const IAEA_Float *z,  /* position in cartesian coordinates*/
const IAEA_Float *u,
const IAEA_Float *v,
const IAEA_Float *w,  /* direction in cartesian coordinates*/
const IAEA_Float *extra_floats,
const IAEA_I32 *extra_ints)
{ iaea_write_particle(id, n_stat, type,
                                E, wt, x, y, z, u, v, w, extra_floats, extra_ints); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_WRITE_PARTICLE__(const IAEA_I32 *id, IAEA_I32 *n_stat,
const IAEA_I32 *type, /* particle type */
const IAEA_Float *E,  /* kinetic energy in MeV */
const IAEA_Float *wt, /* statistical weight */
const IAEA_Float *x,
const IAEA_Float *y,
const IAEA_Float *z,  /* position in cartesian coordinates*/
const IAEA_Float *u,
const IAEA_Float *v,
const IAEA_Float *w,  /* direction in cartesian coordinates*/
const IAEA_Float *extra_floats,
const IAEA_I32 *extra_ints)
{ iaea_write_particle(id, n_stat, type,
                                E, wt, x, y, z, u, v, w, extra_floats, extra_ints); }

/***************************************************************************
* Destroy a source
*
* This function de-initializes the source with Id id, closing all open
* files, deallocating memory, etc. Nothing happens if a source with that
* id does not exist. Header is updated.
****************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_destroy_source(const IAEA_I32 *source_ID, IAEA_I32 *result)
{
   if(*source_ID > MAX_NUM_SOURCES) { *result = -98 ; return;} // Too big phsp ID
   if(*source_ID < 0)               { *result = -97 ; return;} // wrong ID number

   if(p_iaea_header[*source_ID]->fheader == NULL) {*result = -1; return;}

  /* Write an IAEA header */
   // For read-only files nothing happens
   p_iaea_header[*source_ID]->write_header();

   // Closing header file
   fclose(p_iaea_header[*source_ID]->fheader);
   // Deallocating IAEA phsp header
   free(p_iaea_header[*source_ID]);

   // Closing phsp file
   fclose(p_iaea_record[*source_ID]->p_file);
   // Deallocating IAEA record
   free(p_iaea_record[*source_ID]);

   __iaea_source_used[*source_ID] = false;

   *result = 1; // Return OK

   return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_destroy_source_(const IAEA_I32 *source_ID, IAEA_I32 *result)
{ iaea_destroy_source(source_ID, result); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_destroy_source__(const IAEA_I32 *source_ID, IAEA_I32 *result)
{ iaea_destroy_source(source_ID, result); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_DESTROY_SOURCE(const IAEA_I32 *source_ID, IAEA_I32 *result)
{ iaea_destroy_source(source_ID, result); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_DESTROY_SOURCE_(const IAEA_I32 *source_ID, IAEA_I32 *result)
{ iaea_destroy_source(source_ID, result); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_DESTROY_SOURCE__(const IAEA_I32 *source_ID, IAEA_I32 *result)
{ iaea_destroy_source(source_ID, result); }

/***************************************************************************
* Print the current header associated to source id
*
* result is set to negative if phsp source does not exist.
****************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_print_header(const IAEA_I32 *source_ID, IAEA_I32 *result)
{
   if(p_iaea_header[*source_ID]->fheader == NULL) {*result = -1; return;}

  /* Print current IAEA header for source id */
   p_iaea_header[*source_ID]->print_header();

   *result = 1; // Return OK

   return;
}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_print_header_(const IAEA_I32 *source_ID, IAEA_I32 *result)
{ iaea_print_header(source_ID, result); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_print_header__(const IAEA_I32 *source_ID, IAEA_I32 *result)
{ iaea_print_header(source_ID, result); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_PRINT_HEADER(const IAEA_I32 *source_ID, IAEA_I32 *result)
{ iaea_print_header(source_ID, result); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_PRINT_HEADER_(const IAEA_I32 *source_ID, IAEA_I32 *result)
{ iaea_print_header(source_ID, result); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_PRINT_HEADER__(const IAEA_I32 *source_ID, IAEA_I32 *result)
{ iaea_print_header(source_ID, result); }
//IAEA_EXTERN_C IAEA_EXPORT //(MACG) to silent -Wduplicate-decl-specifier (bug?)

/***************************************************************************
* Copy header of the source_id to the header of the destiny_id
*
* result is set to negative if phsp source or destiny do not exist.
****************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_copy_header(const IAEA_I32 *source_ID,
                                const IAEA_I32 *destiny_ID, IAEA_I32 *result)
{
   if(p_iaea_header[*source_ID]->fheader == NULL) {*result = -1; return;}
   if(p_iaea_header[*destiny_ID]->fheader == NULL) {*result = -1; return;}

   // Selective copy of string variables

      p_iaea_header[*destiny_ID]->checksum =
            p_iaea_header[*source_ID]->checksum ;
      p_iaea_header[*destiny_ID]->record_length =
            p_iaea_header[*source_ID]->record_length ;
      p_iaea_header[*destiny_ID]->byte_order =
            p_iaea_header[*source_ID]->byte_order ;

// ******************************************************************************
// 2. Mandatory description of the phsp

      strcpy(p_iaea_header[*destiny_ID]->coordinate_system_description,
            p_iaea_header[*source_ID]->coordinate_system_description) ;

      int file_type = p_iaea_header[*source_ID]->file_type;
      if(file_type == 1)
      {
            // For event generators
            strcpy(p_iaea_header[*destiny_ID]->input_file_for_event_generator,
                  p_iaea_header[*source_ID]->input_file_for_event_generator) ;
            *result = 1; // Return OK
            return;
      }

      p_iaea_header[*destiny_ID]->orig_histories =
            p_iaea_header[*source_ID]->orig_histories ;

// ******************************************************************************
// 3. Mandatory additional information
      /*********************************************/
      strcpy(p_iaea_header[*destiny_ID]->machine_type,
            p_iaea_header[*source_ID]->machine_type);

      strcpy(p_iaea_header[*destiny_ID]->MC_code_and_version,
            p_iaea_header[*source_ID]->MC_code_and_version);

      p_iaea_header[*destiny_ID]->global_photon_energy_cutoff =
            p_iaea_header[*source_ID]->global_photon_energy_cutoff;

      p_iaea_header[*destiny_ID]->global_particle_energy_cutoff =
            p_iaea_header[*source_ID]->global_particle_energy_cutoff;

      strcpy(p_iaea_header[*destiny_ID]->transport_parameters,
            p_iaea_header[*source_ID]->transport_parameters);

// ******************************************************************************
// 4. Optional description
      strcpy(p_iaea_header[*destiny_ID]->beam_name,
            p_iaea_header[*source_ID]->beam_name);
      strcpy(p_iaea_header[*destiny_ID]->field_size,
            p_iaea_header[*source_ID]->field_size);
      strcpy(p_iaea_header[*destiny_ID]->nominal_SSD,
            p_iaea_header[*source_ID]->nominal_SSD);
      strcpy(p_iaea_header[*destiny_ID]->variance_reduction_techniques,
            p_iaea_header[*source_ID]->variance_reduction_techniques);
      strcpy(p_iaea_header[*destiny_ID]->initial_source_description,
            p_iaea_header[*source_ID]->initial_source_description);

      // Documentation sub-section
      /*********************************************/
      strcpy(p_iaea_header[*destiny_ID]->MC_input_filename,
            p_iaea_header[*source_ID]->MC_input_filename);
      strcpy(p_iaea_header[*destiny_ID]->published_reference,
            p_iaea_header[*source_ID]->published_reference);
      strcpy(p_iaea_header[*destiny_ID]->authors,
            p_iaea_header[*source_ID]->authors);
      strcpy(p_iaea_header[*destiny_ID]->institution,
            p_iaea_header[*source_ID]->institution);
      strcpy(p_iaea_header[*destiny_ID]->link_validation,
            p_iaea_header[*source_ID]->link_validation);
      strcpy(p_iaea_header[*destiny_ID]->additional_notes,
            p_iaea_header[*source_ID]->additional_notes);

    *result = 1; // Return OK

    return;

}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_copy_header_(const IAEA_I32 *source_ID,
                                     const IAEA_I32 *destiny_ID,
                                           IAEA_I32 *result)
{ iaea_copy_header(source_ID, destiny_ID, result); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_copy_header__(const IAEA_I32 *source_ID,
                                     const IAEA_I32 *destiny_ID,
                                           IAEA_I32 *result)
{ iaea_copy_header(source_ID, destiny_ID, result); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_COPY_HEADER(const IAEA_I32 *source_ID,
                                     const IAEA_I32 *destiny_ID,
                                           IAEA_I32 *result)
{ iaea_copy_header(source_ID, destiny_ID, result); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_COPY_HEADER_(const IAEA_I32 *source_ID,
                                     const IAEA_I32 *destiny_ID,
                                           IAEA_I32 *result)
{ iaea_copy_header(source_ID, destiny_ID, result); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_COPY_HEADER__(const IAEA_I32 *source_ID,
                                     const IAEA_I32 *destiny_ID,
                                           IAEA_I32 *result)
{ iaea_copy_header(source_ID, destiny_ID, result); }


/***************************************************************************
* Update header of the source_id
*
* result is set to negative if phsp source does not exist.
****************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_update_header(const IAEA_I32 *source_ID, IAEA_I32 *result)
{
   if(p_iaea_header[*source_ID]->fheader == NULL) {*result = -1; return;}

  /* Write an IAEA header */
   // For read-only files nothing happens
   p_iaea_header[*source_ID]->write_header();

   *result = 1; // Return OK
   return;

}
IAEA_EXTERN_C IAEA_EXPORT
void iaea_update_header_(const IAEA_I32 *source_ID, IAEA_I32 *result)
{ iaea_update_header(source_ID, result); }
IAEA_EXTERN_C IAEA_EXPORT
void iaea_update_header__(const IAEA_I32 *source_ID, IAEA_I32 *result)
{ iaea_update_header(source_ID, result); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_UPDATE_HEADER(const IAEA_I32 *source_ID, IAEA_I32 *result)
{ iaea_update_header(source_ID, result); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_UPDATE_HEADER_(const IAEA_I32 *source_ID, IAEA_I32 *result)
{ iaea_update_header(source_ID, result); }
IAEA_EXTERN_C IAEA_EXPORT
void IAEA_UPDATE_HEADER__(const IAEA_I32 *source_ID, IAEA_I32 *result)
{ iaea_update_header(source_ID, result); }
