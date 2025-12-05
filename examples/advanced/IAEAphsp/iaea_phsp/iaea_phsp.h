
/* 
 * INTERFACE FOR IAEA PHSP ROUTINES (CONTAINS ONLY DECLARATIONS)
 *
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

#ifndef IAEA_PHSP
#define IAEA_PHSP

#include "iaea_config.h"

/************************************************************************
* Initialization 
*
* Given a file name header_file of length hf_length, initialize a
* new IAEA particle source, assign an unique source_ID to it and return 
* this Id in result. Dont assume header_file is null-terminated as the
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
IAEA_EXTERN_C IAEA_EXPORT 
void iaea_new_source(IAEA_I32 *source_ID, char *header_file,   
                     const IAEA_I32 *access, IAEA_I32 *result, 
                     int hf_length);

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
                             IAEA_I64 *n_particle);

/************************************************************************
* Maximum energy 
*
* Return the maximum energy of an initialized IAEA source with Id id
* in Emax. Set Emax to negative if a source with that Id does not exist.
************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT 
void iaea_get_maximum_energy(const IAEA_I32 *id, IAEA_Float *Emax);

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
                            IAEA_I32 *n_extra_int);

/*************************************************************************
* Number of additional floats and integers to be stored
*
* Set the number of additional floats in n_extra_float and the number
* of additional integers in n_extra_integer for the source with Id id
* to be stored in the corresponding file.
*************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT 
void iaea_set_extra_numbers(const IAEA_I32 *id, IAEA_I32 *n_extra_float,
                                  IAEA_I32 *n_extra_int);

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
                                            IAEA_I32 *type);

/********************************************************************************
* Set a type type of the extra float variable corresponding to the "index" number
* for a corresponding header of the phsp "id". Index is running from zero.            
*
* The current list of types for extra float variables is:
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
                                             IAEA_I32 *type);

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
      IAEA_I32 extralong_types[], IAEA_I32 extrafloat_types[]);
                                                                              
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
                                IAEA_Float *constant);

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
                                 IAEA_Float *constant, IAEA_I32 *result);

/*****************************************************************************
* Get n_indep_particles number of statistically independent particles read 
* so far from the Source with Id id.                                  
*
* Set n_indep_particles to negative if such source does not exist.
******************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_get_used_original_particles(const IAEA_I32 *id, 
                                            IAEA_I64 *n_indep_particles);

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
                                             IAEA_I64 *number_of_original_particles);

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
                                       IAEA_I64 *number_of_original_particles);

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
* The variable is_ok should be set to 0 if everything went smoothly,
* or to some error code if it didnt.
**************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT 
void iaea_set_parallel(const IAEA_I32 *id, const IAEA_I32 *i_parallel,
                       const IAEA_I32 *i_chunk, const IAEA_I32 *n_chunk, 
                       IAEA_I32 *is_ok);

/**************************************************************************
* setting the pointer to a user-specified record no. in the file
*
* record_num is the user-specified record number passed to the function.
* id is the phase space file identifier.
* The variable result should be set to 0 if everything went smoothly,
* or to some error code if it didnt.
**************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_set_record(const IAEA_I32 *id, const IAEA_I64 *record_num,
                           IAEA_I32 *result);

/**************************************************************************
* check that the file size equals the value of checksum in the header
*
* id is the phase space file identifier.  If the size of the phase space
* file is not equal to checksum, then result returns -1, otherwise result
* is set to 0.
**************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT
void iaea_check_file_size_byte_order(const IAEA_I32 *id, IAEA_I32 *result);

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
IAEA_I32 *extra_ints);

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
const IAEA_I32 *extra_ints);

/***************************************************************************
* Destroy a source 
*
* This function de-initializes the source with Id id, closing all open
* files, deallocating memory, etc. Nothing happens if a source with that
* id does not exist. Header is updated.
****************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT 
void iaea_destroy_source(const IAEA_I32 *source_ID, IAEA_I32 *result);

/***************************************************************************
* Print the current header associated to source id
*
* result is set to negative if phsp source does not exist.
****************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT 
void iaea_print_header(const IAEA_I32 *source_ID, IAEA_I32 *result);

/***************************************************************************
* Copy header of the source_id to the header of the destiny_id 
****************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT 
void iaea_copy_header(const IAEA_I32 *source_ID, const IAEA_I32 *destiny_ID, 
                      IAEA_I32 *result);

/***************************************************************************
* Update header of the source_id 
****************************************************************************/
IAEA_EXTERN_C IAEA_EXPORT 
void iaea_update_header(const IAEA_I32 *source_ID, IAEA_I32 *result);

#endif
