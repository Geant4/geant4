/*
   Copyright 2000-2003 Virginia Commonwealth University


   Advisory:
1. The authors make no claim of accuracy of the information in these files or the 
   results derived from use of these files.
2. You are not allowed to re-distribute these files or the information contained within.
3. This methods and information contained in these files was interpreted by the authors 
   from various sources.  It is the users sole responsibility to verify the accuracy 
   of these files.  
4. If you find a error within these files, we ask that you to contact the authors 
   and the distributor of these files.
5. We ask that you acknowledge the source of these files in publications that use results
   derived from the input 

   Please contact us if you have any questions 
*/
/* Header File for General Utilities for CPP programs
	File Created:
		18-December-1995: Combined ok_check.cpp and some open_file
	Modification History:
		 09-feb-1996: JVS: add pprintf
       07-Nov-1996: JVS: FAIL_SAFE definition changed, NULL=0
       06-Jan-1998: JVS: Add array_read
       07-May-1998: JVS: add myerrno
       11-Sept-1998: JVS: add max, min definitions
      // 28-Oct-1998: PJK: added MAX_NUM_FIELDS for output_path in case_info.h
      03-Dec-1998: JVS: Add #ifndef UTILITIES_H_INCLUDED to ensure single inclusion of the file
                        All comments "c" compliant
       02-Mar-1999: JVS: Add eprintf_mode
       Dec 3, 1999: JVS: Add check_byte_order
       April 24, 2001: JVS: Add cp
       August 24, 2001: JVS: Change OK from 1 to 0...so exit(OK) is unix standard normal...
       Feb 10, 2005: JVS:Add writeLittleEndianBinaryFile()
       April 21, 2005: JVS: Add readBinaryDataFromFile()
*/
#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED

#define OK     0
#define ERROR -1
#define FAIL  -1
#define ON     1
#define OFF    0
#define FAIL_SAFE -999        /* used to be NULL**JVS 11/7/96** */
#define MAX_STR_LEN 512       /* maximum length of a string */
#define MAX_BUFFER_SIZE 16384 /* maximum number of characters in the output buffer */
#define STR_NULL   ((char) 0)
/* #define MAX_NUM_FIELDS 100 */
#ifndef N_DATA
#define N_DATA 4096 /* maximum number of datapoints */
#endif
#ifdef MAIN
int myerrno = OK;
#else
extern int myerrno;
#endif
#ifdef MAIN
int eprintf_mode = ON;
#else
extern int eprintf_mode;
#endif
#ifndef MAIN
extern /* create global pbuffer */
#endif          
char *pbuffer;  /* pointer to buffer for output of run info for failures */

#ifndef	ENDIAN_H_INCLUDED
#define	ENDIAN_H_INCLUDED


/* Definitions for byte order, according to significance of bytes, from low
   addresses to high addresses.  The value is what you get by putting '4'
   in the most significant byte, '3' in the second most significant byte,
   '2' in the second least significant byte, and '1' in the least
   significant byte.  */

#define	__LITTLE_ENDIAN	1234
#define	__BIG_ENDIAN	4321
#define	__PDP_ENDIAN	3412

//(MACG) endian macros defined in macOS 15, silent -Wmacro-redefined warning
#if !defined(LITTLE_ENDIAN) && defined(__LITTLE_ENDIAN)
  #define LITTLE_ENDIAN   __LITTLE_ENDIAN
#endif
#if !defined(BIG_ENDIAN) && defined(__BIG_ENDIAN)
  #define BIG_ENDIAN      __BIG_ENDIAN
#endif
#if !defined(PDP_ENDIAN) && defined(__PDP_ENDIAN)
  #define PDP_ENDIAN      __PDP_ENDIAN
#endif
//(MACG) --- end of changes
#define UNKNOWN_ENDIAN  0000

#endif	/* endian.h */

/* #define max and min */
/* (MACG) using min and max functions through algorithm c++ library
#ifndef MINMAX_DEFINED
#define MINMAX_DEFINED
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#define max(a,b)    (((a) > (b)) ? (a) : (b))
//#define __max       max
//#define __min       min
#endif (MACG)*/ /*MINMAX_DEFINED */

/* *************************************************************************** */
float reverse_float_byte_order(float xold);
short reverse_short_byte_order(short xold);
int reverse_int_byte_order(int xold);
int advance(char *istr, int *sval, int len);
int check_byte_order(void);
int clean_name(char *tmp_path, char *opath);
int clean_name(char *);
int copy(char *SourceFile, char *DestinationFile);
//(MACG) int eprintf(char *fmt, ... );
int eprintf(const char *fmt, ... );
int ok_check(void);
int ok_checks(char *string);
//(MACG) FILE *open_file(char *filename,char *extension, char *access);
FILE *open_file(char *filename,const char *extension,const char *access);
int pprintf(char *fmt, ... );
int latex_string(char *string, char *nstring);
void print_runtime_info(int argc, char *argv[]);
void allocate_pbuffer(void);
float interpolate(float xh, float xl, float xm, float yh, float yl);
int array_read(FILE *istrm, float *array, int max_array);
int array_read(char *in_string, float *array, int max_array);
int view_errors(void);
int writeBigEndianBinaryFile(char *doseFileName,  int nDoseArray, float *doseArray);
int writeLittleEndianBinaryFile(char *doseFileName,  int nDoseArray, float *doseArray);
int writeBinaryFile(char *doseFileName, int nDoseArray, float *doseArray, int swab_flag);
int writeBinaryDataToFile(FILE *outputStream, int nArray, float *array, int swab_flag);
int readBinaryDataFromFile(FILE *iStream, int nItemsToRead, float **arrayToRead, int swab_flag);
int readBinaryDataFromFile(FILE *iStream, int nItemsToRead, float *inputArray, int swab_flag);
// RCN added 
int fget_c_string(char *string, int Max_Str_Len, FILE *fspec);
int get_string(FILE *fspec, char *string);
#endif
