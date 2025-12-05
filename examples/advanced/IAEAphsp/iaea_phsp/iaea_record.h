#include <cstdio>

#ifndef IAEA_RECORD
#define IAEA_RECORD

#include "utilities.h"
#include "iaea_config.h"

/* *********************************************************************** */
// defines

// To use additional float or integers

#ifndef NUM_EXTRA_FLOAT
  #define NUM_EXTRA_FLOAT 10 // Maximum 10 extra float stored
#endif

#ifndef NUM_EXTRA_LONG
  #define NUM_EXTRA_LONG  10 // Maximum 10 extra long stored
#endif

#define MAX_NUM_PARTICLES 5 // 1 photons 
                            // 2 electrons
                            // 3 positrons
                            // 4 neutrons
                            // 5 protons
#define MAX_NUM_SOURCES 30

#define OK     0
#define FAIL  -1

/* *********************************************************************** */
// structures

struct iaea_record_type
{
  FILE *p_file;   // phase space file pointer   

  short particle; // mandatory       (photon:1 electron:2 positron:3 neutron:4 proton:5 ...)
  
  float  energy;  // mandatory
  
  IAEA_I32 IsNewHistory; // coded as sign of energy 
                         //  Type changed from short to IAEA_I32 to store EGS n_stat

  float x;       int ix;       
  float y;       int iy;
  float z;       int iz;
  float u;       int iu;
  float v;       int iv;
  float w;       int iw;      // sign of w coded as sign of code
  float weight;  int iweight;

  short iextrafloat; 
  short iextralong;  

  float extrafloat[NUM_EXTRA_FLOAT];  // (default: no extra float stored)
  IAEA_I32 extralong[NUM_EXTRA_LONG];      // (default: one extra long stored)

public:
      short read_particle();
      short write_particle();
      short initialize();
};

#endif
