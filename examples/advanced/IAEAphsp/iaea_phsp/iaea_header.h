#ifndef IAEA_HEADER
#define IAEA_HEADER

/* *********************************************************************** */
#include "iaea_record.h"

// defines
#define SEGMENT_BEG_TOKEN '$'
#define SEGMENT_END_TOKEN ':'
#ifndef MAX_STR_LEN
  #define MAX_STR_LEN 512   /* maximum length of a string */
#endif
#define MAX_NUMB_LINES 30   /* maximum number of lines in a block */

#define MAX_NUMB_EXTRALONG_TYPES 7 /* maximum number of extra long allowed */
//  0: User defined generic type
//  1: Incremental history number n_hist
//     n_hist = 0 if previous primary particle scored
//     n_hist > 0 indicates how many primary particle read before the current one
//  2: LATCH (EGS)
//  3: ILB5 (PENELOPE) 
//  4: ILB4 (PENELOPE) 
//  5: ILB3 (PENELOPE) 
//  6: ILB2 (PENELOPE) 
//  7: ILB1 (PENELOPE) 
//  more to be defined

#define MAX_NUMB_EXTRAFLOAT_TYPES 3 /* maximum number of extra float allowed */
 //  0: User defined generic type
 //  1: XLAST (x coord. of the last interaction)  
 //  2: YLAST (y coord. of the last interaction)  
 //  3: ZLAST (z coord. of the last interaction)  
 //  more to be defined


struct iaea_header_type
{
  FILE *fheader;
  // ******************************************************************************
  // 1. PHSP format
  
  int file_type;            // 0 = phsp file ;  1 = phsp generator 
  int byte_order;           // as defined by get_byte_order routine
  int record_contents[9];   // record_contents[i] = 1 or 0 (variable or constant)
                            // correspond to the following logical variables :
                            //             ix,iy,iz,iu.iv,iw;
                            //       iweight,iextrafloat,iextralong;  

  float record_constant[7]; // if record_contents[i<7] = 0
                            // then record_constant[i] contents the constant value
                            // extra floats and longs are always variable
                            // so no need to store them

  // contains the keyword describing each stored extrafloat
  int extrafloat_contents[NUM_EXTRA_FLOAT];
  
  // contains the keyword describing each stored extralong
  int extralong_contents[NUM_EXTRA_LONG];   

  int record_length;
  //  record_length = 1 +                                     (particle)
  //                  4 +                                     (energy)
  //                  SUM(i=0;i<3) {record_contents[i]*4}     (ix,iy,iz)
  //                  SUM(i=3;i<6) {record_contents[i]*4}     (iu,iv,iw)
  //                  record_contents[6]*4 +                  (iweigth)
  //                  record_contents[7]*4 +                  (iextrafloat)
  //                  record_contents[8]*4 +                  (iextralong)
  IAEA_I64 checksum;

  // ******************************************************************************
  // 2. Mandatory description of the phsp
  
  char coordinate_system_description[MAX_STR_LEN*MAX_NUMB_LINES+1];

  // Counters for phsp file
  IAEA_I64 orig_histories;  
  IAEA_I64 nParticles;
  IAEA_I64 particle_number[MAX_NUM_PARTICLES];

  // Event generator input file
  char input_file_for_event_generator[MAX_STR_LEN*MAX_NUMB_LINES+1];
  
  // ******************************************************************************
  // 3. Mandatory additional information
  
  unsigned int iaea_index; // Agency ID
  char title[MAX_STR_LEN*MAX_NUMB_LINES+1];
  
  char machine_type[MAX_STR_LEN*MAX_NUMB_LINES+1];
  
  char MC_code_and_version[MAX_STR_LEN*MAX_NUMB_LINES+1];
  
  float global_photon_energy_cutoff;
  
  float global_particle_energy_cutoff;
  
  char transport_parameters[MAX_STR_LEN*MAX_NUMB_LINES+1];

  // ******************************************************************************
  // 4. Optional description
  
  char beam_name[MAX_STR_LEN*MAX_NUMB_LINES+1];
  
  char field_size[MAX_STR_LEN*MAX_NUMB_LINES+1];
  
  char nominal_SSD[MAX_STR_LEN*MAX_NUMB_LINES+1];
  
  char variance_reduction_techniques[MAX_STR_LEN*MAX_NUMB_LINES+1];
  
  char initial_source_description[MAX_STR_LEN*MAX_NUMB_LINES+1];

  // Documentation sub-section
  char MC_input_filename[MAX_STR_LEN*MAX_NUMB_LINES+1];
  
  // Assumed to be the preferred citation
  char published_reference[MAX_STR_LEN*MAX_NUMB_LINES+1]; 
  char authors[MAX_STR_LEN*MAX_NUMB_LINES+1];
  
  char institution[MAX_STR_LEN*MAX_NUMB_LINES+1];
  
  char link_validation[MAX_STR_LEN*MAX_NUMB_LINES+1];

  char additional_notes[MAX_STR_LEN*MAX_NUMB_LINES+1];

  // ******************************************************************************
  // 5. Optional statistical information
  double averageKineticEnergy[MAX_NUM_PARTICLES];
  double sumParticleWeight[MAX_NUM_PARTICLES];
  double maximumKineticEnergy[MAX_NUM_PARTICLES];
  double minimumKineticEnergy[MAX_NUM_PARTICLES];
  double minimumX, maximumX;
  double minimumY, maximumY;
  double minimumZ, maximumZ;  
  double minimumWeight[MAX_NUM_PARTICLES];
  double maximumWeight[MAX_NUM_PARTICLES];  

  IAEA_I64 read_indep_histories;  

// CLASS FUNCTIONS

public:
      int read_header();
      int write_header();
      int print_header();
      int set_record_contents(iaea_record_type *p_iaea_record);
      int get_record_contents(iaea_record_type *p_iaea_record);
      void initialize_counters();
      void update_counters(iaea_record_type *p_iaea_record);

private:
      int read_block(char *lineread, const char *blockname);
      int get_block(char *lineread);
      int get_blockname(char *line, const char *blockname);
      int write_blockname(const char *blockname);

      int check_byte_order();
      void print_statistics();
};

#endif
