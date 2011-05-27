// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                      --- engineIDulong ---
//                      function implementation file
// -----------------------------------------------------------------------
//
// =======================================================================
// Mark Fischler  - Created: Mar. 8, 2005
// =======================================================================

#include <string>
#include <vector>

namespace CLHEP {

static std::vector<unsigned long> gen_crc_table() {
  /* generate the table of CRC remainders for all possible bytes */
  static const unsigned long POLYNOMIAL = 0x04c11db7UL;
  std::vector<unsigned long> crc_table;
  for ( unsigned long i = 0;  i < 256;  ++i ) {
    unsigned long crc = i << 24;
    for ( int j = 0;  j < 8;  j++ ) {
      if ( crc & 0x80000000UL ) {
        crc = ( ( crc << 1 ) ^ POLYNOMIAL ) & 0xffffffffUL;
      } else {
        crc = ( crc << 1 ) & 0xffffffffUL; 
      }
    }
    crc_table.push_back(crc);
  }
  return crc_table;
}

unsigned long crc32ul(const std::string & s) {
  static std::vector<unsigned long> crc_table =  gen_crc_table();
  unsigned long crc = 0;
  int end = s.length();
  for (int j = 0; j != end; ++j) {
    int i = ( (int) ( crc >> 24) ^ s[j] ) & 0xff;
    crc = ( ( crc << 8 ) ^ crc_table[i] ) & 0xffffffffUL;
  }
  return crc;
}

}  // namespace CLHEP

