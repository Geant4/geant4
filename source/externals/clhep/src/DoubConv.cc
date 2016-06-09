// -*- C++ -*-
// $Id:$
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- DoubConv ---
//                      class implementation file
// -----------------------------------------------------------------------

#include "CLHEP/Random/DoubConv.h"

#include <sstream>
#include <iomanip>

namespace CLHEP {

bool DoubConv::byte_order_known = false;
int  DoubConv::byte_order[8];

void DoubConv::fill_byte_order () {
  double x = 1.0;
  int t30 = 1 << 30;
  int t22 = 1 << 22;
  x *= t30;
  x *= t22;
  double y = 1;
  double z = 1;
  x *= z;
  for (int k=0; k<6; k++) {
    x += y*z;
    y += 1;
    z *= 256;
  }
  // x, in IEEE format, would now be 0x4330060504030201
  union DB8 {
    unsigned char b[8];
    double d;
  };
  DB8 xb;
  xb.d = x;
  int n;
  static const int UNSET = -1; 
  for (n=0; n<8; n++) {
    byte_order[n] = UNSET;
  }
  int order;
  for (n=0; n<8; n++) {
    switch ( xb.b[n] ) {
      case 0x43:
        order = 0;
        break;
      case 0x30:
        order = 1;
        break;
      case 0x06:
        order = 2;
        break;
      case 0x05:
        order = 3;
        break;
      case 0x04:
        order = 4;
        break;
      case 0x03:
        order = 5;
        break;
      case 0x02:
        order = 6;
        break;
      case 0x01:
        order = 7;
        break;
      default:
        throw DoubConvException(
		"Cannot determine byte-ordering of doubles on this system");
    } 
    if (byte_order[n] != UNSET) {
        throw DoubConvException(
		"Confusion in byte-ordering of doubles on this system");
    }    
    byte_order[n] = order;
    byte_order_known = true;
  }
  return;
}

std::string DoubConv::d2x(double d) {
  if ( !byte_order_known ) fill_byte_order ();
  DB8 db;
  db.d = d;
  std::ostringstream ss;
  for (int i=0; i<8; ++i) {
    int k = byte_order[i];
    ss << std::hex << std::setw(2) << std::setfill('0') << (int)db.b[k];
  }
  return ss.str();
}

std::vector<unsigned long> DoubConv::dto2longs(double d) {
  std::vector<unsigned long> v(2);
  if ( !byte_order_known ) fill_byte_order ();
  DB8 db;
  db.d = d;
  v[0] =   ((static_cast<unsigned long>(db.b[byte_order[0]])) << 24)
         | ((static_cast<unsigned long>(db.b[byte_order[1]])) << 16)
         | ((static_cast<unsigned long>(db.b[byte_order[2]])) <<  8)
         | ((static_cast<unsigned long>(db.b[byte_order[3]]))      );
  v[1] =   ((static_cast<unsigned long>(db.b[byte_order[4]])) << 24)
         | ((static_cast<unsigned long>(db.b[byte_order[5]])) << 16)
         | ((static_cast<unsigned long>(db.b[byte_order[6]])) <<  8)
         | ((static_cast<unsigned long>(db.b[byte_order[7]]))      );
  return v; 
}

double DoubConv::longs2double (const std::vector<unsigned long> & v) {
  DB8 db;
  unsigned char bytes[8];
  if ( !byte_order_known ) fill_byte_order ();
  bytes[0] = static_cast<unsigned char>((v[0] >> 24) & 0xFF);
  bytes[1] = static_cast<unsigned char>((v[0] >> 16) & 0xFF);
  bytes[2] = static_cast<unsigned char>((v[0] >>  8) & 0xFF);
  bytes[3] = static_cast<unsigned char>((v[0]      ) & 0xFF);
  bytes[4] = static_cast<unsigned char>((v[1] >> 24) & 0xFF);
  bytes[5] = static_cast<unsigned char>((v[1] >> 16) & 0xFF);
  bytes[6] = static_cast<unsigned char>((v[1] >>  8) & 0xFF);
  bytes[7] = static_cast<unsigned char>((v[1]      ) & 0xFF);
  for (int i=0; i<8; ++i) {
    db.b[byte_order[i]] =  bytes[i];
  }  
  return db.d;
}

} // end namespace CLHEP
