// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              VariableLengthParition.cc
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
// CHANGE HISTORY
// --------------
//
// 30 June 1999, P R Truscott, DERA UK
// Version number update 0.b.2 -> 0.b.3, but no functional change.
//
// 28 August 1999, F Lei, DERA UK
// Version number update 0.b.3 -> 0.b.4, but no functional change.
//
// 17 September 1999, F Lei & P R Truscott, DERA UK
// Version 0.b.5
// Now uses STL rather than RW vectors.
//
// 09 March 2000, P R Truscott, DERA UK
// Update 0.b.3 -> 1.0, for compliance with ISO ANSI C++.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "VariableLengthPartition.hh"

#include <iostream.h>
////////////////////////////////////////////////////////////////////////////////
//
VariableLengthPartition::VariableLengthPartition
  (double *ep_list, size_t ep_list_len, 
   side *ep_conv, size_t )
{
  //
  //
  // Transfer over bin-edge information to the STL vector.
  //
  bin.insert(bin.begin(), ep_list, ep_list+ep_list_len);
  nbin = bin.size()-1; 
  conv = *ep_conv;
  return;
}
////////////////////////////////////////////////////////////////////////////////
//
VariableLengthPartition::VariableLengthPartition ()
{
  conv = LEFT;
  return;
}
////////////////////////////////////////////////////////////////////////////////
//
const long VariableLengthPartition::get_elem_bin (double *data_point)
{
  //
  //
  // Find out whether *data_point falls below of exceeds the range of the
  // binning scheme.
  //
  int location (0);
  if ((conv == LEFT && *data_point > bin[nbin]) ||
      (conv == RIGHT && *data_point >= bin[nbin])) {
    location = BIN_OVERFLOW;}
  else if ((conv == LEFT && *data_point <= *bin.begin()) ||
      (conv == RIGHT && *data_point < *bin.begin())) {
    location = BIN_UNDERFLOW;}
  else {
  //
  //
  // *data_point is in range - perform a binary search.
  //
    int upperLimit = bin.size()-1;
    int lowerLimit = 0;
    location = (upperLimit-1)/2;
    while (*data_point < bin[location] || *data_point > bin[location+1]) {
      if (*data_point < bin[location]) {upperLimit = location;}
      else {lowerLimit = location;}
      location = (upperLimit+lowerLimit)/2;
    }
    if (bin[location] == *data_point && conv == RIGHT) {
      location++;}
  }
  return location;
}
////////////////////////////////////////////////////////////////////////////////
//
double VariableLengthPartition::get_bin_position (size_t bin_id)
{
  return bin[bin_id];
}
////////////////////////////////////////////////////////////////////////////////


