//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file biasing/ReverseMC01/include/VariableLengthPartition.hh
/// \brief Definition of the VariableLengthPartition class
//
#ifndef VariableLengthPartition_h
#define VariableLengthPartition_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              VariableLengthParition.hh
//
// Version:                1.0
// Date:                09/03/00
// Author:                P R Truscott
// Organisation:        DERA UK
// Customer:                ESA/ESTEC, NOORDWIJK
// Contract:                12115/96/NL/JG Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------
//
// This class is used to emulate the VariableLengthPartition class in LHC++
// Histoograms.  It is used to contain an arbitrary binnng scheme.  Member
// functions provide information on the locations of the bin-edges, and the
// bin within which a user-provided G4double falls.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
//
// VariableLengthPartition (double *ep_list, size_t ep_list_len,
// side *ep_conv, size_t ep_conv_len)
//    Constructor defining:
//    1. ep_list is a pointer to a double array defining the bin edges;
//    2. ep_list_len is the number of values in ep_list;
//    3. ep_conv the convention of data-points exactly on a bin-edge (if the
//       datapoint is associated with the bin to the left of the bin-edge
//       then ep_conv[i]==LEFT, otherwise ep_conv[i]==RIGHT);
//    4. ep_conv_len is the number of values in ep_conv_len.
//
// VariableLengthPartition ()
//    Default constructor.
//
// ~VariableLengthPartition ()
//    Destructor.
//
// const long get_elem_bin (double *data_point)
//    Returns the bin number within which the value *data_point falls.  If
//    *data_point falls below the lower edge of the lowest bin, then
//    BIN_UNDERFLOW is returned; if *data_point exceeds the upper edge of the
//    highest bin then BIN_OVERFLOW is returned (see Histograms).
//
// double get_bin_position (size_t bin_id)
//    Returns the lower-edge of the bin bin_id.  If bin_id is equal to the
//    number of bins, then the upper edge of the last bin is returned.
//
// size_t total_bins ()
//    Returns the number of bins in the binning scheme.
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
#include "side.hh"
#include "Histograms.hh"
#include "globals.hh"
#include <vector>
//
class VariableLengthPartition
{
  public:
    VariableLengthPartition
      (double *ep_list, size_t ep_list_len,
       side *ep_conv, size_t ep_conv_len);
    VariableLengthPartition ();
    ~VariableLengthPartition () {};
    long get_elem_bin (double *data_point);
    double get_bin_position (size_t bin_id);
    size_t total_bins () {return nbin;}

  private:
    size_t nbin;
    side   conv;
    long   location;

    std::vector<double> bin;
}
;
////////////////////////////////////////////////////////////////////////////////
#endif
