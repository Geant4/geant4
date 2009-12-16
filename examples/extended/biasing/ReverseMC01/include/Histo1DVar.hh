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
#ifndef Histo1DVar_h
#define Histo1DVar_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              Histo1DVar.hh
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
// DESCRIPTION
// -----------
//
// The Histo1DVar class allows definition of a 1-D binning scheme and the
// storage of events according to that scheme.  Member functions permit output
// of the bin contents and standard deviation of the events recorded for each
// bin.  The member functions for this class are only those required to meet
// the needs of the SSAT.  However, where possible, the data-types of the
// arguments and returned values are identical to those used in the Histoograms
// category of LHC++.  Information associated with the binning scheme is
// contained in a VariableLengthPartition object.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
//
// Histo1DVar (G4String name, double *ep_list, 
//             size_t ep_list_len, side conv)
//    Constructor defining the name and binning scheme for the histogram:
//    1. name is the name of the histogram;
//    2. ep_list is a pointer to a double array defining the bin edges;
//    3. ep_list_len is the number of values in ep_list;
//    4. conv the convention of data-points exactly on a bin-edge (if the
//       datapoint is associated with the bin to the left of the bin-edge then
//       conv==LEFT, otherwise conv==RIGHT).
//
// Histo1DVar ()
//    Default constructor.
//
// ~Histo1DVar ()
//    Destructor.
//
// void reset ()
//    Zeros-out the contents of the histogram (binning scheme is unaffected).
//
// void fill (double data_point, double weight = 1.0)
//    Increments the bin associated with the datapoint data_point by the amount
//    weight.
//
// double get_bin_value (HistSpecialBin specialBin)
//    Returns the contents of special bin specialBin (specialBin may have
//    values underflow_bin, inrange or overflow_bin).
//
// double get_bin_error (HistSpecialBin specialBin)
//    Returns the error associated with special bin specialBin (specialBin may
//    have values underflow_bin, inrange or overflow_bin).
// 
// double get_bin_value (int i)
//    Returns the contents of bin i (if i is less than 0 then the underflow_bin
//    is returned; if i is greater than the number of bins then the
//    overflow_bin is returned).
//
// double get_bin_error (int i)
//    Returns the error associated with bin i (if i is less than 0 then the
//    error associated with the underflow_bin is returned; if i is greater than
//    the number of bins then the error associated with the overflow_bin is
//    returned).
//
// double get_all_bins ()
//    Returns the integral under the whole histogram, including underflow_bin
//    and overflow_bin.
//
// void div (double r)
//    Divides the contents of all bins by the value r.
//
// VariableLengthPartition part;
//    Accesses the VariableLengthPartition part defining the histogram binning
//    scheme.
//
// void set_name (G4String new_name)
//    Sets the name of the histogram to new_name.
//
// char const *get_name ()
//    Returns the name of the histogram.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 30 June 1999, P R Truscott, DERA UK
// Version number update 0.b.2 -> 0.b.3, but no functional change.
//
//
// 28 August 1999, F Lei & P R Truscott, DERA UK
// Version number update 0.b.3 -> 0.b.4, but no functional change. 
//
// 17 September 1999, P R Truscott, DERA UK
// Version number update 0.b.4 -> 0.b.5, but no functional change.
//
// 09 March 2000, P R Truscott, DERA UK
// Update 0.b.3 -> 1.0, for compliance with ISO ANSI C++ (no functional change).
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"

#include "Histograms.hh"
#include "VariableLengthPartition.hh"
#include "side.hh"
#include <vector>
//////////////////////////////////////////////////////////////////////////////
//
class Histo1DVar
{
public:
  Histo1DVar (G4String name, double *ep_list, size_t ep_list_len, side conv);
  Histo1DVar ();
  ~Histo1DVar () {};
  void reset ();
  void fill (double data_point, double weight);
  double get_bin_value (HistSpecialBin specialBin);
  double get_bin_error (HistSpecialBin specialBin);
  double get_bin_value (int );
  double get_bin_error (int);
  double get_bin_position (int);
  double get_all_bins ();
  void div (double r);

  VariableLengthPartition part;

private:
  G4String theName;
  size_t length;
  std::vector<double>  totalWeight;
  std::vector<double>  meanPosition;
  double   overflowTotalWeight;
  double   underflowTotalWeight;
  double   meanOverflowPosition;
  double   meanUnderflowPosition;
  std::vector<double>   totalWeightSquared;
  double   overflowTotalWeightSquared;
  double   underflowTotalWeightSquared;
  std::vector<long>     nEvents;
  long     overflownEvents;
  long     underflownEvents;
  long     nAllEvents;

  //
  //
  // INLINE DECLARATIONS/DEFINTIONS:
  //
public:
  inline void set_name (G4String new_name) {theName = new_name;}
  inline const char *get_name () {return theName.data();}
};
////////////////////////////////////////////////////////////////////////////////
#endif





