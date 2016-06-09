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
#ifndef Histo2DVar_h
#define Histo2DVar_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              Histo2DVar.hh
//
// Version:		1.0
// Date:		09/03/00
// Author:		P R Truscott, F Lei
// Organisation:	DERA UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		12115/96/NL/JG Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------
//
// The Histo2DVar class allows definition of a 2-D binning scheme (i.e. using
// two, 1-D binning schemes) and the storage of events according to the binning
// scheme.  Member functions permit output of the bin contents and standard
// deviation of the events recorded for each bin.  The member functions for
// this class are only those required to meet the needs of the SSAT.  However,
// where possible the data-types of the arguments and returned values are
// identical to those used in the Histoograms category of LHC++.  Information
// associated with the binning scheme is contained in VariableLengthPartition
// object.
//
// Since the maximum length of a VariableLengthPartition is 256, 255 bins x255
// bins is the maximum number size for a Histo2DVar object.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
//
// Histo2DVar (G4String name, double *epx_list, size_t epx_list_len,
// side x_conv, double *epy_list, size_t epy_list_len, side y_conv)
//    Constructor defining:
//    1. name is the name of the histogram;
//    2. epx_list is a pointer to a double array defining the x-coordinate bin
//       edges;
//    3. epx_list_len is the number of values in epx_list;
//    4. x_conv the convention of data-points exactly on a bin-edge for the
//       x-coordinate (if the datapoint is associated with the bin to the left
//       of the bin-edge then conv==LEFT, otherwise conv==RIGHT);
//    5. epy_list is a pointer to a double array defining the y-coordinate bin
//       edges;
//    6. epy_list_len is the number of values in epy_list;
//    7. y_conv the convention of data-points exactly on a bin-edge for the
//       y-coordinate.
//
// Histo2DVar ()
//    Default constructor.
//
// ~Histo2DVar ()
//    Destructor.
//
// void reset ()
//    Zeros-out the contents of the histogram (binning scheme is unaffected).
//
// void fill (double x_point, double y_point, double weight = 1.0)
//    Increments the bin associated with the datapoint (x_point,y_point) by
//    the amount weight.
//
// double get_bin_value (HistSpecialBin xSpecialBin,
// HistSpecialBin ySpecialBin)
//    Returns the contents of special bin (xSpecialBin,ySpecialBin) 
//    (xSpecialBin and ySpeciaBin may have values underflow_bin, inrange or
//    overflow_bin).
//
// double get_bin_error (HistSpecialBin xSpecialBin,
// HistSpecialBin ySpecialBin)
//    Returns the error associated with special bin (xSpecialBin,ySpecialBin)
//    (xSpecialBin and ySpeciaBin may have values underflow_bin, inrange or
//    overflow_bin).
//
// double get_bin_value (int i, int j)
//    Returns the contents of bin (i,j).  If i is less than 0 then the value
//    returned is for x==underflow_bin, y==underflow_bin, inrange, or
//    overflow_bin, depending upon whether j is less than 0, in-range or
//    greater than the number of y-bins.  If i is greater than the number of
//    x-bins then the value returned is for x==overflow_bin, y==underflow_bin,
//    inrange, or overflow_bin, depending upon whether j is less than 0,
//    in-range or greater than the number of y-bins.
//
// double get_bin_error (int i, int j)
//    Returns the error associated with bin i.  If i is less than 0 then the
//    value returned is for x==underflow_bin, y==underflow_bin, inrange, or
//    overflow_bin, depending upon whether j is less than 0, in-range or
//    greater than the number of y-bins.  If i is greater than the number of
//    x-bins then the value returned is for x==overflow_bin, y==underflow_bin,
//    inrange, or overflow_bin, depending upon whether j is less than 0,
//    in-range or greater than the number of y-bins.
//
// double get_all_bins ()
//    Returns the integral under the whole histogram, including the 
//    over-/underflow bins.
//
// void div (double r)
//    Divides the contents of all bins by the value r.
//
// VariableLengthPartition x_part;
//    Accesses the VariableLengthPartition x_part defining the histogram
//    x-binning scheme.
//
// VariableLengthPartition y_part;
//    Accesses the VariableLengthPartition y_part defining the histogram
//    y-binning scheme.
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
////////////////////////////////////////////////////////////////////////////////
//
class Histo2DVar
{
  public:
    Histo2DVar (G4String name,
      double *epx_list, size_t epx_list_len, side x_conv, 
      double *epy_list, size_t epy_list_len, side y_conv);
    Histo2DVar ();
    ~Histo2DVar () {};
    void reset ();
    void fill (double x_point, double y_point, double weight);
    double get_bin_value (int , int );
    double get_bin_error (int , int );
    void div (double );

    VariableLengthPartition x_part;
    VariableLengthPartition y_part;

  private:
    G4String theName;
    double   totalWeight[256][256];
    double   uoflowTotalWeight[3][3];
    double   totalWeightSquared[256][256];
    double   uoflowTotalWeightSquared[3][3];
    long     nEvents[256][256];
    long     uoflownEvents[3][3];
    long     nAllEvents;

  //
  //
  // INLINE DECLARATIONS/DEFINTIONS:
  //
  public:
    inline void set_name (G4String new_name) {theName = new_name;}
    inline const char *get_name () {return theName.data();}
    inline double get_bin_value
      (HistSpecialBin xSpecialBin, HistSpecialBin ySpecialBin)
      {return uoflowTotalWeight[xSpecialBin][ySpecialBin];}
    inline double get_bin_error
      (HistSpecialBin xSpecialBin, HistSpecialBin ySpecialBin)
      {return std::sqrt(uoflowTotalWeightSquared[xSpecialBin][ySpecialBin]);}
    inline double get_all_bins () const
      {return uoflowTotalWeight[inrange][inrange];}
};
////////////////////////////////////////////////////////////////////////////////
#endif



