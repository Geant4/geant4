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
/// \file biasing/ReverseMC01/src/Histo2Dvar.cc
/// \brief Implementation of the Histo2Dvar class
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              Histo2DVar.cc
//
// Version:                1.0
// Date:                09/03/00
// Author:                P R Truscott, F Lei
// Organisation:        DERA UK
// Customer:                ESA/ESTEC, NOORDWIJK
// Contract:                12115/96/NL/JG Work Order No. 3
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
// 28 August 1999, F Lei, DERA UK
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
#include "Histo2DVar.hh"
////////////////////////////////////////////////////////////////////////////////
//
Histo2DVar::Histo2DVar (G4String name,
      double *epx_list, size_t epx_list_len, side x_conv, 
      double *epy_list, size_t epy_list_len, side y_conv)
{
  //
  //
  // Set the name of the histogram, define the x and y VariableLengthPartitions
  // and reset the contents of the histogram.
  //
  theName = name;
  side x_conv_list[256] = {LEFT};
  size_t i;
  for (i=0; i<epx_list_len; i++) {x_conv_list[i]=x_conv;}
  x_part = VariableLengthPartition
    (epx_list, epx_list_len, x_conv_list, epx_list_len);

  side y_conv_list[256] = {LEFT};
  for (i=0; i<epy_list_len; i++) {y_conv_list[i]=y_conv;}
  y_part = VariableLengthPartition
    (epy_list, epy_list_len, y_conv_list, epy_list_len);

  reset();
  return;
}
////////////////////////////////////////////////////////////////////////////////
//
Histo2DVar::Histo2DVar ()
{
  //
  //
  // Set default name and partitions, and reset the contents of the histogram.
  //
  theName = "Blank Array";
  x_part = VariableLengthPartition ();
  y_part = VariableLengthPartition ();

  reset();
  return;
}
////////////////////////////////////////////////////////////////////////////////
//
void Histo2DVar::reset ()
{
  int i(0);
  int j(0);
  for (i=0; i<256; i++) {
    for (j=0; j<256; j++) {
      totalWeight[i][j]        = 0.;
      totalWeightSquared[i][j] = 0.;
      nEvents[i][j]            = 0;
    }
  }

  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      uoflowTotalWeight[i][j]        = 0.;
      uoflowTotalWeightSquared[i][j] = 0.;
      uoflownEvents[i][j]            = 0;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void Histo2DVar::fill (double x_point, double y_point, double weight = 1.0)
{
  //
  //
  // Determine the bin numbers for the points x_point and y_point.
  //
  int i = (x_part.get_elem_bin(&x_point));
  int j = (y_part.get_elem_bin(&y_point));

  if (i == BIN_UNDERFLOW || i == BIN_OVERFLOW || 
      j == BIN_UNDERFLOW || j == BIN_OVERFLOW) {

  //
  //
  // Underflow or overflow condition - modify uoflow variables.
  //
    HistSpecialBin xSpecialBin(inrange);
    if (i == BIN_UNDERFLOW) {xSpecialBin = underflow_bin;}
    else if (i == BIN_OVERFLOW) {xSpecialBin = overflow_bin;}

    HistSpecialBin ySpecialBin(inrange);
    if (j == BIN_UNDERFLOW) {ySpecialBin = underflow_bin;}
    else if (j == BIN_OVERFLOW) {ySpecialBin = overflow_bin;}

    uoflowTotalWeight[xSpecialBin][ySpecialBin]        += weight;
    uoflowTotalWeightSquared[xSpecialBin][ySpecialBin] += weight*weight;
    uoflownEvents[xSpecialBin][ySpecialBin]++;
  }
  else {
  //
  //
  // The data point referred to by x_point and y_point is within the
  // histogrammed area - modify conventional histogram variables.
  //    
    totalWeight[i][j]        += weight;
    totalWeightSquared[i][j] += weight*weight;
    nEvents[i][j]++;

    uoflowTotalWeight[inrange][inrange]        += weight;
    uoflowTotalWeightSquared[inrange][inrange] += weight*weight;
    uoflownEvents[inrange][inrange]++;

  }
  nAllEvents++;
}
////////////////////////////////////////////////////////////////////////////////
//
double Histo2DVar::get_bin_value (int i, int j)
{
  //
  //
  // Determine if i lies within the range of the histogram x-axis.
  //
  HistSpecialBin xSpecialBin(inrange);
  if (size_t(i) > (x_part.total_bins()-1)) {xSpecialBin = overflow_bin;}
  else if (i < 0) {xSpecialBin = underflow_bin;}

  //
  //
  // Determine if j lies within the range of the histogram y-axis.
  //
  HistSpecialBin ySpecialBin(inrange);
  if (size_t(j) > (y_part.total_bins()-1)) {ySpecialBin = overflow_bin;}
  else if (j < 0) {ySpecialBin = underflow_bin;}

  //
  //
  // If i and j are within range, output the conventional histogram variables.
  // Otherwise output uoflow variables.
  //
  double value(0.);
  if (xSpecialBin == inrange && ySpecialBin == inrange) {
    value = totalWeight[i][j];}
  else {
    value = uoflowTotalWeight[xSpecialBin][ySpecialBin];}

  return value;
}
////////////////////////////////////////////////////////////////////////////////
//
double Histo2DVar::get_bin_error (int i, int j)
{
  //
  //
  // Determine if i lies within the range of the histogram x-axis.
  //
  HistSpecialBin xSpecialBin(inrange);
  if (size_t(i) > (x_part.total_bins()-1)) {xSpecialBin = overflow_bin;}
  else if (i < 0) {xSpecialBin = underflow_bin;}

  //
  //
  // Determine if j lies within the range of the histogram y-axis.
  //
  HistSpecialBin ySpecialBin(inrange);
  if (size_t(j) > (y_part.total_bins()-1)) {ySpecialBin = overflow_bin;}
  else if (j < 0) {ySpecialBin = underflow_bin;}

  //
  //
  // If i and j are within range, output the conventional histogram variables.
  // Otherwise output uoflow variables.
  //
  double error(0.);
  if (xSpecialBin == inrange && ySpecialBin == inrange)
    //    {error = sqrt(totalWeightSquared[i][j] -
    //  totalWeight[i][j]*totalWeight[i][j]);}
    { 
      if (nEvents[i][j]>0) {
        error = totalWeight[i][j]/std::sqrt((G4double) nEvents[i][j]);
      }else {
        error = 0.;} }
  else
    {error = std::sqrt(uoflowTotalWeightSquared[xSpecialBin][ySpecialBin] -
                  uoflowTotalWeight[xSpecialBin][ySpecialBin]*
                  uoflowTotalWeight[xSpecialBin][ySpecialBin]);}
      return error;
}

////////////////////////////////////////////////////////////////////////////////
//
void Histo2DVar::div (double r)
{
  //
  //
  // Apply division to the conventional histogram variables.
  //
  int i(0);
  int j(0);
  for (i = 0; i < 3; i ++)
  {
    for (j = 0; j < 3; j++)
    {
      uoflowTotalWeight[i][j]        = uoflowTotalWeight[i][j] / r;
      uoflowTotalWeightSquared[i][j] = uoflowTotalWeightSquared[i][j] / r;
    }
  }

  //
  //
  // Apply division to the uoflow variables.
  //
  for (i = 0; i < int(x_part.total_bins()); i ++)
  {
    for (j = 0; j < int(y_part.total_bins()); j++)
    {
      totalWeight[i][j]        = totalWeight[i][j] / r;
      totalWeightSquared[i][j] = totalWeightSquared[i][j] / r /r;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
