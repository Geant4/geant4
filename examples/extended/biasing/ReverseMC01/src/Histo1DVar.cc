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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              Histo1DVar.cc
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
#include "Histo1DVar.hh"
////////////////////////////////////////////////////////////////////////////////
//
Histo1DVar::Histo1DVar (G4String name, double *ep_list, size_t ep_list_len,
  side conv = LEFT)
  : theName(name)
{
  //
  // Set the name of the histogram, define the VariableLengthPartition and
  // reset the contents of the histogram.
  //
  side conv_list[1] = {conv};
  part = VariableLengthPartition
    (ep_list, ep_list_len, conv_list, 1);
  length = part.total_bins();
  reset();
  return;
}
////////////////////////////////////////////////////////////////////////////////
//
Histo1DVar::Histo1DVar ()
{
  //
  //
  // Set default name and partition, and reset the contents of the histogram.
  //
  theName = "Blank Array";
  part = VariableLengthPartition ();
  length= 0;
  reset();
  return;
}
////////////////////////////////////////////////////////////////////////////////
//
void Histo1DVar::reset ()
{
  totalWeight.clear();
  meanPosition.clear();
  totalWeightSquared.clear();
  nEvents.clear();
  for (size_t i=0; i< length; i++) {
    totalWeight.push_back(0.);
    meanPosition.push_back(0.);
    totalWeightSquared.push_back(0.);
    nEvents.push_back(0);
  }
  underflowTotalWeight         = 0.;
  overflowTotalWeight          = 0.;

  meanUnderflowPosition         = 0.;
  meanOverflowPosition          = 0.;
  
  underflowTotalWeightSquared  = 0.;
  overflowTotalWeightSquared   = 0.;
  
  underflownEvents             = 0;
  overflownEvents              = 0;
  
  nAllEvents                   = 0;
}
////////////////////////////////////////////////////////////////////////////////
//
void Histo1DVar::fill (double data_point, double weight = 1.0)
{
  //
  //
  // Determine the bin numbers for the point data_point.
  //
  int i = (part.get_elem_bin(&data_point));

  switch (i) {

  //
  //
  // If an underflow or overflow condition is present, then modify the overflow
  // or underflow variables, otherwise modify the conventional histogram
  // variables.
  //
    case BIN_OVERFLOW : 
      overflowTotalWeight        += weight;
      overflowTotalWeightSquared += weight*weight;
      meanOverflowPosition += data_point*weight; // saved as total for effiency 
      overflownEvents++;
      break;
    case BIN_UNDERFLOW :
      underflowTotalWeight        += weight;
      underflowTotalWeightSquared += weight*weight;
      meanUnderflowPosition += data_point*weight; // saved as total for effiency 
      underflownEvents++;
      break;
    default:
      totalWeight[i]        += weight;
      totalWeightSquared[i] += weight*weight;
      meanPosition[i] += data_point*weight; // saved as total for effiency
      nEvents[i]++;
  }
  nAllEvents++;
}
////////////////////////////////////////////////////////////////////////////////
//
double Histo1DVar::get_bin_value (HistSpecialBin specialBin)
{
  double value(0.);
  switch (specialBin) {

  //
  //
  // Output the contents of the overflow, underflow or inrange bin depending
  // upon the value of specialBin.
  //
    case overflow_bin :
      value = overflowTotalWeight;
      break;
    case underflow_bin :
      value = underflowTotalWeight;
      break;
    case inrange :
      value = get_all_bins();
      break;
  }
  return value;
}
////////////////////////////////////////////////////////////////////////////////
//
double Histo1DVar::get_bin_error (HistSpecialBin specialBin)
{
  double error(0.);
  switch (specialBin) {

  //
  //
  // Output the error of the overflow, underflow or inrange bin depending
  // upon the value of specialBin.
  //
    case overflow_bin :
      if (overflownEvents>0) {
        error = overflowTotalWeight/std::sqrt((G4double) overflownEvents);}
      else {
        error = 0.;}
      break;
    case underflow_bin :
      if (underflownEvents>0) {
        error = underflowTotalWeight/std::sqrt((G4double) underflownEvents);}
      else {
        error = 0.;}
      break;
    case inrange : 
      error = 0.;
      if (nAllEvents-overflownEvents-underflownEvents>0) {
        for (size_t i=0; i<(part.total_bins()); i++) {
          error+=totalWeightSquared[i];}
        error = error/std::sqrt((G4double) nAllEvents-overflownEvents-underflownEvents);
      }
      break;
  }
  return error;
}
////////////////////////////////////////////////////////////////////////////////
//
double Histo1DVar::get_bin_value (int i)
{ 
  //
  //
  // If i is within range, output the conventional histogram variables.
  // Otherwise output overflow or underflow variables.
  //
  double value(0.);
  if (i > int(part.total_bins()-1)) 
    {value = overflowTotalWeight;}
  else if (i < 0) 
    {value = underflowTotalWeight;}
  else 
    {value = totalWeight[i];}
  //  cout << i << " " << value << " " << totalWeight[i] << endl;
  return value;
}
////////////////////////////////////////////////////////////////////////////////
//
double Histo1DVar::get_bin_error (int i)
{
  //
  //
  // If i is within range, output the conventional histogram variables.
  // Otherwise output overflow or underflow variables.
  //
  double error(0.);
  if (i > int(part.total_bins())-1) {
    if (overflownEvents>0) {
      error = overflowTotalWeight/std::sqrt((G4double) overflownEvents);}
    else {
      error = 0.;}
  }
  else if (i < 0) {
    if (underflownEvents>0) {
      error = underflowTotalWeight/std::sqrt((G4double) underflownEvents);}
    else {
      error = 0.;}
  }
  else {
    if (nEvents[i]>0) {
      error = totalWeight[i]/std::sqrt((G4double) nEvents[i]);}
    else {
      error = 0.;}
  }
  return error;
}

////////////////////////////////////////////////////////////////////////////////
//
double Histo1DVar::get_bin_position (int i)
{ 
  //
  //
  // If i is within range, output the conventional histogram variables.
  // Otherwise output overflow or underflow variables.
  //
  double value = 0.;
  if (i > int((part.total_bins()-1))) {
    if (overflownEvents > 0 ) value = meanOverflowPosition/overflowTotalWeight;
  }
  else if (i < 0) {
    if (underflownEvents > 0 ) value = meanUnderflowPosition/underflowTotalWeight;    
  }
  else if (nEvents[i] > 0 ) value = meanPosition[i]/totalWeight[i];
  
  return value;
}

////////////////////////////////////////////////////////////////////////////////
//
double Histo1DVar::get_all_bins ()
{
  //
  //
  // Sum up all bins, including overflow and underflow.
  //
  double sum = underflowTotalWeight + overflowTotalWeight;
  for (size_t i=0; i<(part.total_bins()); i++) {sum += totalWeight[i];}

  return sum;
}
////////////////////////////////////////////////////////////////////////////////
//
void Histo1DVar::div (double r)
{
  //
  //
  // Divide all histogram (including underflow and overflow) variables by r.
  //
  overflowTotalWeight         = overflowTotalWeight / r;
  overflowTotalWeightSquared  = overflowTotalWeightSquared / r;
  meanOverflowPosition        = meanOverflowPosition / r;
  underflowTotalWeight        = underflowTotalWeight / r;
  underflowTotalWeightSquared = underflowTotalWeightSquared / r;
  meanUnderflowPosition       = meanUnderflowPosition / r;
  for (size_t i = 0; i<(part.total_bins()); i++)
  {
    totalWeight[i]        = totalWeight[i] / r;
    totalWeightSquared[i] = totalWeightSquared[i] / r;
    meanPosition[i] = meanPosition[i]/r;
  }
}
////////////////////////////////////////////////////////////////////////////////











