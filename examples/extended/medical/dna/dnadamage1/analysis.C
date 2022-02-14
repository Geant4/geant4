// -------------------------------------------------------------------
// $Id: 
// -------------------------------------------------------------------
//
// *********************************************************************
// To execute this macro under ROOT after your simulation ended, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X analysis.C' at the ROOT session prompt
// *********************************************************************


{
    gROOT->ProcessLine(".L molecule.C");
    gROOT->ProcessLine(".x plot.C");
}
