////////////////////////////////////////////////////////////////////////////////
//
#include "MLHisto1D.hh"
#include <iomanip>
////////////////////////////////////////////////////////////////////////////////
//
CSVofstream & operator << (CSVofstream &s, MLHisto1D &q)
//
//
// Definition of the insertion operator << to provide the histogram output to
// CSVofstream.
//
{
//  s.setf(ios::scientific);
//  s.precision(3);
  s <<std::setiosflags(std::ios::scientific);
  s <<std::setprecision(4);
  for (G4int i=0; i < G4int(q.part.total_bins()); i++) {
    s <<" "
      <<std::setw(10) <<q.part.get_bin_position(i) <<", "
      <<std::setw(10) <<q.part.get_bin_position(i+1) <<", "
      <<std::setw(10) <<q.get_bin_position(i) <<", "
      <<std::setw(10) <<q.get_bin_value(i) *
      q.GetNormalisation() <<", "
      <<std::setw(10) <<q.get_bin_error(i) *
      q.GetNormalisation() <<G4endl;
  }
  return s;
}
////////////////////////////////////////////////////////////////////////////////
//
RPTofstream & operator << (RPTofstream &s, MLHisto1D &q)
//
//
// Definition of the insertion operator << to provide the histogram output to
// RPTofstream.
//
{
  s.precision (3);
  s <<q.get_name() <<G4endl;
  s <<"   Bin Range                Mean bin position   Bin value      Error  "
    <<G4endl;
  s <<"    ["<<q.get_unitx() <<"]" <<"                        ["
    << q.get_unitx() <<"]  "<<"        [" << q.get_unity() <<"]"<< G4endl;
  s <<std::setw(28) <<"UNDERFLOW";
  s.setf (std::ios::scientific);
  s <<std::setw(14) <<q.get_bin_position(-1) << " ";
  s <<std::setw(12) <<q.get_bin_value(underflow_bin) *
    q.GetNormalisation() << " ";
  s <<std::setw(12) <<q.get_bin_error(underflow_bin) *
    q.GetNormalisation() <<G4endl;
  for (G4int i=0; i < G4int(q.part.total_bins()); i++) {
    s <<std::setw(12) <<q.part.get_bin_position(i) <<" "
      <<std::setw(12) <<q.part.get_bin_position(i+1) <<" "
      <<std::setw(14) <<q.get_bin_position(i) << " "
      <<std::setw(12) <<q.get_bin_value(i) *
      q.GetNormalisation() <<" "
      <<std::setw(12) <<q.get_bin_error(i) *
      q.GetNormalisation() <<G4endl;
  }
  s <<std::setw(28) <<"OVERFLOW";
  s <<std::setw(14) <<q.get_bin_position(1000000) << " ";
  s <<std::setw(12) <<q.get_bin_value(overflow_bin) *
    q.GetNormalisation() <<" ";
  s <<std::setw(12) <<q.get_bin_error(overflow_bin) *
    q.GetNormalisation() <<G4endl;
  return s;
}
////////////////////////////////////////////////////////////////////////////////
