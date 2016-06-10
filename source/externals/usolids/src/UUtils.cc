//
// ********************************************************************
// * This Software is part of the AIDA Unified Solids Library package *
// * See: https://aidasoft.web.cern.ch/USolids                        *
// ********************************************************************
//
// $Id:$
//
// --------------------------------------------------------------------
//
// UUtils
//
// 19.10.12 Marek Gayer
// --------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "UVector3.hh"
#include "UTransform3D.hh"
#include "UUtils.hh"
#include "VUSolid.hh"

using namespace std;

//______________________________________________________________________________
void UUtils::TransformLimits(UVector3& min, UVector3& max, const UTransform3D& transformation)
{
  // The goal of this method is to convert the quantities min and max (representing the
  // bounding box of a given solid in its local frame) to the main frame, using
  // "transformation"
  UVector3 vertices[8] =     // Detemination of the vertices thanks to the extension of each solid:
  {
    UVector3(min.x(), min.y(), min.z()), // 1st vertice:
    UVector3(min.x(), max.y(), min.z()), // 2nd vertice:
    UVector3(max.x(), max.y(), min.z()),
    UVector3(max.x(), min.y(), min.z()),
    UVector3(min.x(), min.y(), max.z()),
    UVector3(min.x(), max.y(), max.z()),
    UVector3(max.x(), max.y(), max.z()),
    UVector3(max.x(), min.y(), max.z())
  };

  min.Set(kInfinity);
  max.Set(-kInfinity);

  // Loop on th vertices
  int limit = sizeof(vertices) / sizeof(UVector3);
  for (int i = 0 ; i < limit; i++)
  {
    // From local frame to the gobal one:
    // Current positions on the three axis:
    UVector3 current = transformation.GlobalPoint(vertices[i]);

    // If need be, replacement of the min & max values:
    if (current.x() > max.x()) max.x() = current.x();
    if (current.x() < min.x()) min.x() = current.x();

    if (current.y() > max.y()) max.y() = current.y();
    if (current.y() < min.y()) min.y() = current.y();

    if (current.z() > max.z()) max.z() = current.z();
    if (current.z() < min.z()) min.z() = current.z();
  }
}

double UUtils::Random(double min, double max)
{
  // srand((unsigned)time(NULL));
  double number = (double) rand() / RAND_MAX;
  double res = min + number * (max - min);
  return res;
}

string UUtils::ToString(int number)
{
  std::stringstream ss;
  ss << number;
  return ss.str();
}

string UUtils::ToString(double number)
{
  std::stringstream ss;
  ss << number;
  return ss.str();
}

int UUtils::FileSize(const std::string& filePath)
{
  std::streampos fsize = 0;
  std::ifstream file(filePath.c_str(), std::ios::binary);

  fsize = file.tellg();
  file.seekg(0, std::ios::end);
  fsize = file.tellg() - fsize;
  file.close();

  return fsize;
}


int UUtils::StrPos(const string& haystack, const string& needle)
{
  int sleng = haystack.length();
  int nleng = needle.length();

  if (sleng == 0 || nleng == 0)
    return -1;

  for (int i = 0, j = 0; i < sleng; j = 0, i++)
  {
    while (i + j < sleng && j < nleng && haystack[i + j] == needle[j])
      j++;
    if (j == nleng)
      return i;
  }
  return -1;
}
void UUtils:: Exception(const char* originOfException,
                        const char* exceptionCode,
                        UExceptionSeverity severity,
                        int level,
                        const char* description)

{
  bool toBeAborted = true;
  static const std::string es_banner
    = "\n-------- EEEE ------- UException-START -------- EEEE -------\n";
  static const std::string ee_banner
    = "\n-------- EEEE ------- UException-END --------- EEEE -------\n";
  static const std::string ws_banner
    = "\n-------- WWWW ------- UException-START -------- WWWW -------\n";
  static const std::string we_banner
    = "\n-------- WWWW -------- UException-END --------- WWWW -------\n";
  std::ostringstream message;
  message << "\n*** ExceptionHandler is not defined ***\n"
          << "*** Exception : " << exceptionCode << std::endl
          << "      issued by : " << originOfException << std::endl
          << description << std::endl;
  switch (severity)
  {
    case UFatalError:
      std::cerr << es_banner << message.str() << "*** Fatal Exception ***"
                << ee_banner << std::endl;
      break;
    case UFatalErrorInArguments:
      std::cerr << es_banner << message.str() << "*** Fatal Error In Argument ***"
                << ee_banner << std::endl;
      break;
    case UError:
      std::cerr << es_banner << message.str() << "*** Error ***" << level
                << ee_banner << std::endl;
      break;
    case UWarning:
      std::cerr << ws_banner << message.str() << "*** This is just a warning message ***"
                << we_banner << std::endl;
      toBeAborted = false;
      break;
    default:
      std::cout << ws_banner << message.str()
                << "*** This is just a message for your information. ***"
                << we_banner << std::endl;
      toBeAborted = false;
      break;
  }

  if (toBeAborted)
  {

    std::cerr << std::endl << "*** GException: Aborting execution ***" << std::endl;
    abort();
  }


}


/*
void UTessellatedSolid::ImportFromSTLFile(std::string filename)
{
vector <UTriangularFacet *> fFacets;

USTL::ReadFromBinaryFile(filename, fFacets);

int size = fFacets.size();
for (int i = 0; i < size; ++i)
{
UTriangularFacet *facet = fFacets[i];
AddFacet(facet);
}
SetSolidClosed(true);
}
*/

/*
int size = fFacets.size();
for (int j = 0; j < 100; ++j) //2.418 , 2.511
for (int i = 0; i < size; ++i)
{
UFacet &facet = *facetsi[j];
a += facet.GetNumberOfVertices();
}
if (a % rand() == -1) cout << a;
*/

/*
for (int j = 0; j < 100; ++j) //2.917 3.01
{
int size = fFacets.size();
for (int i = 0; i < size; ++i)
{
UFacet &facet = *fFacets[i];
a += facet.GetNumberOfVertices();
}
}
*/

/*
for (int j = 0; j < 100; ++j) // 2.589
{
std::vector<UFacet *>::const_iterator i, begin = fFacets.begin(), end = fFacets.end();
for (i = begin; i < end; ++i)
{
UFacet &facet = *(*i);
a += facet.GetNumberOfVertices();
}
}

return location;
*/
