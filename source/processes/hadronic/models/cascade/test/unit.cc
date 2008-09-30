// usage demonstrated in unit/GNUmakefile
#include "unit/G4Unit.hh"

testSuite(IllustrationSuite)
namespace {

testCase(FloatingPoint, IllustrationSuite)
{
    using namespace BERT;

    float expectedValue = 1.234f;
    float errorLimit = 0.001f;
    float value = 1.03f + 0.204f;

    // test if two float values were close enough to each other
    assertTrue(closeValue(expectedValue, value, errorLimit));
    assertTrue(closeValue(expectedValue, value, 0.5));
      assertTrue(closeValue(1.0, 2.0, 0.1)); // AH::: fails    

   
    // specific information when an assertion failed: BERT::closeValueInfo()
    // NOTE: closeValueInfo() uses operator << to convert parameters to string description
    // be sure a proper operator << is provided for each parameter
    // it is already defined for float, double and long double
    assertTrue(closeValueInfo(expectedValue, value, errorLimit));
    assertTrue(closeValueInfo(1.66666666, 5.0 / 3.0, 0.001));

    // there is also a function for equality test: BERT::equalValueInfo()
    // NOTE: proper operator << for every parameter is also needed
    const int expectedInt = 8;
    assertTrue(equalValueInfo(expectedInt, 0x01 << 3));
}

}// namespace
