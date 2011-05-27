#include <vector>
#include <numeric>
#include <algorithm>
#include "../G4Unit.hh"

testSuite(TemplateSuite)

template <typename T>
void testAlgorithms()
{
    using namespace std;
    using namespace BERT;
    vector<T> ve(100, 1);
    partial_sum(ve.begin(), ve.end(), ve.begin());

    assertTrue(ve.size() > 0);
    assertTrue(1 == ve[0]);
    for (unsigned i = 1; i < ve.size(); ++i) {
        assertTrue(equalValueInfo(ve[i - 1] + 1, ve[i]));
    }
}

namespace {
testCase(IntCase, TemplateSuite)
{
    testAlgorithms<int>();
}

testCase(UnsignedCase, TemplateSuite)
{
    testAlgorithms<unsigned>();
}

} // namespace
