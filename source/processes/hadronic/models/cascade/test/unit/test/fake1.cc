#include "../G4Unit.hh"
#include "G4HeaderCase.hh"

testSuite(StaticLink1)

namespace {
    unsigned i = 0;
    unsigned x()
    {
        return 1;
    }
} // namespace

static unsigned y()
{
    return 1;
}

namespace {

testCase(IncStatic1, StaticLink1)
{
    using namespace BERT;
    ++i;
    // varibles in anonymous namespace cannot be seen by other files
    assertTrue(equalValueInfo(1u, i));
    assertTrue(1 == x() && "functions in anonymous namespace cannot be redefined from other files");
    assertTrue(1 == y() && "static function has no external linkage");
}

} // namespace
