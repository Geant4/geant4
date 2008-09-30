#include "../G4Unit.hh"
#include "G4HeaderCase.hh"

testSuite(StaticLink1)

namespace {
    unsigned i = 0;
    unsigned x()
    {
        return 2;
    }
} // namespace

static unsigned y()
{
    return 2;
}

namespace {

testCase(IncStatic2, StaticLink1)
{
    using namespace BERT;
    ++i;
    // varibles in anonymous namespace cannot be seen by other files
    assertTrue(equalValueInfo(1u, i));
    assertTrue(2 == x() && "functions in anonymous namespace cannot be redefined from other files");
    assertTrue(2 == y() && "static function has no external linkage");
}

} // namespace
