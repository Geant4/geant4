#include "../G4Unit.hh"

// suite #1
// for convenience, define a test suite
// actually it is attached to BERT::Root, the top level test suite
testSuite(IllustrationSuite)

// Putting test cases in anonymous namespace might sometimes be useful,
// because I can freely define assistant functions and classes,
// without conflicting other developer's assistants.
namespace {

// case #1
// Now let's begin to define the first test case.
// Any test case should be attached to a suite,
// and this EqualCase is attached to IllustrationSuite defines above.
testCase(BoolExpression, IllustrationSuite)
{
    // assertTrue() performs a general assertion,
    // usually we put a bool expression that should be true between the parentheses.
    // When its value is false, a failure will be reported, indicating something wrong
    // with either the production code, or the test code, or both
    assertTrue(1 == 1u && "signed and unsigned int shall be equal");

    bool someComplicatedExpression = true;
    assertTrue(someComplicatedExpression);

    // process that reach here may be considered correct in spite of destructors
}

// case #2
// within anonymous namespace, we can freely define assistant functions.
// this function emulates some error handling code
void rejectPositive(int x)
{
    if (x > 0) throw x;
}

// the second test case attached to the same suite
testCase(NoArriveCase, IllustrationSuite)
{
    try {
        int aPositiveValue = 15;
        // With invalid parameter, we should check if an exception is thrown as expected
        rejectPositive(aPositiveValue);

        // this line is executed only if the exception is not thrown,
        // and a failure will be reported
        assertNoArrive("an invalid parameter is not rejected");
    }
    // try to use the *exact* type of exception
    // do NOT use catch (...), or some failure cannot be caught by unit--
    catch (int) {
    }
}

// case #3
// this case test for floating point type.
// NOTE: avoid using == directly,
// instead, make use of BERT::closeValue() function
testCase(FloatingPoint, IllustrationSuite)
{
    using namespace BERT;

    float expectedValue = 1.234f;
    float errorLimit = 0.001f;
    float value = 1.03f + 0.204f;

    // test if two float values were close enough to each other
    assertTrue(closeValue(expectedValue, value, errorLimit));
    assertTrue(closeValue(expectedValue, value, 0.5));

    // if some specific information is needed when an assertion failed,
    // BERT::closeValueInfo() may help
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

// case #4
// some times it is convenient to have fixtures(i.e. setUp() and tearDown() methods).
// in C++, this functionality already exists in form of constructor and destructor

// Now let's define a simple fixture
struct MyFixture {
    // variables that needed by some test case
    double x;
    double y;
    unsigned n;

    // constructor for setUp() job
    MyFixture()
    {
        n = 10;
        x = 1.3;
        y = x;
        for (unsigned i = 0; i < n; ++i) {
            y += y;
        }
    }
    // destructor for tearDown() job
    ~MyFixture()
    {
        // what ever you need
    }
};

// and this test case make use of MyFixture
testCase(FixtureCase, IllustrationSuite)
{
    // declear a local variable to invoke the constructor for setUp job
    // and (automatically) register the destructor for tearDown
    MyFixture f;

    using namespace BERT;

    double factor = 1.0;
    for (unsigned i = 0; i < f.n; ++i) {
        factor *= 2;
    }
    // check if different algorithm get same answer
    assertTrue(closeValueInfo(f.y, f.x * factor, 0.01));
}

} // namespace
