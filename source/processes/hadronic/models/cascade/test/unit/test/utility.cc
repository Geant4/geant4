#include <string>
#include <sstream>
#include <stdexcept>
#include "../G4Unit.hh"

testSuite(UtilitySuite)

namespace {

// test if AssertionInfo transfer info correctly
testCase(AssertionInfoCase, UtilitySuite)
{
    {
        BERT::AssertionInfo info(true);
        assertTrue(true  == info.successful());
        assertTrue(false == info.hasMessage());
    }

    {
        BERT::AssertionInfo info(false);
        assertTrue(false == info.successful());
        assertTrue(false == info.hasMessage());
        assertTrue("" == info.message());
    }

    {
        const std::string message = "arBitRarY sTriNg";
        BERT::AssertionInfo info(false, message);
        assertTrue(false == info.successful());
        assertTrue(true  == info.hasMessage());
        assertTrue(message == info.message());
    }
}

template <typename T>
inline T strTo(const std::string& str)
{
    std::istringstream ii(str);
    T result;
    if (ii >> result) return result;
    else throw std::invalid_argument("cannot convert (" + str + ")");
}

// find sub string enclosed by < and >, starting at "end"
void getNumber(const std::string& message, unsigned& begin, unsigned& end)
{
    using std::string;
    begin = static_cast<unsigned>(message.find("<", end));
    assertTrue(string::npos != begin);
    ++begin;
    assertTrue(message.length() > begin);
    end = static_cast<unsigned>(message.find(">", begin));
    assertTrue(string::npos != end);
}

// test if equalValueInfo give correct hint about failure
testCase(EqualValueCase, UtilitySuite)
{
    {
        BERT::AssertionInfo info(BERT::equalValue(1 + 1, 2));
        assertTrue(info.successful());
    }

    {
        const unsigned expected = 1;
        const unsigned actual = 2;
        BERT::AssertionInfo info(BERT::equalValueInfo(expected, actual));
        assertTrue(!info.successful());
        assertTrue( info.hasMessage());

        using std::string;
        string message = info.message();
        assertTrue(message.length() > 0);

        unsigned begin = 0;
        unsigned end = 0;
        getNumber(message, begin, end);
        assertTrue(expected == strTo<unsigned>(message.substr(begin, end - begin)));

        getNumber(message, begin, end);
        assertTrue(actual == strTo<unsigned>(message.substr(begin, end - begin)));
    }
}

// test if closeValueInfo give correct hint about failure
testCase(CloseValueCase, UtilitySuite)
{
    {
        BERT::AssertionInfo info(BERT::closeValue(2.0, 2.0, 0.01));
        assertTrue(info.successful());
    }

    {
        const float expected = 0.0f + 3.5f - 1.9f;
        const double actual = 2.0;
        const double prec = 0.01;
        BERT::AssertionInfo info(BERT::closeValueInfo(expected, actual, prec));
        assertTrue(!info.successful());
        assertTrue( info.hasMessage());

        using std::string;
        string message = info.message();
        assertTrue(message.length() > 0);

        unsigned begin = 0;
        unsigned end = 0;
        getNumber(message, begin, end);
        const float expectedFromMessage = strTo<float>(message.substr(begin, end - begin));
        assertTrue(BERT::closeValueInfo(expected, expectedFromMessage, 0.001));

        getNumber(message, begin, end);
        const double actualFromMessage = strTo<double>(message.substr(begin, end - begin));
        assertTrue(BERT::closeValueInfo(actual, actualFromMessage, 0.001));

        getNumber(message, begin, end);
        const double precFromMessage = strTo<double>(message.substr(begin, end - begin));
        assertTrue(BERT::closeValueInfo(prec, precFromMessage, 0.001));
    }
}


} // namespace
