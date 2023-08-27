#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "printing.h"
#include <sstream>
#include <string>
#include <vector>

TEST_CASE("printing") {
    SUBCASE("print_vector") {
        std::vector<double> x{1.0, 2.0, 3.0, 4.0};

        // redirect std::cout to buffer
        std::stringstream buffer;
        std::streambuf *prev_cout_buf = std::cout.rdbuf(buffer.rdbuf());

        // BEGIN: code being tested
        print_vector("x", x);
        // END: code being tested

        // restore original buffer before exiting
        std::cout.rdbuf(prev_cout_buf);

        // use the string value of buffer to compare against expected output
        std::string expected = "x = 1, 2, 3, 4, \n";
        std::string text = buffer.str();
        int result = text.compare(expected);

        CHECK(result == 0);
    }

    SUBCASE("print_scientific") {
        // redirect std::cout to buffer
        std::stringstream buffer;
        std::streambuf *prev_cout_buf = std::cout.rdbuf(buffer.rdbuf());

        // BEGIN: code being tested
        print_scientific(123.456, 11, 5);
        // END: code being tested

        // restore original buffer before exiting
        std::cout.rdbuf(prev_cout_buf);

        // use the string value of buffer to compare against expected output
        std::string expected = "1.23456e+02";
        std::string text = buffer.str();
        int result = text.compare(expected);

        CHECK(result == 0);
    }
}
