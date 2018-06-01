#include <fstream>

#include "gtest/gtest.h"
#include "serialization.hpp"
#include "hashers.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";

TEST(Serialization, Number) {
    uint64_t test_num = 10;
    std::ofstream out(test_dump_basename + "_testser");
    serialization::serializeNumber(out, test_num);
    out.close();

    std::ifstream in(test_dump_basename + "_testser");
    EXPECT_EQ(test_num, serialization::loadNumber(in));
}

TEST(Serialization, String) {
    std::string test_string("test_string");
    std::ofstream out(test_dump_basename + "_testser");
    serialization::serializeString(out, test_string);
    out.close();

    std::ifstream in(test_dump_basename + "_testser");
    EXPECT_EQ(test_string, serialization::loadString(in));
}

TEST(Serialization, Mixed) {
    uint64_t test_num = -2llu;
    std::string test_string("test_string");

    std::ofstream out(test_dump_basename + "_testser");
    serialization::serializeNumber(out, test_num);
    serialization::serializeNumber(out, test_num);
    serialization::serializeString(out, test_string);
    serialization::serializeString(out, test_string);
    serialization::serializeNumber(out, test_num);
    serialization::serializeString(out, "");
    serialization::serializeNumber(out, test_num);
    serialization::serializeString(out, test_string);
    serialization::serializeString(out, test_string);
    serialization::serializeNumber(out, test_num);
    serialization::serializeString(out, test_string);
    out.close();

    std::ifstream in(test_dump_basename + "_testser");
    EXPECT_EQ(test_num, serialization::loadNumber(in));
    EXPECT_EQ(test_num, serialization::loadNumber(in));
    EXPECT_EQ(test_string, serialization::loadString(in));
    EXPECT_EQ(test_string, serialization::loadString(in));
    EXPECT_EQ(test_num, serialization::loadNumber(in));
    EXPECT_EQ(std::string(""), serialization::loadString(in));
    EXPECT_EQ(test_num, serialization::loadNumber(in));
    EXPECT_EQ(test_string, serialization::loadString(in));
    EXPECT_EQ(test_string, serialization::loadString(in));
    EXPECT_EQ(test_num, serialization::loadNumber(in));
    EXPECT_EQ(test_string, serialization::loadString(in));
}

TEST(Serialization, Zero) {
    for (size_t num = 0; num < 1000000; num += 1000) {
        std::ofstream out(test_dump_basename + "_serial");
        serialization::serializeNumber(out, num);
        out.close();

        std::ifstream in(test_dump_basename + "_serial");
        ASSERT_EQ(num, serialization::loadNumber(in));
    }
    std::ofstream out(test_dump_basename + "_serial");
    serialization::serializeNumber(out, -1llu);
    out.close();

    std::ifstream in(test_dump_basename + "_serial");
    ASSERT_EQ(-1llu, serialization::loadNumber(in));
}

