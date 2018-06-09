#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <string>


namespace utils {

    /**
     *  This function checks whether two given strings are identical.
     */
    template <class String>
    bool seq_equal(const String &s1, const String &s2, size_t start = 0) {
        if (s1.size() != s2.size())
            return false;

        for (size_t i = start; i < s1.size(); ++i) {
            if (s1.at(i) != s2.at(i))
                return false;
        }
        return true;
    }

    /**
     *  This function checks whether string s1 is co-lexicographically
     *  greater than s2.
     */
    template <class String>
    bool colexicographically_greater(const String &s1, const String &s2) {
        size_t ss1 = s1.size();
        size_t ss2 = s2.size();
        for (size_t i = 1; i <= std::min(ss1, ss2); ++i) {
            if (s1.at(ss1 - i) != s2.at(ss2 - i))
                return (s1.at(ss1 - i) > s2.at(ss2 - i));
        }
        return ss1 > ss2;
    }

    std::string get_filetype(const std::string &fname);

} // namespace utils

#endif // __UTILS_HPP__
