#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <cstddef>
#include <cstring>
#include <vector>
#include <limits>
#include <map>
#include <cstdint>


class Util {
public:

    static std::string base_name(std::string const & path, std::string const & delims)
    {
        return path.substr(path.find_last_of(delims) + 1);
    }

    static inline size_t skipWhitespace(const char* data) {
        size_t counter = 0;
        while ((data[counter] == ' ' || data[counter] == '\t') == true) {
            counter++;
        }
        return counter;
    }

    static inline size_t skipNoneWhitespace(const char * data){
        //A value different from zero (i.e., true) if indeed c is a white-space character. Zero (i.e., false) otherwise.
        size_t counter = 0;
        while((data[counter] == ' '  || data[counter] == '\t' || data[counter] == '\n' || data[counter] == '\0') == false) {
            counter++;
        }
        return counter;
    }

    static inline size_t getWordsOfLine(const char* data, const char** words, size_t maxElement ){
        size_t elementCounter = 0;
        while(*data !=  '\n' && *data != '\0'){
            data += skipWhitespace(data);
            words[elementCounter] = data;
            elementCounter++;
            if(elementCounter >= maxElement)
                return elementCounter;
            data += skipNoneWhitespace(data);
        }
        if(elementCounter < maxElement)
            words[elementCounter] = data;

        return elementCounter;
    }
};
#endif
