//
// Created by serkh on 1/1/17.
//

#ifndef HYPRG_COMMANDLINE_H
#define HYPRG_COMMANDLINE_H

#include <string>
#include <vector>

using std::string;
using std::vector;
class CommandLine {
public:
    CommandLine(int, char**);

    bool parse_flag(const char*);
    int parse_int(const char*);
    double parse_double(const char*);
    string parse_string(const char*);

    string operator[](unsigned);
    unsigned size();
private:
    vector<string> argv;
};


#endif //HYPRG_COMMANDLINE_H
