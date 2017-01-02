//
// Created by serkh on 1/1/17.
//

#ifndef HYPRG_COMMANDLINE_H
#define HYPRG_COMMANDLINE_H

#include <string>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <type_traits>
#include <exception>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>

using std::string;
using std::vector;
using std::map;
using std::set;
using std::is_same;
using std::find;
using std::cout;
using std::endl;

enum ArgType {STRING, INT, DOUBLE};

class CommandLine {
public:
    CommandLine() {}

    void add_bool_argument(const char * argname, const char* description="") {
        string arg = string(argname);
        bool_values[arg] = false;
        registered_flags.insert(arg);
        descriptions[arg] = string(description);
    }

    void add_string_argument(const char * argname, const char* description="") {
        string arg = string(argname);
        registered_arguments[arg] = STRING;
        descriptions[arg] = string(description);
    }
    void add_string_argument(const char * argname, string default_value, const char* description="") {
        add_string_argument(argname);
        set_value(argname, default_value);
    }

    void add_int_argument(const char * argname, const char* description="") {
        string arg = string(argname);
        registered_arguments[arg] = INT;
        descriptions[arg] = string(description);
    }
    void add_int_argument(const char * argname, int default_value, const char* description="") {
        add_int_argument(argname);
        set_value(argname, default_value);
    }

    void add_double_argument(const char * argname, const char* description="") {
        string arg = string(argname);
        registered_arguments[arg] = DOUBLE;
        descriptions[arg] = string(description);
    }
    void add_double_argument(const char * argname, double default_value, const char* description="") {
        add_double_argument(argname);
        set_value(argname, default_value);
    }

    void parse_arguments(int argc, char **argv) {
        vector<string> commandline;
        if (argc > 1) {
            commandline.assign(argv + 1, argv + argc);
            this->argv.assign(argv + 1, argv + argc);
        }

        // parse help
        auto help_found = find(commandline.begin(), commandline.end(), string("--help"));
        if (help_found == commandline.end()) {
            help_found = find(commandline.begin(), commandline.end(), string("-h"));
        }

        if (help_found != commandline.end()) {
            print_help();
            std::exit(0);
        }
        // if argument has one-character name, parse with leading "-", otherwise parse "--"
        // parse flags
        for (const string& bool_flag : registered_flags) {
            string val = val_to_search(bool_flag);
            auto found = find(commandline.begin(), commandline.end(), val);
            if (found != commandline.end()) {
                bool_values[bool_flag] = true;
                commandline.erase(found);
            }
        }
        if (commandline.size() % 2 != 0) {
            throw std::logic_error("incorrect commandline");
        }
        for (auto arg : registered_arguments) {
            const char* argname = arg.first.c_str();
            string val = val_to_search(arg.first);
            ArgType type = arg.second;
            auto iter_found_argname = find(commandline.begin(), commandline.end(), val);
            if (iter_found_argname != commandline.end()) {
                auto iter_found_argval = iter_found_argname + 1;
                if (iter_found_argval == commandline.end() or (*iter_found_argval).find("-") == 0) {
                    throw std::logic_error("incorrect commandline");
                }
                if (type == STRING) {
                    set_value(argname, *iter_found_argval);
                } else if (type == INT) {
                    set_value(argname, std::stoi(*iter_found_argval));
                } else if (type == DOUBLE) {
                    set_value(argname, std::stod(*iter_found_argval));
                }
                commandline.erase(iter_found_argval);
                commandline.erase(iter_found_argname);
            }
        }
        if (commandline.size() > 0) {
            for (const string& c : commandline) {
                cout << c << " ";
            }
            cout << endl;
            throw std::logic_error("unrecognized commandline part") ;
        }
    }

    bool has_bool_argument(const char* argname) {
        string arg = string(argname);
        if(registered_flags.find(arg) == registered_flags.end()) {
            throw std::domain_error("unknown argument checked");
        }
        return bool_values.find(arg) != bool_values.end();
    }
    bool has_string_argument(const char* argname) {
        string arg = string(argname);
        if(registered_arguments.find(arg) == registered_arguments.end()) {
            throw std::domain_error("unknown argument checked");
        }
        return string_values.find(arg) != string_values.end();
    }
    bool has_int_argument(const char* argname) {
        string arg = string(argname);
        if(registered_arguments.find(arg) == registered_arguments.end()) {
            throw std::domain_error("unknown argument checked");
        }
        return int_values.find(arg) != int_values.end();
    }
    bool has_double_argument(const char* argname) {
        string arg = string(argname);
        if(registered_arguments.find(arg) == registered_arguments.end()) {
            throw std::domain_error("unknown argument checked");
        }
        return double_values.find(arg) != double_values.end();
    }

    bool get_bool_argument(const char * argname) {
        string arg = string(argname);
        return bool_values[arg];
    }
    int get_int_argument(const char* argname) {
        string arg = string(argname);
        return int_values[arg];
    }
    double get_double_argument(const char* argname) {
        string arg = string(argname);
        return double_values[arg];
    }
    string get_string_argument(const char* argname) {
        string arg = string(argname);
        return string_values[arg];
    }

    void print_help() {
        cout << "Optional arguments:" << endl;
        for (auto d : descriptions) {
            cout << "\t" << val_to_search(d.first) << "\t\t\t" << d.second << endl;
        }
    }
private:
    map<string, ArgType> registered_arguments;
    set<string> registered_flags;

    map<string, bool> bool_values;
    map<string, int> int_values;
    map<string, double> double_values;
    map<string, string> string_values;

    map<string, string> descriptions;
    vector<string> argv;

    void set_value(const char* argname, string value) {
        string_values[string(argname)] = value;
    }

    void set_value(const char* argname, int value) {
        int_values[string(argname)] = int(value);
    }

    void set_value(const char* argname, double value) {
        double_values[string(argname)] = double(value);
    }

    string val_to_search(string argname) {
        if (argname.size() == 1) {
            return "-" + argname;
        } else {
            return "--" + argname;
        }
    }
};

#endif //HYPRG_COMMANDLINE_H
