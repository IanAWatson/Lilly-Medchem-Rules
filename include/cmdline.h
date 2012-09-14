/**************************************************************************

    Copyright (C) 2011  Eli Lilly and Company

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

**************************************************************************/
#ifndef IW_CMDLINE_H
#define IW_CMDLINE_H

#include <iostream>

#include "iwstring.h"

typedef int clov_magic_number_t;

class Option_and_Value
{
  friend
    std::ostream &
    operator << (std::ostream &, const Option_and_Value &);

  private:
    int _o;
    const char * _value;
    int _valid_as_int;
    int _int_val;
    double _double_val;
    int _valid_as_double;
    clov_magic_number_t _magic;

  public:
    Option_and_Value (int, const char * = NULL);
    ~Option_and_Value ();

    int ok () const;

    char   option () const { return _o;}
    const char * value  () const { return _value;}

//  int int_val (int &);
//  int double_val (double &);

    int value (int &);
    int value (unsigned int &);
    int value (long int &);
    int value (long long &);
    int value (unsigned long &);
    int value (float &);
    int value (double &);
    int value (IWString &);
    int value (const_IWSubstring &);
    int value (char *);
};


class Command_Line : public resizable_array<const char *>
{
  friend
    std::ostream &
      operator << (const Command_Line &, std::ostream &);

  private:
    resizable_array_p<Option_and_Value> _options;
    int _some_options_start_with_dash;    // perhaps indicative of an error
    int _unrecognised_options_encountered;
    clov_magic_number_t _magic;

  public:
    Command_Line (int, char **, const char *);
    ~Command_Line ();

    int ok () const;
    int debug_print (ostream &) const;

    int some_options_start_with_dash () const { return _some_options_start_with_dash;}
    int unrecognised_options_encountered () const { return _unrecognised_options_encountered;}

    Option_and_Value * ov (const char, int = 0);

    int option_present (const char) const;
    int option_count (const char) const;

    const char * option_value (const char, int = 0) const;

    int value (const char, int &, int = 0) const;
    int value (const char, unsigned int &, int = 0) const;
    int value (const char, long int &, int = 0) const;
    int value (const char, long long &, int = 0) const;
    int value (const char, unsigned long &, int = 0) const;
    int value (const char, double &, int = 0) const;
    int value (const char, float &, int = 0) const;
    int value (const char, IWString &, int = 0) const;
    int value (const char, const_IWSubstring &, int = 0) const;
    int value (const char, char *, int = 0) const;

//#ifdef IW_STD_STRING_DEFINED
    int value (const char, std::string &, int = 0) const;
    std::string std_string_value(const char, int = 0) const;
//#endif

    const const_IWSubstring string_value (const char, int = 0) const;

    int all_values (const char, resizable_array<const char *> &) const;
    int all_values (const char, resizable_array_p<IWString> &, int = 0) const;

//  int position (char c, int f = 0) const { return option_present (c, f) - 1;}
};

#endif
