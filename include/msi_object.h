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
#ifndef IW_MSI_OBJECT_H
#define IW_MSI_OBJECT_H 1

#include <iostream>

using namespace std;

#include "iwstring.h"
#include "iwstring_data_source.h"

#define MSI_ATTRIBUTE_TYPE_INT     3
#define MSI_ATTRIBUTE_TYPE_FLOAT   4
#define MSI_ATTRIBUTE_TYPE_DOUBLE  5
#define MSI_ATTRIBUTE_TYPE_STRING  6
#define MSI_ATTRIBUTE_TYPE_OBJECT  8

class msi_attribute
{
  friend
    ostream &
      operator << (ostream &, const msi_attribute &);
  private:
    IWString _name;
    int    _type;    // bad design, yes....

    IWString _string;

    resizable_array_p<IWString> _string_values;
    resizable_array<int>        _int_values;
    resizable_array<double>     _double_values;

    int    _valid;

    void _default_values ();

  public:
    msi_attribute (); //char *);
    ~msi_attribute ();

    int ok () const;

    int build (const IWString &);

    int valid_as_int    () const { return _int_values.number_elements ();}
    int valid_as_double () const { return _double_values.number_elements ();}

    int number_double_values   () const { return _double_values.number_elements ();}
    int number_int_values      () const { return _int_values.number_elements ();}
    int number_string_values   () const;

    const IWString & name () const { return _name;}
    const IWString & stringval () const { return _string;}

    int value (int &) const;
    int value (unsigned int &) const;
    int value (float &) const;
    int value (double &) const;
    int value (IWString &) const;
    int value (const_IWSubstring &) const;

    int  fetch_int_multiple_values (int, int *) const;
    int  int_multi_value (int i) const { return _int_values[i];}

    int  fetch_double_multiple_values (int, double *) const;
    int  fetch_double_multiple_values (resizable_array<double> & result) const
           { return result.copy (_double_values);}
    double double_multi_value (int i) const { return _double_values[i];}

    const IWString * string_multi_value (int i) const;

//  First arg is the value to be returned, the 2nd arg is just an index

    int next_value (int &, int &) const;
    int next_value (double &, int &) const;
    int next_value (IWString &, int &) const;
};

class msi_object : public resizable_array_p<msi_object>
{
  friend
    ostream & operator << (ostream &, const msi_object &);

  private:
    resizable_array_p<msi_attribute> _attributes;
    int      _object_id;
    IWString _name;
    int      _valid;

//  For efficiency, we can ignore various strings when reading attributes
//  ideally, these would be regular expressions, but I don't trust the rx
//  class yet.

    resizable_array_p<IWString> _attributes_to_ignore;

//  Private functions

    int _matches_a_rejection (const IWString &) const;
    
  public:
    msi_object ();
    ~msi_object ();

    int ok () const;
    int print (ostream &, int = 0) const;

    int object_id () const { return _object_id;}

    int active () const { return object_count () + attribute_count ();}

    int object_count () const;
    int highest_object_id () const;

    const IWString & name () const { return _name;}

    int ignore (const char *);

    int read (iwstring_data_source &);

    void names_to_lowercase ();

    const msi_attribute * attribute (const char *, int = 0) const;
    const msi_attribute * attribute (const IWString &, int = 0) const;

    const msi_attribute * attribute (int) const;

    const msi_object * component (const const_IWSubstring &, int = 0) const;

    int  add_attribute (msi_attribute * extra) {return _attributes.add (extra);}

    int  attribute_count () const { return _attributes.number_elements ();};
    int  number_attributes () const { return _attributes.number_elements ();}
    int  attribute_count (const char *) const;

    int  object_count    (int) const;
    int  object_count    (const char *) const;

    int  string_value_for_attribute (const char *, IWString &) const;
};

extern void set_convert_tags_to_lowercase (int);

#endif
