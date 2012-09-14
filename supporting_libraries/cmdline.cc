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
#include "iwconfig.h"
#include <stdlib.h>
#include <stdio.h>

#include <iostream>

#ifdef _WIN32
#include "getopt.h"

/* Global Exportable */
extern int optind;
extern char *optarg;

#else
#include <unistd.h>
#include <getopt.h>
#endif

#include "cmdline.h"
#include "iwstring.h"

#define CL_MAGIC 97531
#define OV_MAGIC 97532

#define UNKNOWN_VALIDITY -1

Option_and_Value::Option_and_Value (int o, const char * val) : _o (o), _value (val)
{
  _valid_as_int = UNKNOWN_VALIDITY;
  _valid_as_double = UNKNOWN_VALIDITY;

  _magic = OV_MAGIC;

  return;
}

int
Option_and_Value::ok () const
{
  if (OV_MAGIC != _magic)
    return 0;

  return 1;
}

Option_and_Value::~Option_and_Value ()
{
  assert (ok ());

  _magic = -5;
}

int
Option_and_Value::value (int & rc)
{
  if (NULL == _value)
    return 0;

  if (UNKNOWN_VALIDITY == _valid_as_int)
    _valid_as_int = is_int (_value, &_int_val);

  if (_valid_as_int)
  {
    _double_val = double (_int_val);
    _valid_as_double = 1;
    rc = _int_val;
    return 1;
  }

  return 0;
}

int
Option_and_Value::value (unsigned int & rc)
{
  int tmp;
  if (! value (tmp))
    return 0;

  if (tmp < 0) 
    return 0;

  rc = static_cast<unsigned int> (tmp);

  return 1;
}

int
Option_and_Value::value (long int & rc)
{
  if (NULL == _value)
    return 0;

  if (UNKNOWN_VALIDITY == _valid_as_int)
    _valid_as_int = is_int (_value, &_int_val);

  if (_valid_as_int)
  {
    _valid_as_int = 1;
    _valid_as_double = 1;
    _double_val = double (_int_val);
    rc = (long int) _int_val;
    return 1;
  }

  return 0;
}

/*
  Pretty dangerous stuff here, no checking for overflow
*/

int
Option_and_Value::value (long long & rc)
{
  if (NULL == _value)
    return 0;

  if (UNKNOWN_VALIDITY == _valid_as_int)
    _valid_as_int = is_int (_value, &_int_val);

  if (_valid_as_int)
  {
    _valid_as_int = 1;
    _valid_as_double = 1;
    _double_val = double (_int_val);
    rc = (long long) _int_val;
    return 1;
  }

  return 0;
}

int
Option_and_Value::value (unsigned long & rc)
{
  if (NULL == _value)
    return 0;

  if (UNKNOWN_VALIDITY == _valid_as_int)
    _valid_as_int = is_int (_value, &_int_val);

  if (_valid_as_int)
  {
    _valid_as_int = 1;
    _valid_as_double = 1;
    _double_val = double (_int_val);
    rc = (unsigned long) _int_val;
    return 1;
  }

  return 0;
}

int
Option_and_Value::value (double & rc)
{
  if (NULL == _value)
    return 0;

  if (UNKNOWN_VALIDITY == _valid_as_double)
    _valid_as_double = is_double (_value, &_double_val);
  
  if (_valid_as_double)
  {
    rc = _double_val;
    return 1;
  }
  else
    return 0;
}

int
Option_and_Value::value (float & rc)
{
  if (NULL == _value)
    return 0;

  if (UNKNOWN_VALIDITY == _valid_as_double)
    _valid_as_double = is_double (_value, &_double_val);
  
  if (_valid_as_double)
  {
    rc = float (_double_val);    // no checks for over/under flow!!!
    return 1;
  }
  else
    return 0;
}

int
Option_and_Value::value (IWString & result)
{
  if (NULL == _value)
    return 0;

  result = _value;

  return 1;
}

int
Option_and_Value::value (const_IWSubstring & result)
{
  if (NULL == _value)
    return 0;

  result = _value;

  return 1;
}

int
Option_and_Value::value (char * buffer)
{
  if (NULL == _value)
    return 0;

  strcpy (buffer, _value);

  return 1;
}

ostream &
operator << (ostream & os, const Option_and_Value & ov)
{
  return os << "Option '" << ov.option () << "', value '" << ov.value () << "'";
}

Command_Line::Command_Line (int argc, char ** argv, const char * options)
{
  _magic = CL_MAGIC;

//argptr = NULL;
  optarg = NULL;

#ifdef _WIN32
  optind = 0;
#else
  optind = 1;   // reinitialise in case of multiple invocations
  opterr = 0;     // suppress error messages
#endif

  _options.resize (argc);     // yes, this is too many

  _some_options_start_with_dash = 0;
  _unrecognised_options_encountered = 0;

  int o;


  while ((o = getopt (argc, argv, options)) != EOF)
  {
    if ('?' == o)
    {
      cerr << "Command_Line: unrecognised option '" << argv[optind - 1] << "'\n";
      _unrecognised_options_encountered++;
    }
    else
    {
      Option_and_Value * tmp = new Option_and_Value (o, optarg);
      _options.add (tmp);
    }
  }

  resize (argc - optind);

  for (int i = optind; i < argc; i++)
  {
    if ('-' == *(argv[i]))
      _some_options_start_with_dash++;
    add (argv[i]);
  }

  return;
}

Command_Line::~Command_Line ()
{
  assert (ok ());

  _magic = -9;
}

int
Command_Line::ok () const
{
  if (_magic != CL_MAGIC)
    return 0;

  for (int i = 0; i < _options.number_elements (); i++)
    if (! _options[i]->ok ())
      return 0;

  return resizable_array<const char *>::ok ();
}

int
Command_Line::debug_print (ostream & os) const
{
  assert (os.good ());

  os << "Command line object contains " << _options.number_elements () << " options and " <<
        _number_elements << " values\n";

  for (int i = 0; i < _options.number_elements (); i++)
    os << *(_options[i]) << endl;

  for (int i = 0; i < _number_elements; i++)
    os << "Value '" << _things[i] << "'\n";

  return 1;
}

int
Command_Line::option_present (const char c) const
{
  for (int i = 0; i < _options.number_elements (); i++)
  {
    const Option_and_Value * oo = _options[i];
    if (c == oo->option ())
      return i + 1;
  }

  return 0;
}

const char *
Command_Line::option_value (const char c, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements (); i++)
  {
    const Option_and_Value * oo = _options[i];
    if (c == oo->option () && occurrence == nfound++)
      return oo->value ();
  }

  return 0;
}

int
Command_Line::all_values (const char c, resizable_array<const char *> & values) const
{
  int rc = 0;

  for (int i = 0; i < _options.number_elements (); i++)
  {
    const Option_and_Value * oo = _options[i];
    if (c == oo->option ())
    {
      values.add (oo->value ());
      rc++;
    }
  }
  return rc;
}

/*
  When returning all instances of an option, we can optionally
  further tokenise these instances. 

  prog -r a -r 'c d'

  would return 3 elements (if split_tokens is specified)
*/

int
Command_Line::all_values (const char c, resizable_array_p<IWString> & values,
                          int split_tokens) const
{
  int rc = 0;

  for (int i = 0; i < _options.number_elements (); i++)
  {
    const Option_and_Value * oo = _options[i];
    if (c != oo->option ())
      continue;


    if (0 == split_tokens || 0 == ::ccount (oo->value (), ' '))    // no need to split, either just one token, or not requested
    {
      rc++;
      IWString * tmp = new IWString (oo->value ());
      values.add (tmp);
      continue;
    }

    const_IWSubstring v (oo->value ());
    int j = 0;
    const_IWSubstring token;
    while (v.nextword (token, j))
    {
      IWString * tmp = new IWString (token);
      values.add (tmp);
      rc++;
    }
  }

  return rc;
}

int
Command_Line::option_count (const char c) const
{
  int rc = 0;

  for (int i = 0; i < _options.number_elements (); i++)
  {
    const Option_and_Value * oo = _options[i];
    if (c == oo->option ())
      rc++;
  }

  return rc;
}

int
Command_Line::value (const char c, int & result, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements (); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option ())
    {
      if (nfound == occurrence)
        return oo->value (result);
      nfound++;
    }
  }

  return 0;
}

int
Command_Line::value (const char c, unsigned int & result, int occurrence) const
{
  int tmp;
  if (! value (c, tmp, occurrence))
    return 0;

  if (tmp < 0)
    return 0;

  result = static_cast<unsigned int> (tmp);

  return 1;
}

int
Command_Line::value (const char c, long int & result, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements (); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option ())
    {
      if (nfound == occurrence)
        return oo->value (result);
      nfound++;
    }
  }

  return 0;
}

int
Command_Line::value (const char c, long long & result, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements (); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option ())
    {
      if (nfound == occurrence)
        return oo->value (result);
      nfound++;
    }
  }

  return 0;
}

int
Command_Line::value (const char c, unsigned long & result, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements (); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option ())
    {
      if (nfound == occurrence)
        return oo->value (result);
      nfound++;
    }
  }

  return 0;
}

/*int
Command_Line::double_val (const char c, double & value, int occurrence)
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements (); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option ())
    {
      if (nfound == occurrence)
        return oo->double_val (value);
      nfound++;
    }
  }

  return 0;
}*/

int
Command_Line::value (const char c, double & result, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements (); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option ())
    {
      if (nfound == occurrence)
        return oo->value (result);
      nfound++;
    }
  }

  return 0;
}

int
Command_Line::value (const char c, float & result, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements (); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option ())
    {
      if (nfound == occurrence)
        return oo->value (result);
      nfound++;
    }
  }

  return 0;
}

int
Command_Line::value (const char c, IWString & result, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements (); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option ())
    {
      if (nfound == occurrence)
        return oo->value (result);
      nfound++;
    }
  }

  return 0;
}

int
Command_Line::value (const char c, const_IWSubstring & result, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements (); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option ())
    {
      if (nfound == occurrence)
        return oo->value (result);
      nfound++;
    }
  }

  return 0;
}

int
Command_Line::value (const char c, char * result, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements (); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option ())
    {
      if (nfound == occurrence)
        return oo->value (result);
      nfound++;
    }
  }

  return 0;
}

const const_IWSubstring
Command_Line::string_value (const char c, int occurrence) const
{
  int nfound = 0;
  for (int i = 0; i < _options.number_elements (); i++)
  {
    Option_and_Value * oo = _options[i];
    if (c == oo->option ())
    {
      if (nfound == occurrence)
      {
        const char * v = oo->value ();
        if (NULL == v)
          return "";
        else
          return oo->value ();
      }

      nfound++;
    }
  }

  return "";
}

ostream &
operator << (const Command_Line & cl, ostream & os)
{
  return os;
}
