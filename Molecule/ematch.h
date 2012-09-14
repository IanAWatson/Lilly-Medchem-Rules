/**************************************************************************

    Copyright (C) 2012  Eli Lilly and Company

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
#ifndef E_MATCH_H
#define E_MATCH_H

#include "iwaray.h"
#include "iwcrex.h"

#include "element.h"

/*
  Atom matching capability. 
*/

class Element_Matcher
{
  private:
    const Element * _e;

    int _isotope;

    int _match_organic_only;
    int _match_non_organic_only;
    int _match_non_periodic_only;

    IW_Regular_Expression _symbol_rx;

//  private functions

    void _default_values ();

  public:
    Element_Matcher ();
    Element_Matcher (atomic_number_t);
    Element_Matcher (const Element *);
    Element_Matcher (const char *);
    Element_Matcher (const IWString &);

    int ok () const;
    int debug_print (ostream &) const;

    void set_element (const Element *);
    int  construct_from_string (const char *, int);
    int  construct_from_string (const char *);
    int  construct_from_string (const const_IWSubstring &);
    int  construct_from_string (const IWString &);

    int operator_less_less (ostream & os) const;

    const Element * element () const { return _e;}

    int isotope () const { return _isotope;}

    int matches (const Element *, int = 0);    // = 0 parameter is isotope
};

extern ostream &
operator<< (ostream &, const Element_Matcher &);

class Command_Line;

class Set_of_Element_Matches : public resizable_array_p<Element_Matcher>
{
  private:
  public:

    int construct_from_command_line (Command_Line &, int, char);

    int matches (const Element *, int = 0);    // = 0 parameter is isotope
};

extern void display_element_matcher_syntax (ostream & os);

#endif
