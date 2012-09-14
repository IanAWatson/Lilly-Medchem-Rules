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
#ifndef ASMARTS_COMPONENT_H
#define ASMARTS_COMPONENT_H

#include "iwstring.h"

class Atomic_Smarts_Component : public const_IWSubstring
{
  private:
    int _unary_operator;
    int _op;              // the operator following the token
    Atomic_Smarts_Component * _next;

//  private functions

    int _parse (const_IWSubstring &);

  public:
    Atomic_Smarts_Component ();
    ~Atomic_Smarts_Component ();

    int ok () const;
    int debug_print (ostream &) const;

    Atomic_Smarts_Component * next () const { return _next;}

    int op () const { return _op;}
    int unary_operator () const { return _unary_operator;}

    int parse (const_IWSubstring);
};

ostream & operator << (ostream &, const Atomic_Smarts_Component &);


#endif
