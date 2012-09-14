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
#ifndef QRY_AND_DEMERIT_H
#define QRY_AND_DEMERIT_H

#include "iwaray.h"

#include "substructure.h"

class IWString;
class Demerit;
class Molecule_to_Match;

class Query_and_Demerit_Value: public Substructure_Query
{
  private:
    IWString _description;
    int _reject;
    int _demerit;
    int _demerit_each_occurrence;

    int _molecules_examined;
    int _molecules_demerited;

//  private functions

    void _default_values ();

  public:
    Query_and_Demerit_Value ();
    Query_and_Demerit_Value (const const_IWSubstring &);

    int debug_print (ostream &) const;
    int ok () const;

    int set_rejection (int);
    int determine_action (const IWString &);

    int evaluate (Molecule_to_Match &, Demerit &);
};

#endif
