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
#ifndef IW_DEMERIT_H
#define IW_DEMERIT_H

#include <iostream>
#include "iwstring.h"

#define DEFAULT_REJECTION_THRESHOLD 100

/*
  The variable _rejected keep track of whether or not a molecule
  has been specifically rejected by any one rule
*/

class Demerit
{
  private:
    int _score;
    int _number_different_demerits_applied;
    IWString _types;

//  private functions

    void _increment (int);
//  void _add_hit_type (const char *);
//  void _add_hit_type (const IWString &);
    void _add_hit_type (int, const const_IWSubstring &);

  public:
    Demerit ();

    int ok () const;
    int debug_print (ostream &) const;

    int score () const { return _score;};
    int number_different_demerits_applied () const { return _number_different_demerits_applied;}

    const IWString & types () const { return _types;}

//  int extra (int, const char *);
//  int extra (int, const IWString &);
    int extra (int, const const_IWSubstring);
    int reject (const const_IWSubstring);
//  int reject (const char *);
//  int reject (const IWString &);
    int rejected () const;
    int rejected_by_single_rule () const { return 1 == _number_different_demerits_applied;}

    int write_in_tdt_form (ostream &) const;
};

extern void set_rejection_threshold (int);
extern int  rejection_threshold ();

/*
  Sept 2004. Dan Robertson needs to get all reasons separated with
  their individual demerit values. If that's the case, we build the
  _TYPE variable differently
*/

extern void set_demerit_reason_contains_individual_demerits (int);
extern void set_store_demerit_reasons_like_tsubstructure (int);

#endif
