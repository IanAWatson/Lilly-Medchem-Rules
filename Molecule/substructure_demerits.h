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
#ifndef SUBSTRUCTURE_DEMERIT_H
#define SUBSTRUCTURE_DEMERIT_H

class Molecule;
class Demerit;
class Charge_Assigner;

namespace substructure_demerits
{
extern void set_verbose (int);

extern void set_keep_going_after_rejection (int);

/*
  Feb 2005. People want to be able to apply just the rejections from here
*/

extern void set_only_apply_rejection_rules ();

extern int hard_coded_queries_statistics (ostream &);
extern int initialise_hard_coded_queries_to_do (IWString &);

extern int hard_coded_queries (Molecule &, Demerit &);

/*
  We need a means of passing the charge assigner in substructure_demerits.cc back
  to the calling programme so it can be initialised from the command line. Awful!
*/

extern Charge_Assigner & charge_assigner ();
};

#endif
