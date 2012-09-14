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
#ifndef IS_ACTUALLY_CHIRAL_H
#define IS_ACTUALLY_CHIRAL_H

#include "iwmtypes.h"
#include "path_scoring.h"

class Molecule;

extern int is_actually_chiral (Molecule & m, atom_number_t);
extern int is_actually_chiral (Molecule & m, atom_number_t, resizable_array_p<Path_Scoring> &);

// Returns the number of invalid chiral centres removed

extern int do_remove_invalid_chiral_centres (Molecule & m);

extern void set_max_iterations (int m);

#endif
