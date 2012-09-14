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
#ifndef RWMOLECULE_H
#define RWMOLECULE_H

#define EXTRA_STRING_RECORD(ds, b, c) \
  if (! (ds).next_record ((b)))\
  {\
    cerr << (c) << " eof\n";\
    return 0;\
  }

extern int rwmolecule_error (const char *, iwstring_data_source &);

extern int write_coordinates (ostream &, const Atom * a, int = 0);

extern int append_sybyl_atom_type (ostream & os, int atype, const IWString & asymbol);

extern void reset_rwmolecule_file_scope_variables ();

#endif
