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
#include <stdlib.h>

#include "substructure.h"
#include "target.h"
#include "misc2.h"

Elements_Needed::Elements_Needed ()
{
  _z = INVALID_ATOMIC_NUMBER;

  return;
}

Elements_Needed::Elements_Needed (atomic_number_t s) : _z (s)
{
}

int
Elements_Needed::ok () const
{
  return 1;
}

int
Elements_Needed::debug_print (ostream & os, const IWString & ind) const
{
//os << ind << "Elements_Needed: z = " << _z << ' ';
  return Min_Max_Specifier<int>::debug_print (os);
}

int
Elements_Needed::matches (Query_Atoms_Matched & qam) const
{
  if (INVALID_ATOMIC_NUMBER == _z)
  {
    cerr << "Elements_Needed::matches: cannot match invalid atom number\n";
    debug_print (cerr, "");
    iwabort ();
    return 0;
  }

  int esize = qam.number_elements ();

  int hits = 0;
  for (int i = 0; i < esize; i++)
  {
    const Substructure_Atom * a = qam[i];

    if (_z == a->current_hold_atom ()->atomic_number ())
      hits++;
  }

  return Min_Max_Specifier<int>::matches (hits);
}

int
Elements_Needed::matches (Molecule_to_Match & target_molecule) const
{
  int count = target_molecule.atoms_with_atomic_number (_z);

#ifdef DEBUG_ELEMENTS_NEEDED_MATCHES
  cerr << "Elements_Needed::matches: target contains " << count << " returning " << Min_Max_Specifier<int>::matches (count) << endl;
#endif

  return Min_Max_Specifier<int>::matches (count);
}
