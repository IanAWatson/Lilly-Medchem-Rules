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

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#include "temp_detach_atoms.h"

#include "molecule.h"

Temp_Detach_Atoms::Temp_Detach_Atoms ()
{
  _active = 1;

  _remove_hydrogens_no_longer_needed = 1;

  _matoms = 0;

  _connection = NULL;

  _bt = SINGLE_BOND;

  return;
}

Temp_Detach_Atoms::~Temp_Detach_Atoms ()
{
  if (NULL != _connection)
    delete [] _connection;

  return;
}

int
Temp_Detach_Atoms::recognised_directive (const const_IWSubstring & token)
{
  if ("noremove" == token)
  {
    _remove_hydrogens_no_longer_needed = 0;
    return 1;
  }
  
  if ("nodetach" == token)
  {
    _active = 0;
    return 1;
  }

  return 0;
}

void
Temp_Detach_Atoms::do_not_reattach_to_atom (atom_number_t a)
{
  assert (a >= 0 && a < _matoms);
  assert (NULL != _connection );

  _connection[a] = -1;

  return;
}

/*
  Break any bond between a singly bonded Z and its attached atom. Record
  the identify of the neighbouring atom in the _CONNECTION array
*/

int
Temp_Detach_Atoms::detach_atoms (Molecule & m, atomic_number_t z)
{
  int matoms = m.natoms ();

  if (0 == matoms)
    return 0;

  if (matoms > _matoms)
  {
    if (NULL != _connection)
      delete [] _connection;

    _connection = new int[matoms];
  }

  _matoms = matoms;

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (z == m.atomic_number (i) && 1 == m.ncon (i))
    {
      atom_number_t o = m.other (i, 0);

      _connection[i] = o;

      m.remove_bond_between_atoms (i, o);

      rc++;
    }
    else
      _connection[i] = INVALID_ATOM_NUMBER;
  }

  _need_to_reattach = rc;

  return rc;
}

int
Temp_Detach_Atoms::reattach_atoms (Molecule & m)
{
  assert (_matoms == m.natoms ());

  if (0 == _need_to_reattach)
    return 1;      // don't need to do anything

  Set_of_Atoms atoms_to_be_removed;
  atoms_to_be_removed.resize (_matoms);

  for (int i = 0; i < _matoms; i++)
  {
    atom_number_t o = _connection[i];

    if (o < 0)
      continue;

    if (m.implicit_hydrogens (o))
      m.add_bond (i, o, _bt);
    else
      atoms_to_be_removed.add (i);
  }

  if (_remove_hydrogens_no_longer_needed)
    m.remove_atoms (atoms_to_be_removed);

  return 1;
}

