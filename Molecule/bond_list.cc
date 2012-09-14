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
#include <assert.h>
#include <memory>

using namespace std;

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#include "bond_list.h"

Bond_list::Bond_list ()
{
  _magic = BOND_LIST_MAGIC;
}

Bond_list::~Bond_list ()
{
  _magic = 0;
}

int
Bond_list::ok () const
{
  if (BOND_LIST_MAGIC != _magic)
    return 0;

  if (! resizable_array_p<Bond>::ok ())
    return 0;

  return 1;
}

int
Bond_list::debug_print (ostream & os) const
{
  assert (os.good ());

  os << "Bond list contains " << _number_elements  << " members\n";

  for (int i = 0; i < _number_elements; i++)
  {
    os << " Bond " << i << " ";
    const Bond *b = _things[i];
    b->debug_print (os);
    os << endl;
  }

  return 1;
}

int
Bond_list::which_bond (atom_number_t a1, atom_number_t a2) const
{
  assert (ok ());
  assert (a1 >= 0 && a2 >= 0 && a1 != a2);

  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->involves (a1, a2))
      return i;
  }

  return -1;
}

int
Bond_list::remove_bond_between_atoms (atom_number_t a1, atom_number_t a2)
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->involves (a1, a2))
    {
      remove_item (i);
      return 1;
    }
  }

  cerr << "Bond_list::remove_bond_between_atoms: no bond between atoms " << a1 << " and " << a2 << endl;
  assert (NULL == "this should not happen");
  return 0;
}

/*
  Change the bond list to reflect the fact that atoms A1 and A2
  have been swapped. 
*/

int
Bond_list::swap_atoms (atom_number_t a1, atom_number_t a2)
{
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    Bond * b = _things[i];
    rc += b->swap_atoms (a1, a2);
  }

  return rc;
}

int
Bond_list:: move_atom_to_end_of_atom_list (atom_number_t zatom, int atoms_in_molecule)
{
  for (int i = 0; i < _number_elements; i++)
  {
    Bond * b = _things[i];

    atom_number_t a1 = b->a1 ();
    atom_number_t a2 = b->a2 ();

    if (a1 == zatom)
      a1 = atoms_in_molecule - 1;
    else if (a1 > zatom)
      a1--;

    if (a2 == zatom)
      a2 = atoms_in_molecule - 1;
    else if (a2 > zatom)
      a2--;

    b->set_a1a2 (a1, a2);
  }

  return 1;
}

#ifdef BONDS_KNOW_RING_MEMBERSHIP

int
Bond_list::invalidate_ring_info ()
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->invalidate_nrings ();
  }

  return 1;
}

#endif

int
Bond_list::copy_bond_types (bond_type_t * b) const
{
  assert (NULL != b);

  for (int i = 0; i < _number_elements; i++)
  {
    b[i] = BOND_TYPE_ONLY (_things[i]->btype ());
  }

  return _number_elements;
}

/*
  Optimisation!! If the first bond doesn't have a number assigned, we assume
  that none do. Dangerous, but for efficiency
*/

void
Bond_list::invalidate_bond_numbers ()
{
  if (0 == _number_elements)
    return;

  if (! _things[0]->bond_number_assigned ())
    return;

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->invalidate_bond_number ();
  }

  return;
}

int
Bond_list::set_all_bond_types (bond_type_t bt)
{
  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_not_directional();   // if all bonds are being set to one type, directionality will be lost

    if (bt == _things[i]->btype ())
      continue;

    _things[i]->set_bond_type (bt);
    rc++;
  }

  return rc;
}

int
Bond_list::cis_trans_bonds_present() const
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->is_directional())
      return 1;
  }

  return 0;
}
