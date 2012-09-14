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
#include "assert.h"
#include <iomanip>
using namespace std;

#include "iwminmax.h"
#include "misc.h"

// get the private functions

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#define COMPILING_MOLECULED
#define COMPILING_SMILES_CC

class Smiles_First_Atom;
class CRDM_args;

#include "molecule.h"
#include "pearlman.h"
#include "path.h"
#include "misc2.h"

int
Molecule::_initialise_distance_matrix ()
{
  assert (NULL == _distance_matrix);
  assert (_number_elements > 0);

  _distance_matrix = new_int (_number_elements * _number_elements);

  return 1;
}

/*
  We store the distance matrix in "mostly" upper triangular form,
  so this central routine assumes a2 > a1
*/

int
Molecule::_bonds_between (atom_number_t a1, atom_number_t a2)
{
  assert (a1 < a2);

  if (NULL == _distance_matrix)
    _initialise_distance_matrix ();

  int * row = &_distance_matrix[_number_elements * a1];

//#define DEBUG_BONDS_BETWEEN
#ifdef DEBUG_BONDS_BETWEEN
  cerr << "Qbonds_between: between " << a1 << " and " << a2 << endl;

  int precision = 2;
  if (_number_elements > 9)
    precision = 3;
  else
    precision = 4;

  cerr << "MX is ";
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << setw (precision) << row[i];
  }
  cerr << endl;
#endif

// If it is already known, or is on the leading edge, grab it.

//cerr << "Atom " << a1 << " to " << a2 << " row = " << row[a2] << endl;

  if (row[a2] > 0)
    return row[a2];
  else if (row[a2] < 0)
    return - row[a2];

//cerr << "Yipes, DM incomplete, atoms " << a1 << " and " << a2 << ", d = " << row[a2] << " continuing...\n";

// The distance has not been computed. Work it out. Identify the
// most positive, negative distance along the row (if present).

  const int invalid_dist_value = - (nedges () + 1);   // longer than longest path in molecule

  int dist = invalid_dist_value;
  for (int i = 0; i < _number_elements; i++)
  {
    if (row[i] < 0 && row[i] > dist)
      dist = row[i];
  }

// If there were no negative numbers along the row, initialise some

  if (invalid_dist_value == dist)
  {
    const Atom * a = _things[a1];

    int a1con = a->ncon ();
    for (int i = 0; i < a1con; i++)
    {
      atom_number_t j = a->other (a1, i);
      row[j] = -1;
    }
    dist = -1;

    if (-1 == row[a2])
      return 1;
  }

  int nb = nedges ();

  while (1)
  {
    int positive_dist = - dist;
    for (int i = 0; i < _number_elements; i++)
    {
      if (row[i] != dist)
        continue;

#ifdef DEBUG_BONDS_BETWEEN
      cerr << "At distance " << dist << " processing " << i << endl;
#endif

      const Atom * a = _things[i];

      int icon = a->ncon ();
      for (int j = 0; j < icon; j++)
      {
        atom_number_t k = a->other (i, j);

#ifdef DEBUG_BONDS_BETWEEN
        cerr << "  Attached to atom " << j << " current = " << row[k] << endl;
#endif       
        if (0 == row[k])
          row[k] = dist - 1;
        else if (row[k] > positive_dist + 1)
          row[k] = dist - 1;
        else if (row[k] < dist - 1)
          row[k] = dist - 1;

        if (k == a2)
          return positive_dist + 1;
      }

      row[i] = positive_dist;    // only change it when all connections processed.
    }

    dist--;
    if (- dist > nb)
    {
      cerr << "Fatal, error, cannot find dist " << a1 << " to " << a2 << endl;
      iwabort ();
    }
  }
}

int
Molecule::bonds_between (atom_number_t a1, atom_number_t a2)
{
//assert (ok_2_atoms (a1, a2));

//cerr << "Molecule::bonds_between: atoms " << a1 << " and " << a2 << " dm = " << _distance_matrix << endl;

  if (NULL == _distance_matrix)
    _initialise_distance_matrix ();

// The atoms must be in the same fragment

  if (! _fragment_information.contains_valid_data ())
    (void) number_fragments ();

//assert (_fragment_information.fragment_membership (a1) == _fragment_information.fragment_membership (a2));
  if( _fragment_information.fragment_membership (a1) != _fragment_information.fragment_membership (a2) )
    return ATOMS_NOT_BONDED;

  if (a1 > a2)
    return _bonds_between (a2, a1);
  else
    return _bonds_between (a1, a2);
}

//#define DEBUG_ATOMS_BETWEEN

int
Molecule::atoms_between (atom_number_t a1,
                         atom_number_t a2,
                         Set_of_Atoms & s)
{
  int d = bonds_between (a1, a2);

  if (1 == d)
  {
    s.resize (0);
    return 0;
  }

  if (s.number_elements ())
    s.resize_keep_storage (0);
  else
    s.resize (d - 1);

#ifdef DEBUG_ATOMS_BETWEEN
  cerr << "Molecule::atoms_between:atoms " << a1 << " '" << smarts_equivalent_for_atom (a1) << "' and " << a2 << " '" << smarts_equivalent_for_atom (a2) << " are " << d << " bonds apart\n";
#endif

  return _atoms_between (a1, a2, d - 1, s);
}

int
Molecule::_atoms_between (atom_number_t a1,
                          atom_number_t a2,
                          int distance_needed,
                          Set_of_Atoms & s)
{
#ifdef DEBUG_ATOMS_BETWEEN
  cerr << "Molecule::_atoms_between:contine to atom " << a1 << " '" << smarts_equivalent_for_atom (a1) << "' and " << a2 << " '" << smarts_equivalent_for_atom (a2) << "' d = " << distance_needed << endl;
#endif

  assert (NULL != _distance_matrix);

  const Atom * a = _things[a1];

  int acon = a->ncon ();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other (a1, i);

    if (j == a2)     // done, we got to A2
      return 1;

#ifdef DEBUG_ATOMS_BETWEEN
    cerr << "Distance between " << j << " and " << a2 << " is " << _distance_matrix[j * _number_elements + a2] << endl;
#endif

    int d;
    if (j < a2)
      d = _bonds_between (j, a2);
    else
      d = _bonds_between (a2, j);

    if (d != distance_needed)
      continue;

    s.add (j);

    return 1 + _atoms_between (j, a2, distance_needed - 1, s);
  }

  cerr << "Molecule::_atoms_between:yipes, from " << a1 << " '" << smarts_equivalent_for_atom (a1) << "' nothing " << distance_needed << " bonds to " << a2 << " '" << smarts_equivalent_for_atom (a2) << "'\n";
  iwabort ();

  return 0;
}

int
Molecule::longest_path ()
{
  iwmax<int> rc (0);

  for (int i = 0; i < _number_elements; i++)
  {
    for (int j = i + 1; j < _number_elements; j++)
    {
      rc.extra (_bonds_between (i, j));
    }
  }

  return rc.maxval ();
}

//#define DEBUG_COMPUTE_ROW_DM





/*
  Fastest version uses the bond list
*/

void
Molecule::_compute_distance_matrix ()
{
  int nb = _bond_list.number_elements();

  set_vector(_distance_matrix, _number_elements * _number_elements, _number_elements + _number_elements);

  const Bond * const * allbonds = _bond_list.rawdata();

  for (int i = 0; i < _number_elements; i++)   // diagonals are zero
  {
    _distance_matrix[i * _number_elements + i] = 0;
  }

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = allbonds[i];

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    int * row1 = _distance_matrix + (a1 * _number_elements);
    int * row2 = _distance_matrix + (a2 * _number_elements);

    row1[a2] = 1;
    row2[a1] = 1;
  }

  while (1)
  {
    int keep_going = 0;

    for (int i = 0; i < nb; i++)
    {
      const Bond * b = allbonds[i];

      atom_number_t a1 = b->a1();
      atom_number_t a2 = b->a2();

      int * row1 = _distance_matrix + (a1 * _number_elements);
      int * row2 = _distance_matrix + (a2 * _number_elements);

      int tmp;
      for (int j = 0; j < _number_elements; j++)
      {
        tmp = row2[j] + 1;
        if (tmp < row1[j])
        {
          row1[j] = tmp;
          keep_going = 1;
        }
        tmp = row1[j] + 1;
        if (tmp < row2[j])
        {
          row2[j] = tmp;
          keep_going = 1;
        }
      }
    }

#ifdef DEBUG_DISTANCE_MATRIX
    cerr << "End of cycle, kg " << keep_going << endl;
    for (int j = 0; j < _number_elements; j++)
    {
      cerr << " Atom " << j << ':';
      const int * r = _distance_matrix + (j * _number_elements);

      for (int k = 0; k < _number_elements; k++)
      {
        cerr << ' ' << r[k];
      }
      cerr << endl;
    }
#endif

    if (0 == keep_going)
      return;
  }
}

//#define CHECK_DISTANCE_MATRIX

int
Molecule::recompute_distance_matrix ()
{
  if (! _fragment_information.contains_valid_data ())
    (void) number_fragments ();

  if (NULL == _distance_matrix)
    _distance_matrix = new_int (_number_elements * _number_elements, _number_elements + 9);
  else
    set_vector (_distance_matrix, _number_elements * _number_elements, _number_elements + 9);

  if (1 == _number_elements)
    return 1;

  if (2 == _number_elements)
  {
    _distance_matrix[1] = 1;
    _distance_matrix[3] = 1;

    return 1;
  }

  _compute_distance_matrix();

  int rc = 1;

#ifdef CHECK_DISTANCE_MATRIX
  for (int i = 0; i < _number_elements; i++)
  {
    if (0 != _distance_matrix[_number_elements * i + i])
    {
      cerr << "Non zero diagonal on distance matrix, i = " << i << " value = " << _distance_matrix[_number_elements * i + i] << endl;
      rc = 0;
    }

    for (int j = i + 1; j < _number_elements; j++)
    {
      if (_distance_matrix[_number_elements * i + j] != _distance_matrix[_number_elements * j + i])
      {
        cerr << "Distance matrix error, from " << i << " to " << j << " is " << _distance_matrix[_number_elements * i + j] << endl;
        cerr << "Distance matrix error, from " << j << " to " << i << " is " << _distance_matrix[_number_elements * j + i] << endl;
        rc = 0;
      }
      else if (0 == _distance_matrix[_number_elements * i + j])
      {
        cerr << "Zero distance matrix entry " << i << ", " << j << endl;
        rc = 0;
      }
      else if (_distance_matrix[_number_elements * i + j] < 0)
      {
        cerr << "Incomplete distance matrix determinations, atoms " << i << " and " << j << endl;
      }
      else
      {
//      cerr << "Bonds between " << i << " and " << j << " is " << _distance_matrix[_number_elements * i + j] << endl;
      }
    }
  }

  if (0 == rc)
    iwabort ();
#endif

  return rc;
}

// arch-tag: c8079359-e517-4d71-8303-1697d2ab7377
