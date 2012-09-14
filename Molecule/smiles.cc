/**************************************************************************

    Copyright (C) 2012  Eli Lilly and Company

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

#include "misc.h"
#include "iw_auto_array.h"
#include "iwbits.h"
#include "iwrandom.h"


// Define this to get the molecule private functions defined here

#define COMPILING_MOLECULER_CC
#define COMPILING_SMILES_CC
#define COMPILING_CTB

#include "molecule.h"
#include "path.h"
#include "smiles.h"
#include "aromatic.h"
#include "misc2.h"
#include "iwrnm.h"
#include "iwrcb.h"

//#define USE_IWMALLOC
#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

/*
  I need an object to describe how to choose the first atom in a smiles
*/

#define SMILES_FIRST_ATOM_DEFAULT 1
#define SMILES_FIRST_ATOM_NICE 2
#define SMILES_FIRST_ATOM_RANDOM 3
#define SMILES_FIRST_ATOM_UNIQUE 4

/*
  Note that if 
*/

class Smiles_First_Atom
{
  private:
    int _how_to_choose;    // one of the defined values above

    atom_number_t _a;

  public:
    Smiles_First_Atom ();

    int smdeflt () const { return SMILES_FIRST_ATOM_DEFAULT == _how_to_choose;}
    int nice    () const { return SMILES_FIRST_ATOM_NICE    == _how_to_choose;}
    int random  () const { return SMILES_FIRST_ATOM_RANDOM  == _how_to_choose;}
    int unique  () const { return SMILES_FIRST_ATOM_UNIQUE  == _how_to_choose;}

    void set_build_type (int s) { _how_to_choose = s;}

    void set_atom_number (atom_number_t a) { _a = a;}

    int atom_specified (atom_number_t &) const;
    int atom_specified () const { return INVALID_ATOM_NUMBER != _a;}
    void unset_first_atom () { _a = INVALID_ATOM_NUMBER; _how_to_choose = SMILES_FIRST_ATOM_DEFAULT;}
};

Smiles_First_Atom::Smiles_First_Atom ()
{
  _how_to_choose = SMILES_FIRST_ATOM_DEFAULT;

  _a = INVALID_ATOM_NUMBER;

  return;
}

int
Smiles_First_Atom::atom_specified (atom_number_t & first_atom) const
{
  if (INVALID_ATOM_NUMBER == _a)
    return 0;

  first_atom = _a;

  return 1;
}

//#define DEBUG_SMILES_CHOOSE_FIRST_ATOM
#ifdef DEBUG_SMILES_CHOOSE_FIRST_ATOM
static ostream &
operator << (ostream & os, const Smiles_First_Atom & smfa)
{
  os << "SMFA: ";
  if (smfa.atom_specified ())
  {
    atom_number_t a;
    (void) smfa.atom_specified (a);
    os << "atom " << a;
  }
  else if (smfa.smdeflt ())
    os << "default";
  else if (smfa.nice())
    os << "nice";
  else if (smfa.unique())
    os << "unique";
  else if (smfa.random())
    os << "random";
  else 
    os << " HUH";

  return os;
}
#endif

static int _write_aromatic_bonds_as_colons = 0;

void
set_write_smiles_aromatic_bonds_as_colons (int s)
{
  _write_aromatic_bonds_as_colons = s;
}

int
write_smiles_aromatic_bonds_as_colons ()
{
  return _write_aromatic_bonds_as_colons;
}

void
Smiles_Information::_default_values ()
{
  _smiles_order_type = INVALID_SMILES_ORDER_TYPE;

  _smiles_order = NULL;

  _smiles_is_smarts = 0;

  _create_smarts_embedding = NULL;

  _user_specified_atomic_smarts = NULL;

  return;
}

Smiles_Information::Smiles_Information (int natoms) : _natoms(natoms)
{
  _default_values();

  return;
};

Smiles_Information::Smiles_Information () : _natoms(-1)
{
  _default_values();

  return;
}

Smiles_Information::~Smiles_Information ()
{
  if (NULL != _smiles_order)
    delete [] _smiles_order;

  if (NULL != _create_smarts_embedding)
    delete [] _create_smarts_embedding;

  if (NULL != _user_specified_atomic_smarts)
    delete [] _user_specified_atomic_smarts;

  return;
}

int
Smiles_Information::debug_print (ostream & os) const
{
  if (NULL != _smiles_order)
  {
    if (UNIQUE_SMILES_ORDER_TYPE == _smiles_order_type)
      os << "Unique smiles order computed\n";
    else if (RANDOM_SMILES_ORDER_TYPE == _smiles_order_type)
      os << "Random smiles order computed\n";
    else if (DEFAULT_SMILES_ORDER_TYPE == _smiles_order_type)
      os << "Default smiles order computed\n";
    else if (SUBSET_SMILES_ORDER_TYPE == _smiles_order_type)
      os << "Subset smiles order computed\n";
    else
      os << "Hmmm, smiles order type is " << _smiles_order_type << endl;
  }

  _ring_closure_bonds.write_bonds (os);

  if (_smiles_start_atom.number_elements() > 0)
  {
    for (int i = 0; i < _smiles_start_atom.number_elements(); i++)
    {
      os << "smiles in fragment " << i << " starts with atom " << _smiles_start_atom[i] << endl;
    }
  }

  if (_smiles.length())
    os << "Smiles is '" << _smiles << "'\n";

  return os.good();
}

void
Smiles_Information::make_empty ()
{
  _smiles = EMPTY_MOLECULE_SMILES;

  if (NULL != _smiles_order)
  {
    delete [] _smiles_order;
    _smiles_order = NULL;
  }

  return;
}

int
Smiles_Information::prepare_to_build_ordering (int matoms)
{
  _smiles_order_type = INVALID_SMILES_ORDER_TYPE;

  if (NULL == _smiles_order)
    _smiles_order = new_int(matoms, -1);
  else
    set_vector (_smiles_order, matoms, -1);

  _natoms = matoms;

  if (! _ring_closure_bonds.activate(matoms))
    return 0;

  _smiles_start_atom.resize_keep_storage(0);

  return 1;
}


int
Smiles_Information::prepare_to_build_smiles (int matoms)
{
  _smiles.resize_keep_storage(0);

  if (_smiles.elements_allocated() < 3 * matoms)
    _smiles.resize(3 * matoms);

  if (_atom_order_in_smiles.number_elements() > 0)
    _atom_order_in_smiles.resize_keep_storage(0);
  else
    _atom_order_in_smiles.resize(matoms);

  _natoms = matoms;

  return 1;
}

void
Smiles_Information::invalidate ()
{
  _smiles_order_type = INVALID_SMILES_ORDER_TYPE;

  _smiles.resize_keep_storage(0);

  if (NULL != _smiles_order)
  {
    delete [] _smiles_order;
    _smiles_order = NULL;
  }

  return;
}

int
Smiles_Information::create_smarts_embedding (atom_number_t zatom) const
{
  if (NULL == _create_smarts_embedding)  // should be a fatal error
    return 0;

  return _create_smarts_embedding[zatom];
}

int
Smiles_Information::set_create_smarts_embedding (int s)
{
  assert (_natoms > 0);

  if (NULL == _create_smarts_embedding)
    _create_smarts_embedding = new_int(_natoms, s);
  else
    set_vector(_create_smarts_embedding, _natoms, s);

  return 1;
}

int
Smiles_Information::set_create_smarts_embedding (atom_number_t zatom,
                                                 int s)
{
  assert (_natoms > 0 && zatom >= 0 && zatom < _natoms);

  if (NULL == _create_smarts_embedding)
    _create_smarts_embedding = new_int(_natoms);

  _create_smarts_embedding[zatom] = s;

  return 1;
}

/*
*/

int
Molecule::_smiles_choose_first_atom (const int * zorder,
                                     Smiles_First_Atom & smfa,
                                     atom_number_t & first_atom,
                                     const int * include_atom)
{
#ifdef DEBUG_SMILES_CHOOSE_FIRST_ATOM
  cerr << "INto _smiles_choose_first_atom " << smfa << endl;
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << "Atom " << i << " order " << zorder[i] << endl;
  }
#endif

  if (smfa.atom_specified(first_atom) && zorder[first_atom] < 0)
  {
    smfa.unset_first_atom();
    return 1;
  }

  if (smfa.smdeflt())
    return _smiles_choose_first_atom(zorder, first_atom, include_atom);

  if (smfa.nice())
    return _smiles_choose_nice_first_atom(zorder, first_atom, include_atom);

  if (smfa.random())
    return _smiles_choose_random_first_atom(zorder, first_atom, include_atom);

  if (smfa.unique())
  {
    if (NULL == include_atom)
      return _smiles_choose_unique_first_atom(zorder, first_atom);
    else
      return _smiles_choose_unique_first_atom(zorder, first_atom, include_atom);
  }

  cerr << "Molecule::_smiles_choose_first_atom: no method specified\n";
  iwabort();

  return 0;
}

/*
  Note that this is incorrect, a chiral atom can start a smiles, but
  not sure if that works or not...

  Mar 2004. Ran into problems with a subset. Every member of the molecule subset
  was a chiral centre. Therefore save a possible match in ATOM_TO_RETURN_IF_NOTHING_ELSE_FOUND
*/

int
Molecule::_smiles_choose_first_atom (const int * zorder,
                                     atom_number_t & first_atom,
                                     const int * include_atom)
{
  int include_chiral_info = include_chiral_info_in_smiles();

  atom_number_t atom_to_return_if_nothing_else_found = INVALID_ATOM_NUMBER;

  for (int i = 0; i < _number_elements; i++)
  {
    if (zorder[i] >= 0)    // already done
      continue;

    if (NULL != include_atom && 0 == include_atom[i])
      continue;

//  A smiles cannot START at a chiral atom

#ifdef DEBUG_SMILES_CHOOSE_FIRST_ATOM
    cerr << "Can the smiles start with atom " << i << endl;
    cerr << "include_chiral_info " << include_chiral_info << endl;
    cerr << "ncon " << _things[i]->ncon() << endl;
    cerr << "chiral " << chiral_centre_at_atom(i) << endl;
    if (include_chiral_info && _things[i]->ncon() > 2 && chiral_centre_at_atom(i))
      cerr << "Nope, that looks chiral\n";
#endif

    if (include_chiral_info && _things[i]->ncon() > 2 && chiral_centre_at_atom(i))
    {
      atom_to_return_if_nothing_else_found = i;
      continue;
    }

    first_atom = i;
    return 1;
  }

  if (INVALID_ATOM_NUMBER == atom_to_return_if_nothing_else_found)
    return 0;

  first_atom = atom_to_return_if_nothing_else_found;

  return 1;
}

int
Molecule::_smiles_choose_nice_first_atom (const int * zorder,
                               atom_number_t & first_atom,
                               const int * include_atom)
{
  int include_chiral_info = include_chiral_info_in_smiles();

  atom_number_t zdefault = INVALID_ATOM_NUMBER;
  atom_number_t best_singly_connected = INVALID_ATOM_NUMBER;
  atom_number_t best_singly_connected_non_ring = INVALID_ATOM_NUMBER;
  atom_number_t best_doubly_connected = INVALID_ATOM_NUMBER;
  atom_number_t best_doubly_connected_non_ring = INVALID_ATOM_NUMBER;
  atom_number_t best_non_ring = INVALID_ATOM_NUMBER;

  for (int i = 0; i < _number_elements; i++)
  {
    if (zorder[i] >= 0)    // already done
      continue;

    if (NULL != include_atom && 0 == include_atom[i])
      continue;

//  A smiles cannot START at a chiral atom

    if (include_chiral_info && chiral_centre_at_atom(i))
      continue;

    if (INVALID_ATOM_NUMBER == zdefault)
      zdefault = i;

    const Atom * a = _things[i];

    int icon = a->ncon();

    int nr = is_ring_atom(i);

    if (0 == nr && INVALID_ATOM_NUMBER == best_non_ring)
      best_non_ring = i;

    if (1 == icon)
    {
      if (0 == nr && INVALID_ATOM_NUMBER == best_singly_connected_non_ring)
        best_singly_connected_non_ring = i;
      else if (INVALID_ATOM_NUMBER == best_singly_connected)
       best_singly_connected = i;
    }
    else if (2 == icon)
    {
      if (0 == nr && INVALID_ATOM_NUMBER == best_doubly_connected_non_ring)
        best_doubly_connected_non_ring = i;
      else if (INVALID_ATOM_NUMBER == best_doubly_connected)
       best_doubly_connected = i;
    }
  }

  if (INVALID_ATOM_NUMBER != best_singly_connected_non_ring)
  {
    first_atom = best_singly_connected_non_ring;
    return 1;
  }

  if (INVALID_ATOM_NUMBER != best_doubly_connected_non_ring)
  {
    first_atom = best_doubly_connected_non_ring;
    return 1;
  }

  if (INVALID_ATOM_NUMBER != best_non_ring)
  {
    first_atom = best_non_ring;
    return 1;
  }

  if (INVALID_ATOM_NUMBER != best_doubly_connected_non_ring)
  {
    first_atom = best_doubly_connected_non_ring;
    return 1;
  }

  if (INVALID_ATOM_NUMBER != best_doubly_connected)
  {
    first_atom = best_doubly_connected;
    return 1;
  }

  if (INVALID_ATOM_NUMBER != zdefault)
  {
    first_atom = zdefault;
    return 1;
  }

  return 0;
}

static Random_Number_Working_Storage smiles_random_number_stream;

/*
  External entry point for setting the random number used for
  random smiles generation
*/

void
set_smiles_random_number_seed (random_number_seed_t seed)
{
  smiles_random_number_stream.set_seed(seed);
}

random_number_seed_t
set_smiles_random_number_seed_random ()
{
  return smiles_random_number_stream.choose_random_seed();
}

/*
  As you can see, this is really not truly random, in that the
  first atoms will be favoured.
*/

int
Molecule::_smiles_choose_random_first_atom (const int * zorder,
                    atom_number_t & first_atom,
                    const int * include_atom)
{
  int include_chiral_info = include_chiral_info_in_smiles();

  int istart = smiles_random_number_stream.intbtwij(0, _number_elements);

  atom_number_t atom_to_return_if_nothing_else_found = INVALID_ATOM_NUMBER;

  for (int i = istart; i < _number_elements; i++)
  {
    if (zorder[i] >= 0)    // already done
      continue;

    if (NULL != include_atom && 0 == include_atom[i])
      continue;

//  A smiles cannot START at a chiral atom

    if (include_chiral_info && chiral_centre_at_atom(i))
    {
      atom_to_return_if_nothing_else_found = i;
      continue;
    }

    first_atom = i;
    return 1;
  }

// If we come to here, we did not find a suitable atom in the
// range [istart..matoms). How about (0..istart)

  for (int i = 0; i < istart; i++)
  {
//  cerr << "Checking random start atom " << i << " zorder = " << zorder[i] << endl;
    if (zorder[i] >= 0)    // already processed
      continue;

//  A smiles cannot START at a chiral atom

    if (include_chiral_info_in_smiles() && chiral_centre_at_atom(i))
    {
      atom_to_return_if_nothing_else_found = i;
      continue;
    }

    first_atom = i;
    return 1;
  }

  if (INVALID_ATOM_NUMBER != atom_to_return_if_nothing_else_found)
  {
    first_atom = atom_to_return_if_nothing_else_found;
    return 1;
  }

//cerr << "No random start atom available\n";

// Nothing possible

  return 0;
}

/*
  For unique smiles, we choose the highest ranked, most lowly connected atom
*/

int
Molecule::_smiles_choose_unique_first_atom (const int * zorder,
                                            atom_number_t & first_atom)
{
  const int * canonical_rank = _symmetry_class_and_canonical_rank.canonical_rank();

  assert (NULL != canonical_rank);

  int include_chiral_info = include_chiral_info_in_smiles();

  int min_ncon = _number_elements;    // greater than all ncon() values

  atom_number_t zdefault = INVALID_ATOM_NUMBER;
  int rsave = 0;

  atom_number_t chiral_atom_if_we_need_to = INVALID_ATOM_NUMBER;   // in case everything chiral [P@]12[P@]3[P@@]1[P@@]23

  for (int i = 0; i < _number_elements; i++)
  {
    if (zorder[i] >= 0)    // already processed
      continue;

//  A smiles cannot START at a chiral atom

    if (include_chiral_info && chiral_centre_at_atom(i))
    {
      if (chiral_atom_if_we_need_to < 0)
        chiral_atom_if_we_need_to = i;
      continue;
    }

    const Atom * a = _things[i];

    int ncon = a->ncon();

    if (ncon < min_ncon)
    {
      min_ncon = ncon;
      zdefault = i;
      rsave = canonical_rank[i];
    }
    else if (ncon == min_ncon && canonical_rank[i] > rsave)
    {
      zdefault = i;
      rsave = canonical_rank[i];
    }
  }

  if (INVALID_ATOM_NUMBER != zdefault)
  {
    first_atom = zdefault;
    return 1;
  }

  if (chiral_atom_if_we_need_to >= 0)
  {
    first_atom = chiral_atom_if_we_need_to;
    return 1;
  }

  return 0;
}

/*
  When dealing with a subset, we need to take into account only the number
  of atoms connected in the subset
*/

int
Molecule::_smiles_choose_unique_first_atom (const int * zorder,
                                            atom_number_t & first_atom,
                                            const int * include_atom)
{
  const int * canonical_rank = _symmetry_class_and_canonical_rank.canonical_rank();

  assert (NULL != canonical_rank);

  assert (NULL != include_atom);

  int include_chiral_info = include_chiral_info_in_smiles();

  int min_ncon = _number_elements;    // greater than all ncon() values

  atom_number_t zdefault = INVALID_ATOM_NUMBER;
  int rsave = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (zorder[i] >= 0)    // already processed
      continue;

    if (0 == include_atom[i])
      continue;

//  A smiles cannot START at a chiral atom

    if (include_chiral_info && chiral_centre_at_atom(i))
      continue;

    const Atom * a = _things[i];

    int ncon = a->ncon(i, include_atom);

    if (ncon < min_ncon)
    {
      min_ncon = ncon;
      zdefault = i;
      rsave = canonical_rank[i];
    }
    else if (ncon == min_ncon && canonical_rank[i] > rsave)
    {
      zdefault = i;
      rsave = canonical_rank[i];
    }
  }

  if (INVALID_ATOM_NUMBER != zdefault)
  {
    first_atom = zdefault;
    return 1;
  }

  return 0;
}


/*
  The default (fast) behaviour is to find the first singly connected
  and follow that. Otherwise follow the first multiple bond
*/

int
Molecule::_smiles_choose_next_atom (const int * zorder,
                         atom_number_t current_atom,
                         atom_number_t & next_atom,
                         const int * include_atom)
{
  atom_number_t zdefault = INVALID_ATOM_NUMBER;

  const Atom * a = _things[current_atom];

  int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j;
    bond_type_t   bt;
    a->other_and_type(current_atom, i, j, bt);

    if (zorder[j] >= 0)    // atom already classified
      continue;

    if (NULL != include_atom && 0 == include_atom[j])
      continue;

    if (1 == _things[j]->ncon())
    {
      next_atom = j;
      return 1;
    }

    if (! IS_SINGLE_BOND(bt))
    {
      next_atom = j;
      return 1;
    }
    else if (INVALID_ATOM_NUMBER == zdefault)
      zdefault = j;
  }

  if (INVALID_ATOM_NUMBER != zdefault)
  {
    next_atom = zdefault;
    return 1;
  }

  return 0;
}

/*
  Same as smiles_choose_next_atom except that we preferentially follow
  any non ring path before following paths in a ring
*/

int
Molecule::_smiles_choose_nice_next_atom (const int * zorder,
                         atom_number_t current_atom,
                         atom_number_t & next_atom,
                         const int * include_atom)
{
  atom_number_t zdefault = INVALID_ATOM_NUMBER;
  atom_number_t non_ring_atom = INVALID_ATOM_NUMBER;
  atom_number_t multiple_bond = INVALID_ATOM_NUMBER;

  const Atom * a = _things[current_atom];

  int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j;
    bond_type_t   bt;
    a->other_and_type(current_atom, i, j, bt);

    if (zorder[j] >= 0)     // atom already classified
      continue;

    if (NULL != include_atom && 0 == include_atom[j])
      continue;

    if (1 == _things[j]->ncon())
    {
      next_atom = j;
      return 1;
    }

    if (INVALID_ATOM_NUMBER == zdefault)
      zdefault = j;

    if (is_non_ring_atom(j) && INVALID_ATOM_NUMBER == non_ring_atom)
      non_ring_atom = j;

    if (! IS_SINGLE_BOND (bt))    // always prefer multiple bonds
      multiple_bond = j;
  }

  if (INVALID_ATOM_NUMBER != non_ring_atom)
  {
    next_atom = non_ring_atom;
    return 1;
  }

  if (INVALID_ATOM_NUMBER != multiple_bond)
  {
    next_atom = multiple_bond;
    return 1;
  }

  if (INVALID_ATOM_NUMBER != zdefault)
  {
    next_atom = zdefault;
    return 1;
  }

  return 0;
}

int
Molecule::_smiles_choose_random_next_atom (const int * zorder,
                                atom_number_t current_atom,
                                atom_number_t & next_atom,
                                const int * include_atom)
{
  const Atom * c = _things[current_atom];

  int acon = c->ncon();

  atom_number_t zdefault = INVALID_ATOM_NUMBER;
  atom_number_t multiple_bond = INVALID_ATOM_NUMBER;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j;
    bond_type_t   bt;
    c->other_and_type(current_atom, i, j, bt);

    if (zorder[j] >= 0)
      continue;

    if (NULL != include_atom && 0 == include_atom[j])
      continue;

    if (! IS_SINGLE_BOND (bt))
    {
      if (INVALID_ATOM_NUMBER == multiple_bond)
        multiple_bond = j;
      else if (smiles_random_number_stream.random_one_or_zero())
        multiple_bond = j;
    }
    else if (INVALID_ATOM_NUMBER == zdefault)
      zdefault = j;
    else if (smiles_random_number_stream.random_one_or_zero())
      zdefault = j;
  }

  if (INVALID_ATOM_NUMBER != multiple_bond)
  {
    next_atom = multiple_bond;
    return 1;
  }

  if (INVALID_ATOM_NUMBER != zdefault)
  {
    next_atom = zdefault;
    return 1;
  }

  return 0;
}

//#define DEBUG_UNIQUE_SMILES_ORDERING

/*
  The rule for unique order is connection with the highest bond order
  first, with ties broken by unique ordering.
*/

int
Molecule::_smiles_choose_unique_next_atom (const int * zorder,
                                atom_number_t current_atom,
                                atom_number_t & next_atom,
                                const int * include_atom)
{
#ifdef DEBUG_UNIQUE_SMILES_ORDERING
  cerr << "Choosing next unique atom from " << current_atom << endl;
#endif

// Variables for bond order decisions

  atom_number_t zdefault = INVALID_ATOM_NUMBER;
  int hbc_save = 0;                 // initialised to shut gcc up
  int highest_bond_count = 0;     // definitely initialised to 0

  const Atom * a = _things[current_atom];

  int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(current_atom);

    if (zorder[j] >= 0)
      continue;

    if (NULL != include_atom && 0 == include_atom[j])
      continue;

    int bcount;

//  Note that the counting here is different from in the function number_of_bonds ().
//  In that function, single and double bonds are counted first, here we want to catch
//  aromatic bonds first as contributing one to the bond count

    if (b->is_aromatic())
      bcount = 1;
    else
      bcount = b->number_of_bonds();

#ifdef DEBUG_UNIQUE_SMILES_ORDERING
    cerr << "Choosing next unique atom, " << j << " (bcount = " << bcount << ", rj = " << canonical_rank(j) << ")\n";
#endif

    if (bcount < highest_bond_count)   // definitely not
      continue;

    int rj = canonical_rank(j);

    if (bcount > highest_bond_count)
    {
      highest_bond_count = bcount;
      zdefault = j;
      hbc_save = rj;
    }
    else if (rj > hbc_save)    // bcount == highest_bond_count
    {
#ifdef DEBUG_UNIQUE_SMILES_ORDERING
      cerr << "Bonds equal, rj = " << rj << " hbc_save = " << hbc_save << endl;
#endif

      zdefault = j;
      hbc_save = rj;
    }
  }

  if (INVALID_ATOM_NUMBER != zdefault)
  {
    next_atom = zdefault;

#ifdef DEBUG_UNIQUE_SMILES_ORDERING
    cerr << "Next unique atom is atom " << next_atom << endl;
#endif

    return 1;
  }

  return 0;
}

//#define DEBUG_BUILD_SMILES_ORDERING

/*
  Look for ring closure bonds attached to the atom, then continue
  the search.
*/

int
Molecule::_build_smiles_ordering (int (Molecule::*identify_next_atom) (const int *, atom_number_t, atom_number_t &, const int *),
                                  const atom_number_t previous_atom,
                                  const atom_number_t zatom,
                                  int & icounter,
                                  const int * include_atom,
                                  Smiles_Information & smi_info)
{
  int * zorder = smi_info.smiles_order();

  assert (zatom >= 0 && zatom < _number_elements && zorder[zatom] < 0);

#ifdef DEBUG_BUILD_SMILES_ORDERING
  cerr << "_build_smiles_ordering continues with atom " << zatom << endl;
#endif

  zorder[zatom] = icounter;

  Atom * a = _things[zatom];

  icounter++;

  int acon = a->ncon();

  if (INVALID_ATOM_NUMBER != previous_atom && 1 == acon)  // got to a terminal atom
    return 0;

  Set_of_Atoms unprocessed_connections;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);;

    if (previous_atom == j)      // the atom from which we came
      continue;

#ifdef DEBUG_BUILD_SMILES_ORDERING
  cerr << "  Atom " << j << " is connected";
  if (zorder[j] >= 0)
    cerr << ". Ring closure detected";
  cerr << endl;
#endif

    if (zorder[j] >= 0)     // we have found a ring closure
      smi_info.add_ring_closure_bond(zatom, j);
    else if (NULL != include_atom && 0 == include_atom[j])
      ;
    else
      unprocessed_connections.add(j);
  }

  int nu = unprocessed_connections.number_elements();

  if (0 == nu)
    return 0;

  if (1 == nu)
    return 1 + _build_smiles_ordering(identify_next_atom,
                 zatom, unprocessed_connections[0], icounter, include_atom, smi_info);

  atom_number_t b;
  int rc = 0;
  while ((this->*identify_next_atom)(zorder, zatom, b, include_atom))
  {
    rc += _build_smiles_ordering(identify_next_atom,
                                  zatom, b, icounter, include_atom, smi_info);
  }

  return rc;
}

/*
  All smiles building comes through here.
*/

int
Molecule::_build_smiles_ordering (Smiles_First_Atom & smfa,
                                  int (Molecule::* identify_next_atom) (const int *, atom_number_t, atom_number_t &, const int *),
                                  const int * include_atom,
                                  Smiles_Information & smi_info)
{
  if (! _fragment_information.contains_valid_data())
    (void) number_fragments();

//cerr << "Molecule::_build_smiles_ordering:there are " << _fragment_information.number_fragments() << " fragments in full molecule with " << _number_elements << " atoms\n";

  smi_info.prepare_to_build_ordering (_number_elements);

  int zcounter = 0;    // allocated in increasing order
  atom_number_t a;     // start atom within each fragment

  int * zorder = smi_info.smiles_order();

  int frag = 0;

  while (_smiles_choose_first_atom(zorder, smfa, a, include_atom))
  {
#ifdef DEBUG_BUILD_SMILES_ORDERING
    cerr << "Starting fragment  with atom " << a << " " << const_smarts_equivalent_for_atom(a) << endl;
#endif

    smi_info.add_start_atom(a);

    _build_smiles_ordering(identify_next_atom, INVALID_ATOM_NUMBER, a, zcounter, include_atom, smi_info);

    frag++;
  }

// If we built a subset, we may have fewer or more fragments detected

  if (0 == frag)
  {
    cerr << "Molecule::_build_smiles_ordering:no atoms selected!\n";
    return 0;
  }

  assert (NULL == include_atom ? frag == _fragment_information.number_fragments() : frag > 0);

#ifdef DEBUG_BUILD_SMILES_ORDERING
  cerr << "Smiles order array constructed, " << frag << " fragments\n";
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << "Atom " << i << " order is " << zorder[i];
    if (_fragment_information.number_fragments() > 1)
      cerr << " fragment " << _fragment_information.fragment_membership(i);
    cerr << endl;
  }

//if (_ring_closure_bonds.number_elements ())
//{
//  for (int i = 0; i < _ring_closure_bonds.number_elements (); i++) 
//  {
//    const Bond * b = _ring_closure_bonds[i];
//    cerr << "Ring closure bond " << i << " " << *b << endl;
//  }
//}
#endif

  return 1;
}

int
Molecule::_build_smiles_ordering (Smiles_Information & smi_info,
                                  const int * include_atom)
{
  if (0 == _number_elements)
    return 1;

  if (! _fragment_information.contains_valid_data())
    (void) number_fragments();

  smi_info.invalidate();

  smi_info.prepare_to_build_ordering(_number_elements);

  Smiles_First_Atom smfa;

  int rc = _build_smiles_ordering(smfa,
                              &Molecule::_smiles_choose_next_atom,
                              include_atom,
                              smi_info);

  if (0 == rc)
    return 0;

  if (NULL == include_atom)
    smi_info.set_smiles_order_type(DEFAULT_SMILES_ORDER_TYPE);
  else
    smi_info.set_smiles_order_type(SUBSET_SMILES_ORDER_TYPE);

  return rc;
}

/*
  Should an atom be included in a smiles
*/

int
Molecule::_include_atom_in_smiles (atom_number_t zatom) const
{
  const Atom * a = _things[zatom];

  const Element * e = a->element();

  if (1 != e->atomic_number())   // non H atoms always included
    return 1;

  if (0 != a->formal_charge())    // formally charged species always included
    return 1;

  if (a->isotope()) // isotopes always included
    return 1;

  return 0;
}

/*
  In a subset, is a given bond still in a ring
*/

int
Molecule::_ring_bond_in_subset (const int * include_atom,
                                atom_number_t a1,
                                atom_number_t a2)
{
  int nr = nrings();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = ringi(i);

    if (! ri->contains_bond(a1, a2))
      continue;

    if (ri->all_members_set_in_array(include_atom, 1))
      return 1;
  }

  return 0;
}

/*
  Identify the lowest order atom in each fragment
*/

/*void
Molecule::_find_smiles_start_atoms (const int * zorder, 
                                    resizable_array<int> & start_atom) const
{
  const int nf = _number_fragments;

  int * lowest_score = new_int (nf); iw_auto_array<int> free_lowest_score (lowest_score);

  start_atom.extend (nf, numeric_limits<int>::max());

  for (int i = 0; i < _number_elements; i++)
  {
    int j = _things[i]->fragment_membership ();
    if (zorder[i] < lowest_score[j])
    {
      lowest_score[j] = zorder[i];
      start_atom[j] = i;
    }
  }

  return;
}*/

/*
  When deciding from which atom the smiles should be constructed,
  we need to identify the unprocessed atom with the lowest ZORDER value
*/

/*static int
find_unprocessed_atom (const int * zorder, const int * already_done,
                       int n,
                       atom_number_t & result)
{
  int min_order = 9 * n;     // some large positive number
  int isave = -1;

  for (int i = 0; i < n; i++)
  {
    if (already_done[i])
      continue;

    if (zorder[i] < min_order)
    {
      min_order = zorder[i];
      isave = i;
    }
  }

  if (isave >= 0)
  {
    result = isave;
    return 1;
  }

  return 0;
}*/

/*
  The caller has identified a new atom to be added to the smiles. We
  must insert it in the arrays according to the order in zorder
*/

/*static void
insert_atom_and_bond (const int * zorder,
                      resizable_array<atom_number_t> & atoms,
                      atom_number_t a,
                      resizable_array<bond_type_t> &   bonds,
                      bond_type_t bt)
{
  int acon = atoms.number_elements ();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t b = atoms[i];
    if (zorder[a] < zorder[b])
    {
      atoms.insert_before (i, a);
      bonds.insert_before (i, bt);
      return;
    }
  }

// If we come out here, either the list is empty, or the new atom
// is ordered larger than all present. Just append

  atoms.add (a);
  bonds.add (bt);

  return;
}*/

/*
  We are processing the connections to ANCHOR, and need to build
  a sorted array of the adjacent bonds which need to be processed.
  We do an insertion sort, based on the ZORDER value of the atom
  at the other end of the bond
*/

static void
insert_bond (atom_number_t anchor,
             const int * zorder,
             resizable_array<const Bond *> & bonds,
             const Bond * b)
{
  int nb = bonds.number_elements ();
  atom_number_t zatom = b->other (anchor);

  for (int i = 0; i < nb; i++)
  {
    atom_number_t ai = bonds[i]->other (anchor);
    if (zorder[zatom] < zorder[ai])
    {
      bonds.insert_before (i, b);
      return;
    }
  }

// If we come out here, either the list is empty, or the new atom at
// the end of the bond is ordered higher than all others

  bonds.add(b);
  return;
}

//#define DEBUG_SMILES_FORMATION

/*
  Ran into problems with deciding which atoms should be in the smiles
  and which should be omitted. For now, if it is in the CT, it will
  be processed. 

  Very strange stuff when dealing with directional bonds.
  Consider

     2
      \
       1==3
      /    \
     0      4

  And consider starting with atom 2. The smiles is 2\1(\0)=3\4
  BUT, the 0-1 bond will be an / bond.
  Similarly, when starting with atom 0, the smiles is 0/1(/2)=3\4
  BUT, the 1-2 bond will be a \ bond
  When we detect that, we must react appropriately
*/

int
Molecule::_construct_smiles_for_fragment (Smiles_Formation_Info & sfi,
                                          Smiles_Information & smi_info)
{
  int * already_done = sfi.already_done();

  atom_number_t previous_atom = sfi.previous_atom();
  atom_number_t zatom = sfi.zatom();

#ifdef DEBUG_SMILES_FORMATION
  cerr << "_construct_smiles_for_fragment:continues with atom " << zatom << ", previous atom is " << previous_atom << endl;
#endif

  already_done[zatom] = 1;

  smi_info.add_atom(zatom);

  const Atom * a = _things[zatom];

  int acon = a->ncon();

  IWString & smiles = smi_info.smiles();

  const int * zorder = smi_info.smiles_order();

// If there are no other connections from here, no need to go looking
// for ring openings or prioritising connections

  if ((0 == acon) || (1 == acon && INVALID_ATOM_NUMBER != previous_atom))
  {
    (void) _process_atom_for_smiles(sfi, smiles);
    return 1;
  }

  resizable_array<const Bond *>  ring_opening_bonds;
  resizable_array<atom_number_t> ring_closures_found;

  resizable_array<const Bond *> process_these_bonds;

  const int * include_atom = sfi.include_atom();

// Build a sorted list of the connections to be processed, identifying
// any ring openings or closings.
// Note that we duplicate a lot of code here to avoid repeatedly testing INCLUDE_ATOM

  if (NULL == include_atom)
  {
    for (int i = 0; i < acon; i++)
    {
      const Bond * b = a->item(i);

      atom_number_t j = b->other(zatom);

      if (previous_atom == j)
        continue;

      if (already_done[j])              // closing a ring
        ring_closures_found.add(j);
      else if (smi_info.contains_ring_closure_bond(zatom, j))
        insert_bond (zatom, zorder, ring_opening_bonds, b);
      else                                      // just a regular connection
        insert_bond(zatom, zorder, process_these_bonds, b);
    }
  }
  else
  {
    for (int i = 0; i < acon; i++)
    {
      const Bond * b = a->item(i);

      atom_number_t j = b->other(zatom);

      if (previous_atom == j)
        continue;

      if (0 == include_atom[j])
        continue;

      if (already_done[j])              // closing a ring
        ring_closures_found.add(j);
      else if (smi_info.contains_ring_closure_bond(zatom, j))
        insert_bond(zatom, zorder, ring_opening_bonds, b);
      else                                      // just a regular connection
        insert_bond(zatom, zorder, process_these_bonds, b);
    }
  }

#ifdef DEBUG_INSERT_BOND
  int npb = process_these_bonds.number_elements();
  for (int i = 0; i < npb; i++)
  {
    const Bond * b = process_these_bonds[i];
    cerr << " i = " << i << " bond to " << b->other(zatom) << " order " << zorder[b->other(zatom)] << endl;
  }
#endif

  const Chiral_Centre * c = NULL;
  if (include_chiral_info_in_smiles())
    c = chiral_centre_at_atom(zatom);     // will be NULL if atom A is not a chiral centre

// In the case of ring closures to a chiral atom, we need to arrange the
// ring closures found in their ZORDER ordering.

  if (NULL == c)
    ;
  else if (2 == ring_closures_found.number_elements())
  {
    atom_number_t a0 = ring_closures_found[0];
    atom_number_t a1 = ring_closures_found[1];

    if (zorder[a0] > zorder[a1])
      ring_closures_found.swap_elements(0, 1);
  }
  else if (3 == ring_closures_found.number_elements())   // poor-man's sort
  {
    atom_number_t a0 = ring_closures_found[0];
    atom_number_t a1 = ring_closures_found[1];
    atom_number_t a2 = ring_closures_found[2];

    if (zorder[a0] > zorder[a1])
    {
      std::swap(a0, a1);
      ring_closures_found.swap_elements(0, 1);
    }
    if (zorder[a1] > zorder[a2])
    {
      std::swap(a1, a2);
      ring_closures_found.swap_elements(1, 2);
    }
    if (zorder[a0] > zorder[a1])
      ring_closures_found.swap_elements(0, 1);
  }

// Now that we have determined any ring openings, we can append the
// smiles symbol. We must wait until ring openings are determined for
// chiral atoms

  (void) _process_atom_for_smiles(sfi, zorder, ring_opening_bonds, ring_closures_found, c, smiles);

  Ring_Number_Manager & rnm = sfi.rnm();
  assert (rnm.ok());

// Our implementation of chiral atoms requires ring closures to be added
// before ring openings. BUT, we really don't want the same ring number
// appearing twice on one atom (once as a ring closure, and then as a
// new ring opening). So, we get the ring number manager to process the
// ring openings first, but put the results into a temporary

  IWString ring_opening_chars;
  int nro = ring_opening_bonds.number_elements();
  if (nro)
  {
    ring_opening_chars.resize(nro + 2);
    for (int i = 0; i < nro; i++)
    {
      rnm.store_ring(ring_opening_chars, ring_opening_bonds[i], zatom);
    }
  }

  assert (rnm.ok());

  if (ring_closures_found.number_elements())
  {
    rnm.append_ring_closures_for_atom(smiles, zatom, ring_closures_found, c);
  }

  assert (rnm.ok());

  if (nro)
    smiles += ring_opening_chars;

  assert (ring_opening_chars.ok());

#ifdef DEBUG_SMILES_FORMATION
  cerr << "After atom " << zatom << " smiles is now '" << smiles << "'\n";
  if (ring_closures_found.number_elements())
    cerr << "Found " << ring_closures_found.number_elements() << " ring closures\n";
  if (ring_opening_chars.number_elements())
  {
    cerr << "Ring openings to";
    for (int i = 0; i < ring_opening_chars.number_elements(); i++)
    {
      cerr << ' ' << ring_opening_chars[i] << ' ';
    }
    cerr << endl;
  }
#endif

// Handle the case of no further connections
// If there are no connections, there should not be any ring openings

  acon = process_these_bonds.number_elements();
  if (0 == acon)
  {
    if (ring_opening_bonds.number_elements())
      cerr << "No connections, but " << ring_opening_bonds.number_elements() << " ring openings\n";
    assert ( (NULL == include_atom) ? (0 == ring_opening_bonds.number_elements()) : 1);
    return 1;
  }

  int inc_ctb = include_cis_trans_in_smiles();    // do the call once for efficiency

  unsigned int inc_arom;    // aromatic bonds or not
  if (! sfi.write_smiles())
    inc_arom = 1;
  else
    inc_arom = get_include_aromaticity_in_smiles();

  if (write_single_bonds_in_smiles())
    inc_arom |= 2;

  int rc = 0;
  for (int i = 0; i < acon; i++)
  {
    const Bond * b = process_these_bonds[i];
//  cerr << "Bond " << i << " from atom " << zatom << " nrings = " << b->nrings() << " arom " << b->is_aromatic() << " inc_arom " << inc_arom << endl;

    if (i < acon - 1)
      smiles+= '(';

    atom_number_t j = b->other(zatom);

    if (inc_ctb && b->is_directional())
      _process_directional_bond_for_smiles(smiles, b, j);
    else if (_write_aromatic_bonds_as_colons && b->is_aromatic())
      smiles += ':';
    else
      b->append_bond_type(smiles, j, inc_arom);

    sfi.set_zatom(j);

    rc += _construct_smiles_for_fragment(sfi, smi_info);

    sfi.set_previous_atom(previous_atom);
    sfi.set_zatom(zatom);

    if (i < acon - 1)
      smiles += ')';
  }

  return rc;
}

/*
  A smiles for a fragment will begin with atom A
  We need to identify the bond down which the smiles will start
*/

const Bond *
Molecule::_identify_first_smiles_bond (atom_number_t zatom,
                                       const int * zorder)
{
  const Atom * a = _things[zatom];
  int acon = a->ncon();

  int smallest_zorder = -1;
  const Bond * rc = NULL;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);
    atom_number_t j = b->other(zatom);

    if (smallest_zorder < 0)
    {
      smallest_zorder = zorder[j];
      rc = b;
    }
    else if (zorder[j] < smallest_zorder)
    {
      smallest_zorder = zorder[j];
      rc = b;
    }
  }

  return rc;
}

int
Molecule::_construct_smiles (const Fragment_Information & frag_info,
                             Smiles_Information & smi_info,
                             const int * include_atom)
{
  int * already_done = new_int(_number_elements); iw_auto_array<int> free_already_done(already_done);

  smi_info.prepare_to_build_smiles(_number_elements);

  const int * zorder = smi_info.smiles_order();

#ifdef DEBUG_SMILES_FORMATION
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << "Atom " << i << " (" << atomic_symbol(i) << ") order is " << zorder[i];
    if (NULL != include_atom)
      cerr << " include_atom " << include_atom[i];
    cerr << endl;
  }
#endif

  const Set_of_Atoms & smiles_start_atom = smi_info.smiles_start_atom();

  int n = smiles_start_atom.number_elements();

  assert (n == frag_info.number_fragments());   // must be a smiles start atom in each fragment

  int rc = 0;

  int need_dot = 0;

#ifdef DEBUG_SMILES_FORMATION
  cerr << "Generating smiles for molecule with " << frag_info.number_fragments() << " fragments\n";
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << " atom " << i << " '" << smarts_equivalent_for_atom(i) << " in fragment " << frag_info.fragment_membership(i);
    if (NULL != include_atom)
      cerr << " include " << include_atom[i];
    cerr << endl;
  }
#endif

  IWString & smiles = smi_info.smiles();

  for (int i = 0; i < n; i++)
  {
    atom_number_t astart = smiles_start_atom[i];

    int f = frag_info.fragment_membership(astart);

    if (NULL == include_atom)    // no need to check anything
      ;
    else if (include_atom[astart])    // great, that atom is being processed
      ;
    else             // see if we can find something in fragment F
    {
      astart = _choose_highest_canonical_order_in_fragment(f, zorder, include_atom);
      if (INVALID_ATOM_NUMBER == astart)
        continue;
    }

    if (need_dot)
      smiles += '.';

    int nr = frag_info.rings_in_fragment(f);

#ifdef DEBUG_SMILES_FORMATION
    cerr << "Fragment " << f << " contains " << frag_info.bonds_in_fragment(f) << " bonds and " << frag_info.atoms_in_fragment(f) << " atoms, nr = " << nr << ", start atom " << astart << endl;
#endif

    Smiles_Formation_Info sfi(_number_elements, nr);
//  cerr << "Doing smarts? " << smi_info.smiles_is_smarts() << endl;

    if (NULL != smi_info.user_specified_atomic_smarts())
      sfi.set_user_specified_atomic_smarts(smi_info.user_specified_atomic_smarts());

    if (smi_info.smiles_is_smarts())
    {
      sfi.set_make_smarts_embedding(smi_info.create_smarts_embedding());
//    cerr << "Adding user specified atomic smarts " << smi_info.user_specified_atomic_smarts() << endl;
    }

    sfi.set_zatom(astart);
    sfi.set_already_done(already_done);
    sfi.set_include_atom(include_atom);

    if (write_smiles_with_smarts_atoms())
      sfi.set_write_smiles(0);

    rc += _construct_smiles_for_fragment(sfi, smi_info);

    need_dot = 1;
  }

  return rc;
}

/*
  The array _smiles_order can be ordered in a number of ways.
  _smiles_order_type must be updated to correspond with the
  ordering type
*/

const IWString &
Molecule::smiles ()
{
  assert (ok());

  if (! _smiles_information.contains_smiles())   // need to recompute
    ;
  else if (_smiles_information.smiles_is_smarts())
  {
    _smiles_information.make_empty();
    _smiles_information.set_smiles_is_smarts(0);
  }
  else
    return _smiles_information.smiles();

  if (0 == _number_elements)
  {
    _smiles_information.make_empty();
    return _smiles_information.smiles();
  }

  (void) number_fragments();

  if (! _smiles_information.contains_valid_ordering())
  {
    if (! _build_smiles_ordering(_smiles_information, NULL))
    {
      cerr << "Molecule::smiles: cannot construct ordering\n";
      _smiles_information.set_error();
      return _smiles_information.smiles();
    }

    _smiles_information.set_smiles_order_type(DEFAULT_SMILES_ORDER_TYPE);
  }

  _construct_smiles(_fragment_information, _smiles_information, NULL);

  return _smiles_information.smiles();
}

const IWString &
Molecule::smiles (Smiles_Information & smi_info,
                  const int * include_atom)
{
  if (0 == _number_elements)
  {
    smi_info.make_empty();
    return smi_info.smiles();
  }

  if (! _build_smiles_ordering(smi_info, include_atom))
  {
    cerr << "Molecule::smiles:cannot build subset smiles info\n";
    smi_info.set_error();
    return smi_info.smiles();
  }

  Fragment_Information frag_info;

  if (! compute_fragment_information(frag_info, include_atom))
  {
    cerr << "Molecule::smiles:cannot compute fragment info of subset\n";
    smi_info.set_error();
    return smi_info.smiles();
  }

  smi_info.set_smiles_order_type (SUBSET_SMILES_ORDER_TYPE);

  _construct_smiles(frag_info, smi_info, include_atom);

  return smi_info.smiles();
}

const IWString &
Molecule::random_smiles ()
{
  if (0 == _number_elements)
  {
    _smiles_information.make_empty();
    return _smiles_information.smiles();
  }

  invalidate_smiles();

  Smiles_First_Atom smfa;

  smfa.set_build_type(SMILES_FIRST_ATOM_RANDOM);

  (void) _build_smiles_ordering(smfa,
                                 &Molecule::_smiles_choose_random_next_atom,
                                 NULL,
                                 _smiles_information);

  _smiles_information.set_smiles_order_type(RANDOM_SMILES_ORDER_TYPE);

  _construct_smiles(_fragment_information, _smiles_information, NULL);

  return _smiles_information.smiles();
}

/*
  For helping people with text based programmes, we have the ability to
  start a smiles with any atom. This is a kludge, because this doesn't
  really fit well with how the smiles are built.
*/

const IWString &
Molecule::smiles_starting_with_atom (atom_number_t astart)
{
  return smiles_starting_with_atom(astart, _smiles_information, NULL);
}

const IWString &
Molecule::smiles_starting_with_atom (atom_number_t astart,
                                     Smiles_Information & smi_info,
                                     const int * include_atom)
{
  if (0 == _number_elements)
  {
    smi_info.make_empty();
    return smi_info.smiles();
  }

  assert (NULL == include_atom ? 1 : 0 != include_atom[astart]);

  smi_info.invalidate();

  smi_info.prepare_to_build_ordering(_number_elements);

  Smiles_First_Atom smfa;

  smfa.set_atom_number(astart);

  (void) _build_smiles_ordering(smfa,
                                 &Molecule::_smiles_choose_next_atom,
                                 include_atom,
                                 smi_info);

  Fragment_Information frag_info;
  if (! compute_fragment_information(frag_info, include_atom))
  {
    cerr << "Molecule::smiles_starting_with_atom:cannot find fragment info for subset\n";
    smi_info.set_error();
    return smi_info.smiles();
  }

  smi_info.set_smiles_order_type(RANDOM_SMILES_ORDER_TYPE);

  _construct_smiles(frag_info, smi_info, include_atom);

  return smi_info.smiles();
}

/*
  Someone may need to know the order of the atoms in the smiles
*/

int
Molecule::smiles_atom_order (int * s)
{
  (void) smiles();     // will force construction of the array(s)

  const int * so = _smiles_information.smiles_order();

  copy_vector(s, so, _number_elements);

  return 1;
}

/*
  Don't use nice_smiles, I don't think it works properly - if ring membership
  hasn't been computed, it will call _build_smiles_ordering, and that seems
  sub-optimal.

  I've never needed nice_smiles(), so maybe it should be removed.
*/

const IWString &
Molecule::nice_smiles ()
{
  if (0 == _number_elements)
  {
    _smiles_information.make_empty();
    return _smiles_information.smiles();
  }

  _smiles_information.set_smiles_is_smarts(0);

  if (NICE_SMILES_ORDER_TYPE == _smiles_information.smiles_order_type())
    return _smiles_information.smiles();

  return nice_smiles(_smiles_information, NULL);
}

const IWString &
Molecule::nice_smiles (Smiles_Information & smi_info,
                       const int * include_atom)
{
  if (0 == _number_elements)
  {
    _smiles_information.make_empty();
    return _smiles_information.smiles();
  }

  if (NULL != include_atom || NICE_SMILES_ORDER_TYPE != smi_info.smiles_order_type())
  {
    smi_info.prepare_to_build_ordering(_number_elements);

    Smiles_First_Atom smfa;

    smfa.set_build_type(SMILES_FIRST_ATOM_NICE);

    (void) _build_smiles_ordering( smfa,
                                   &Molecule::_smiles_choose_nice_next_atom,
                                   include_atom,
                                   smi_info);
  }

  if (NULL == include_atom)
    smi_info.set_smiles_order_type(NICE_SMILES_ORDER_TYPE);
  else
    smi_info.set_smiles_order_type(SUBSET_SMILES_ORDER_TYPE);

  Fragment_Information frag_info;

  if (! compute_fragment_information(frag_info, include_atom))
  {
    cerr << "Molecule::nice_smiles:cannot compute fragment info for subset\n";
    smi_info.set_error();
    return smi_info.smiles();
  }

  _construct_smiles(frag_info, smi_info, include_atom);

  return smi_info.smiles();
}

//#define DEBUG_UNIQUE_SMILES

const IWString &
Molecule::_unique_smiles (const Fragment_Information & frag_info,
                          Smiles_Information & smi_info,
                          Symmetry_Class_and_Canonical_Rank & sccr,
                          const int * include_atom)
{
//cerr << "Allocated? " << sccr.arrays_allocated() << endl;

  if (0 == _number_elements)
  {
    smi_info.make_empty();
    return smi_info.smiles();
  }

  smi_info.prepare_to_build_ordering(_number_elements);

  compute_aromaticity_if_needed();

#ifdef DEBUG_UNIQUE_SMILES
  cerr << "Computing unique smiles for " << _number_elements << " atoms\n";
  if (! sccr.arrays_allocated())
    cerr << "Computing canonical rank\n";
  else
    cerr << "Canonical rank already computed\n";
#endif

  if (! sccr.arrays_allocated())
    compute_canonical_ranking(sccr, include_atom);

#ifdef DEBUG_UNIQUE_SMILES
  cerr << "Canonical rank computed\n";
  const int * c = sccr.canonical_rank();
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << "Atom " << i << " type " << _things[i]->atomic_symbol();
    if (NULL != include_atom)
      cerr << " include " << include_atom[i];
    cerr << " rank " << c[i] << endl;
  }

  cerr << "Order type " << smi_info.smiles_order_type() << '\n';
#endif

  if (include_cis_trans_in_smiles())
    _adjust_cis_trans_bonds_to_canonical_form(sccr.canonical_rank());

  assert (NULL != _aromaticity);    // aromaticity computed in compute_canonical_ranking

  Smiles_First_Atom smfa;

  smfa.set_build_type(SMILES_FIRST_ATOM_UNIQUE);

  (void) _build_smiles_ordering(smfa,
                                 &Molecule::_smiles_choose_unique_next_atom,
                                 include_atom,
                                 smi_info);

  if (NULL == include_atom)
    smi_info.set_smiles_order_type(UNIQUE_SMILES_ORDER_TYPE);
  else
    smi_info.set_smiles_order_type(SUBSET_SMILES_ORDER_TYPE);

#ifdef DEBUG_UNIQUE_SMILES
  cerr << "Smiles unique order is\n";
  const int * s = smi_info.smiles_order();

  for (int i = 0; i < _number_elements; i++)
  {
    cerr << " i = " << i << " order = " << s[i] << endl;
  }
#endif

  _construct_smiles(frag_info, smi_info, include_atom);

  return smi_info.smiles();
}

/*
  Be careful using unique_smiles() and non_aromatic_unique_smiles()
  If you call unique_smiles() and then non_aromatic_unique_smiles()
  you will likely get the same smiles back. To force a recomputation,
  call invalidate_smiles() in between...

  Jul 2000. Introduce a default aromaticity that is used for all
  unique smiles
*/

static int default_unique_smiles_aromaticity = Daylight;

int
set_default_unique_smiles_aromaticity (int a)
{
  default_unique_smiles_aromaticity = a;

  return 1;
}

/*
  We have some common tasks that need to happen when doing unique smiles determinations.
  Some things need to be set and then restored.
*/

class Hold_and_Restore_Global_Settings
{
  private:
    int _aromsave;     // save the global aromaticity definition

    int _incaromsave;  // save the include aromaticity in smiles definition

  public:
    Hold_and_Restore_Global_Settings (int incarom);
    ~Hold_and_Restore_Global_Settings ();

    int aromaticity_changed () const { return _aromsave != default_unique_smiles_aromaticity;}
};

Hold_and_Restore_Global_Settings::Hold_and_Restore_Global_Settings (int incarom)
{
  _aromsave = global_aromaticity_type();

  set_global_aromaticity_type(default_unique_smiles_aromaticity);

  _incaromsave = get_include_aromaticity_in_smiles();

  set_include_aromaticity_in_smiles(incarom);

  return;
}

Hold_and_Restore_Global_Settings::~Hold_and_Restore_Global_Settings ()
{
  set_global_aromaticity_type(_aromsave);

  set_include_aromaticity_in_smiles(_incaromsave);

  return;
}

const IWString &
Molecule::unique_smiles ()
{
  if (UNIQUE_SMILES_ORDER_TYPE == _smiles_information.smiles_order_type())
    return _smiles_information.smiles();

  Hold_and_Restore_Global_Settings hrgs(1);

  if (hrgs.aromaticity_changed())    // we may have Pearlman aromaticity, but need Daylight for unique smiles
    compute_aromaticity();
  else 
    compute_aromaticity_if_needed();

  _smiles_information.set_smiles_is_smarts(0);

  return _unique_smiles(_fragment_information, _smiles_information, _symmetry_class_and_canonical_rank, NULL);
}

const IWString &
Molecule::unique_smiles (Smiles_Information & smi_info,
                         const int * include_atom)
{
  assert (NULL != include_atom);

  if (0 == _number_elements)
  {
    _smiles_information.make_empty();
    return _smiles_information.smiles();
  }

  _smiles_information.set_smiles_is_smarts(0);

//_fragment_information.debug_print(cerr);

  Hold_and_Restore_Global_Settings hrgs(1);

  if (hrgs.aromaticity_changed())
    compute_aromaticity();

  Fragment_Information frag_info;

  if (! compute_fragment_information(frag_info, include_atom))
  {
    cerr << "Molecule::unique_smiles:cannot compute fragment info for subset\n";
    smi_info.set_error();
    return smi_info.smiles();
  }

  Symmetry_Class_and_Canonical_Rank sccr;

  if (! sccr.allocate_arrays(_number_elements))
    return smi_info.set_error();

  compute_canonical_ranking(sccr, include_atom);

  if (include_cis_trans_in_smiles())
    _adjust_cis_trans_bonds_to_canonical_form(sccr.canonical_rank());

// _smiles_choose_unique_*_atom need to have the molecule's canonical order fixed. Make
// a copy of any existing data in _symmetry_class_and_canonical_rank and store the
// values from sccr into _symmetry_class_and_canonical_rank

  Symmetry_Class_and_Canonical_Rank sccr_save;
  sccr_save.store_values_from(_symmetry_class_and_canonical_rank, _number_elements);

  _symmetry_class_and_canonical_rank.store_values_from(sccr, _number_elements);

  _unique_smiles(frag_info, smi_info, sccr, include_atom);

  _symmetry_class_and_canonical_rank.store_values_from(sccr_save, _number_elements);

  return smi_info.smiles();
}

const IWString &
Molecule::non_aromatic_unique_smiles ()
{
  if (0 == _number_elements)
  {
    _smiles_information.make_empty();
    return _smiles_information.smiles();
  }

  _smiles_information.set_smiles_is_smarts(0);

  Hold_and_Restore_Global_Settings hrgs(0);

  if (hrgs.aromaticity_changed())
    compute_aromaticity();

  return _unique_smiles(_fragment_information, _smiles_information, _symmetry_class_and_canonical_rank, NULL);
}

/*
  In processing a fused system, it has been discerned that some rings
  which had previously been assigned separate fused system identifiers
  are in fact part of the same fused system.
  We examine rings in RINGS, starting with ring RSTART.
  The common fused system identifier is FUSED_SYSTEM_IDENTIFIER.
  The fused system identifiers which need to be changed are in FUSED_SYS_IDS_TO_BE_CHANGES
*/

int
Molecule::_merge_fused_system_identifiers (resizable_array<Ring *> & rings,
                                           int rstart,
                                           int fused_system_identifier,
                                           resizable_array<int> & fused_sys_ids_to_be_changed)
{
  int rc = 0;
  int nr = rings.number_elements();
  for (int i = rstart; i < nr; i++)
  {
    Ring * r = rings[i];
    if (! r->is_fused())
      continue;

    int rfsysid = r->fused_system_identifier();
    if (rfsysid == fused_system_identifier)
      continue;

    if (fused_sys_ids_to_be_changed.contains(rfsysid))
    {
      r->set_fused_system_identifier(fused_system_identifier);
      rc++;
    }
  }

  return rc;
}

int
Molecule::_find_raw_rings (const atom_number_t previous_atom,
                           const atom_number_t current_atom,
                           resizable_array<Ring *> & rings,
                           resizable_array<atom_number_t> & active_rings,
                           int * already_done)
{
  assert (0 == already_done[current_atom]);
  assert (rings.number_elements() == active_rings.number_elements());

  already_done[current_atom] = 1;
//cerr << "GBORETN: processing atom " << current_atom << endl;

  const Atom * c = _things[current_atom];

  int acon = c->ncon();
  if (0 == acon)
    return rings.number_elements();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = c->other(current_atom, i);
    if (previous_atom == j)
      continue;

    if (already_done[j])
    {
//    cerr << "Found new ring to atom " << j << endl;
      Ring * tmp = new Ring;
      tmp->resize(8);
      tmp->add(j);
      tmp->add(current_atom);
      rings.add(tmp);
      tmp->set_fragment_membership(_fragment_information.fragment_membership(current_atom));

      active_rings.add(j);
    }
  }

// Recursively call this function for each bond attached to CURRENT_ATOM

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = c->other(current_atom, i);
    if (already_done[j])
      continue;

    int rstart = rings.number_elements();
    (void) _find_raw_rings(current_atom, j, rings, active_rings, already_done);

    int nrings_now = rings.number_elements();
    if (rstart == nrings_now)   // no new rings found down this bond
      continue;

//  Count the number of new rings down this bond, and try to determine
//  any existing fusion specifications. Note that we may find different
//  fused sys identifiers in the list of new rings. These numbers must
//  all be consolidated, as they now are known to belong to the same
//  fused system,

    int number_new_rings = 0;
    int fused_system_identifier = -99;

//  cerr << "Current = " << current_atom << " to " << j << 
//          " now " << rings.number_elements() << " rings\n";

    resizable_array<int> fused_sys_ids_to_be_changed;

    for (int k = rstart; k < nrings_now; k++)
    {
      if (INVALID_ATOM_NUMBER == active_rings[k])    // ring has been processed
        continue;

      number_new_rings++;
      Ring * r = rings[k];

      if (! r->is_fused())    // not interested in isolated rings
        continue;

      int rsysid = r->fused_system_identifier();    

      if (fused_system_identifier < 0)        // first fused ring found here
        fused_system_identifier = rsysid;
      else if (rsysid == fused_system_identifier)    // already has same id, no change
        ;
      else     // different systems need to be merged
        fused_sys_ids_to_be_changed.add(rsysid);

//    if (r->is_fused())
//      cerr << "FBLOGD current = " << current_atom << " con = " << i << " atom " << j << 
//              " k = " << k << " found fused sys identifier " << r->fused_system_identifier() << 
//              " to " << active_rings[k] << endl;
    }

    if (0 == number_new_rings)
      continue;

    if (-99 == fused_system_identifier)
      fused_system_identifier = current_atom;
    else
      _merge_fused_system_identifiers(rings, rstart, fused_system_identifier, fused_sys_ids_to_be_changed);

//  cerr << "Will assign fused system identifier " << fused_system_identifier << endl;

    for (int k = rstart; k < nrings_now; k++)
    {
      if (INVALID_ATOM_NUMBER == active_rings[k])
        continue;

      Ring * r = rings[k];
      if (number_new_rings > 1)
        r->set_fused_system_identifier(fused_system_identifier);

      if (active_rings[k] == current_atom)    // ring is complete, terminate it
        active_rings[k] = INVALID_ATOM_NUMBER;
      else
        r->add(current_atom);
    }
  }

  assert (rings.number_elements() == active_rings.number_elements());

  return rings.number_elements();
}

//#define DEBUG_FIND_RAW_RINGS_FOR_FRAGMENT

/*
  Process the raw rings for atoms with _fragment_membership == id

  We need to be careful with systems with spiro fusions to ring systems,

C1CCCC2C1CC(C2)1OC(OC1)1OCCO1 3 spiro

  for example. Since the spiro rings are not fused, we can correclty set
  their final ring membership.
  By default, we set the ring membership of fused systems to IW_RING_MEMBERSHIP_IS_A_RING_ATOM,
  but that would overwrite the values correctly found for the spiro rings. So, if we have
  the case of spiro rings joined to a ring system, force a complete SSSR
*/

int
Molecule::_find_raw_rings_for_fragment (int id, int * already_done)
{
  if (0 == nrings())
    return 1;

  if (id < 0 || id >= _fragment_information.number_fragments())
  {
    cerr << "Molecule::_find_raw_rings_for_fragment:finding rings in fragment " << id << " but only " << _fragment_information.number_fragments() << " fragments\n";
    debug_print(cerr);
    assert (NULL == "This is very bad");
  }

  if (NULL == _ring_membership)
    _initialise_ring_membership();

// Initialise all these atoms as 0 ring membership

  int atoms_being_processed = 0;

  atom_number_t start_atom = INVALID_ATOM_NUMBER;   // first atom in the fragment

  const int * fragment_membership = _fragment_information.fragment_membership();

  for (int i = 0; i < _number_elements; i++)
  {
    if (id == fragment_membership[i])
    {
      _ring_membership[i] = 0;
      atoms_being_processed++;
      if (INVALID_ATOM_NUMBER == start_atom)
        start_atom = i;
    }
  }

  int nr = _fragment_information.rings_in_fragment(id);

  assert (nr >= 0);

  if (0 == nr)      // no rings in this fragment
    return 1;

  assert (atoms_being_processed > 2);

  resizable_array<Ring *> rings;

  resizable_array<atom_number_t> active_rings;
  active_rings.resize(nrings());

  _find_raw_rings(INVALID_ATOM_NUMBER, start_atom, rings, active_rings, already_done);

  nr = rings.number_elements();

  assert (nr > 0);

#ifdef DEBUG_FIND_RAW_RINGS_FOR_FRAGMENT
  cerr << "Found " << nr << " rings for fragment " << id << endl;
  for (int i = 0; i < nr; i++)
  {
    cerr << "Ring " << i << " ";
    const Ring * ri = rings[i];
    cerr << (*ri) << endl;
    if (! ok_ring(ri))
    {
      cerr << "Very bad news, not a valid ring\n";
      iwabort();
    }
  }
#endif

// Update ring membership with the details

  int number_fused_rings = 0;

  for (int i = 0; i < nr; i++)
  {
    Ring * ri = rings[i];

    if (ri->is_fused())
      number_fused_rings++;
    else
    {
      ri->increment_vector(_ring_membership, 1);
      _add_ring_to_sssr(ri);
    }
  }

  if (0 == number_fused_rings)
    return nr;

  assert (nr > 1);     // can't be just one fused ring

// Accumulate the fused_system_identifiers of the fused systems that need to be processed with their attached spiro rings

  resizable_array<int> spiro_between_isolated_and_fused;

  if (number_fused_rings < nr)    // must be both fused and isolated rings present - check for spiro fusions between them
  {
    for (int i = 0; i < nr; i++)
    {
      const Ring * ri = rings[i];

      if (! ri->is_fused())
        continue;

      if (ri->fused_ring_check_for_spiro_fusion(_ring_membership))
        spiro_between_isolated_and_fused.add_if_not_already_present(ri->fused_system_identifier());
    }
  }

#ifdef DEBUG_FIND_RAW_RINGS_FOR_FRAGMENT
  cerr << "After examining rings, spiro between isolated and fused = " << spiro_between_isolated_and_fused.number_elements() << ", nr = " << nr << endl;
#endif

// Fused rings that are not bonded to a spiro ring to go raw rings

  for (int i = 0; i < nr; i++)
  {
    Ring * ri = rings[i];

    if (! ri->is_fused())
      continue;

    int fsid = ri->fused_system_identifier();
    if (spiro_between_isolated_and_fused.contains(fsid))
      continue;

    ri->set_vector(_ring_membership, IW_RING_MEMBERSHIP_IS_A_RING_ATOM);

    _raw_rings.add(ri);
//  _experimental_raw_rings.add(ri);
  }

  if (spiro_between_isolated_and_fused.number_elements())
  {
//  cerr << "Spiro fusion to fused system\n";
    for (int i = 0; i < spiro_between_isolated_and_fused.number_elements(); i++)
    {
      int fsid = spiro_between_isolated_and_fused[i];

      _handle_spiro_between_isolated_and_fused(rings, fsid, already_done);
    }

    for (int i = 0; i < nr; i++)
    {
      Ring * ri = rings[i];

      if (! ri->is_fused())
        continue;

      int fsid = ri->fused_system_identifier();
      if (spiro_between_isolated_and_fused.contains(fsid))
        delete ri;
    }
  }

  return nr;
}

/*
  A new isolated ring has been found.  Update the _ring_membership,
  but only values greater than 0 - the others are undetermined yet.

  Sept 97. Previously this set ring membership to 1 only if the
  existing value of _ring_membership was zero. This failed in
  the case of spiro fused rings, so now we check for >= 0 values
  and always increment.

  also tell the bonds about the new ring
*/

int
Molecule::_update_ring_membership (const Ring * r)
{
  int ring_size = r->number_elements();
  atom_number_t prev_atom = r->last_item();

  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = r->item(i);

    cerr << "Molecule::_update_ring_membership: atom " << j << " rm = " << _ring_membership[j] << endl;

    if (IW_RING_MEMBERSHIP_IS_A_RING_ATOM  == _ring_membership[j])   // probably part of a fused system - spiro fusion
      _ring_membership[j] = 1;
    else if (_ring_membership[j] >= 0)
      _ring_membership[j]++;

    cerr << "Molecule::_update_ring_membership: atom " << j << " rm = " << _ring_membership[j] << endl;
    const Atom * aj = _things[j];

    Bond * b = const_cast<Bond *>(aj->bond_to_atom(prev_atom));   // loss of const OK
    assert (b);

    b->set_nrings(1);

    prev_atom = j;
  }

  return 1;
}

int
Molecule::_find_raw_rings (int * already_done)
{
  int nf = _fragment_information.number_fragments();

//cerr << "Finding raw rings for " << nf << " fragments\n";

  for (int i = 0; i < nf; i++)
  {
//  cerr << "Finding raw rings for fragment " << i << endl;
    if (! _find_raw_rings_for_fragment(i, already_done))
    {
      cerr << "Molecule::_find_raw_rings: Bad news, cannot get raw rings for fragment " << i << endl;
      debug_print(cerr);
      iwabort();
    }
  }

  return _raw_rings.number_elements();
}

int
Molecule::_find_raw_rings ()
{
  assert (0 == _raw_rings.number_elements());

  (void) number_fragments();

  if (0 == nrings())
    return 1;

  int * tmp = new_int(_number_elements); iw_auto_array<int> free_tmp(tmp);

  return _find_raw_rings(tmp);
}

int
smiles_error_message (const char * smiles,
                      int length_of_smiles,
                      int characters_processed,
                      const char * message)
{
  if (! display_smiles_interpretation_error_messages())
    return 1;

  assert (message);

  cerr << message << endl;

  int smiles_chars_to_print = characters_processed + 10;
  if (smiles_chars_to_print > length_of_smiles || smiles_chars_to_print > 80)
    smiles_chars_to_print = length_of_smiles;

  for (int i = 0; i < smiles_chars_to_print; i++)
  {
    cerr << smiles[i];
  }
  cerr << endl;

//cerr << "                     ";
  for (int i = 0; i < characters_processed; i++)
    cerr << ' ';
  cerr << "^\n";

  return 1;
}

void
Molecule::_add_ring_to_sssr (Ring * ri)
{
  ri->set_ring_number(_sssr_rings.number_elements());
  _sssr_rings.add(ri);

// _experimental_sssr_rings.add(ri);

// Update the bond ring membership

  atom_number_t aprev = ri->last_item();

  int ring_size = ri->number_elements();

  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = ri->item(i);

    _things[j]->in_another_ring();

    Bond * b = const_cast<Bond *>(_things[j]->bond_to_atom(aprev));

    b->in_another_ring();

    aprev = j;
  }

  return;
}

//#define DEBUG_HANDLE_SPIRO_BETWEEN_ISOLATED_AND_FUSED

/*
  Spiro fusions between isolated and fused systems are problematic. We force SSSR
  determination of all the fused rings in the fragment
*/

int
Molecule::_handle_spiro_between_isolated_and_fused (const resizable_array<Ring *> & rings,
                                                    int fsid,
                                                    int * atmp)
{
  set_vector(atmp, _number_elements, 0);

  int nr = rings.number_elements();

#ifdef DEBUG_HANDLE_SPIRO_BETWEEN_ISOLATED_AND_FUSED
  cerr << "Handling " << nr << " rings for possible spiro/fused systems, _sssr_rings.number_elements  = " << _sssr_rings.number_elements() << endl;
#endif

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = rings[i];

    if (! ri->is_fused())
      continue;

    if (fsid != ri->fused_system_identifier())
      continue;

    ri->set_vector(atmp, 1);
  }

  return _pearlman_sssr(atmp, 1);
}

int
Molecule::_all_atoms_are_chain_atoms (const int * process_these_atoms)
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == process_these_atoms[i])
      continue;

    if (is_ring_atom(i))
      return 0;
  }

  return 1;
}

/*int
Molecule::_determine_ring_closure_bonds (const int * zorder,
                                         const int * include_atom)
{
  assert (NULL != include_atom);

  int * already_done = new_int (_number_elements); iw_auto_array<int> free_already_done (already_done);

  int n = _smiles_start_atom.number_elements ();

  for (int i = 0; i < n; i++)
  {
    atom_number_t astart = _smiles_start_atom[i];

    if (NULL == include_atom)
      ;
    else if (include_atom[astart])
      ;
    else
    {
      astart = _choose_highest_canonical_order_in_fragment (i, zorder, include_atom);
      if (INVALID_ATOM_NUMBER == astart)
        continue;
    }

    _determine_ring_closure_bonds (INVALID_ATOM_NUMBER, astart, zorder, include_atom, already_done);
  }

  return 1;
}*/

/*
  We need to re-determine the ring closure bonds when dealing with a subset
  The canonical order is already known
*/

/*int
Molecule::_determine_ring_closure_bonds (atom_number_t aprev,
                                         atom_number_t zatom,
                                         const int * zorder,
                                         const int * include_atom,
                                         int * already_done)
{
  assert (NULL != include_atom);    // subset only

  already_done[zatom] = 1;

  const Atom * a = _things[zatom];

  int acon = a->ncon ();

  resizable_array<const Bond *> process_these_bonds;

  int rc = 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item (i);

    atom_number_t j = b->other (zatom);

    if (aprev == j)
      continue;

    if (0 == include_atom[j])
      continue;

    if (already_done[j])
    {
      _ring_closure_bonds.add (zatom, j);
      rc++;
    }
    else
    {
      insert_bond (zatom, zorder, process_these_bonds, b);
    }
  }

  int n = process_these_bonds.number_elements ();

  for (int i = 0; i < n; i++)
  {
    const Bond * b = process_these_bonds[i];

    atom_number_t j = b->other (zatom);

    if (already_done[j])
      continue;

    rc += _determine_ring_closure_bonds (zatom, j, zorder, include_atom, already_done);
  }

  return rc;
}*/

/*
  We are doing the smiles for a fragment, but there is a subset, and the molecule's default
  _smiles_start_atom was not part of the subset. Find another suitable starting point
*/

atom_number_t
Molecule::_choose_highest_canonical_order_in_fragment (int f,
                                                       const int * zorder,
                                                       const int * include_atom) const
{
  atom_number_t rc = INVALID_ATOM_NUMBER;
  int z;

  const int * fragment_membership = _fragment_information.fragment_membership ();

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == include_atom[i])
      continue;

    if (f != fragment_membership[i])
      continue;

    if (INVALID_ATOM_NUMBER == rc || z > zorder[i])
    {
      rc = i;
      z = zorder[i];
    }
  }

  return rc;
}

int
Smiles_Information::allocate_user_specified_atomic_smarts()
{
  assert (NULL == _user_specified_atomic_smarts);

  assert (_natoms > 0);

  _user_specified_atomic_smarts = new IWString[_natoms];

  assert (NULL != _user_specified_atomic_smarts);

  return 1;
}

IWString &
Smiles_Information::user_specified_atomic_smarts(atom_number_t zatom)
{
  return _user_specified_atomic_smarts[zatom];
}

const IWString &
Smiles_Information::user_specified_atomic_smarts(atom_number_t zatom) const
{
  return _user_specified_atomic_smarts[zatom];
}

void
Smiles_Information::set_user_specified_atomic_smarts(atom_number_t zatom,
                                                const IWString & s)
{
  if (NULL != _user_specified_atomic_smarts)
    ;
  else if (_natoms <= 0)
  {
    cerr << "Smiles_Information::set_user_specified_atomic_smarts:atom count unknown\n";
    return;
  }
  else
    _user_specified_atomic_smarts = new IWString[_natoms];

  _user_specified_atomic_smarts[zatom] = s;

  return;
}

IWString
Molecule::isotopically_labelled_smiles()
{
  if (0 == _number_elements)
    return (".");

  int * isave = NULL;

  for (int i = 0; i < _number_elements; i++)
  {
    int iso = _things[i]->isotope();

    if (0 == iso)
    {
      _things[i]->set_isotope(i);
      continue;
    }

    if (NULL == isave)
      isave = new_int(_number_elements);

    isave[i] = iso;

    _things[i]->set_isotope(i);
  }

  Smiles_Information sminfo;

  (void) number_fragments();

  if (! _build_smiles_ordering(sminfo, NULL))
  {
    cerr << "Molecule::isotopically_labelled_smiles: cannot construct ordering\n";
    _smiles_information.set_error();
    return sminfo.smiles();
  }

  sminfo.set_smiles_order_type(DEFAULT_SMILES_ORDER_TYPE);

  _construct_smiles(_fragment_information, sminfo, NULL);

  if (NULL != isave)
  {
    for (int i = 0; i < _number_elements; i++)
    {
      _things[i]->set_isotope(isave[i]);
    }

    delete [] isave;
  }
  else
  {
    for (int i = 0; i < _number_elements; i++)
    {
      _things[i]->set_isotope(0);
    }
  }

  return (sminfo.smiles());
}
