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
#include <iomanip>

#include "molecule.h"
#include "smiles.h"
#include "misc2.h"

/*
  If a connection is constructed without any info, all info is invalid
*/

Connection::Connection ()
{
  _a2 = INVALID_ATOM_NUMBER;
  _btype = NOT_A_BOND;

  return;
}

Connection::Connection (atom_number_t a2, bond_type_t btype)
{
  assert (a2 >= 0);
  assert (OK_BOND_TYPE (btype));

  _a2 = a2;
  _btype = btype;
  return;
}

Connection::~Connection ()
{
  _a2 = INVALID_ATOM_NUMBER;
  _btype = NOT_A_BOND;

  return;
}

void
Bond::adjust_for_loss_of_atom (atom_number_t i)
{
  assert (ok ());

  if (_a1 > i)
    _a1--;
  else if (_a1 == i)
  {
    cerr << "Bond::adjust_for_loss_of_atom: involves atom " << i << endl;
    debug_print (cerr);
    assert (NULL == "This should not happen");
  }

  if (_a2 > i)
    _a2--;
  else if (_a2 == i)
  {
    cerr << "Bond::adjust_for_loss_of_atom: involves atom " << i << endl;
    debug_print (cerr);
    assert (NULL == "This should not happen");
  }

  return;
}

void
Connection::set_bond_type (bond_type_t bt)
{
  _btype = bt;

  return;
}

void
Connection::set_aromatic ()
{
  SET_AROMATIC_BOND (_btype);
}

void
Connection::set_non_aromatic ()
{
  SET_NON_AROMATIC_BOND (_btype);
}

void
Connection::set_permanent_aromatic (int s)
{
  if (s)
    _btype = (_btype | PERMANENT_AROMATIC_BOND);
  else if (_btype & PERMANENT_AROMATIC_BOND)
    _btype = (_btype ^ PERMANENT_AROMATIC_BOND);

  return;
}

void
Bond::_default_values ()
{
#ifdef BONDS_KNOW_RING_MEMBERSHIP
//_nrings =  UNKNOWN_BOND_NRINGS;
  _nrings = 0;     // too hard otherwise, but this is dangerous
#endif

  _directional = 0;

  _bond_number = -1;
}

/*
  There are two constructors for a bond. One done in the presence of a
  molecule, and one done without a molecule.
  The constructor in the presence of a molecule audits its arguments much
  more than the one without.
*/

Bond::Bond (atom_number_t a1, atom_number_t a2, bond_type_t btype)
        : Connection (a2, btype)
{
  assert (a1 >= 0 && a2 >= 0 && a1 != a2);
  assert (OK_BOND_TYPE (btype));

  _default_values ();

  _a1 = a1;

  return;
}

Bond::Bond (const Bond & rhs)
{
  _a2 = rhs._a2;
  _btype = rhs._btype;
  _a1 = rhs._a1;
  _directional = rhs._directional;

#ifdef BONDS_KNOW_RING_MEMBERSHIP
  _nrings = rhs._nrings;
#endif

  return;
}

Bond::Bond (const Molecule *m, atom_number_t a1, atom_number_t a2, bond_type_t btype)
       : Connection (a2, btype)
{
  assert (OK_2_ATOMS (m, a1, a2));
  assert (OK_BOND_TYPE (btype));

  _default_values ();

  _a1 = a1;

  return;
}

Bond::Bond()
{
  _default_values();

  _a1 = INVALID_ATOM_NUMBER;
  _a2 = INVALID_ATOM_NUMBER;
  _btype = INVALID_BOND_TYPE;

  return;
}

Bond::~Bond ()
{
  assert (ok ());

  _a1 = INVALID_ATOM_NUMBER;
}

int
Bond::ok () const
{
//return 1;    // remove checking for VDOM

  if (INVALID_ATOM_NUMBER == _a1 || 
      INVALID_ATOM_NUMBER == _a2 ||
      ! OK_BOND_TYPE (_btype) )
    return 0;

  if (_a1 == _a2)
    return 0;

  return 1;
}

/*
  Does this bond consist of atoms Z1 and Z2 (in any order)
*/

int
Bond::involves (atom_number_t z1, atom_number_t z2) const
{
  assert (ok ());

  if (z1 == _a1)
    return z2 == _a2;
  if (z1 == _a2)
    return z2 == _a1;

  return 0;
}

int
Bond::involves_and_what_is_other (atom_number_t z1, atom_number_t & z2) const
{
  assert (ok ());

  if (_a1 == z1)
  {
    z2 = _a2;
    return 1;
  }
  else if (_a2 == z1)
  {
    z2 = _a1;
    return 1;
  }
  else
    return 0;
}

int
Bond::joins (const Bond * b, atom_number_t & common_atom) const
{
  if (_a1 == b->_a1 || _a1 == b->_a2)
  {
    common_atom = _a1;
    return 1;
  }

  if (_a2 == b->_a1 || _a2 == b->_a2)
  {
    common_atom = _a2;
    return 1;
  }

  return 0;     // no atom in common
}

void
Bond::set_bond_type (bond_type_t bt)
{
  assert (OK_BOND_TYPE (bt));

  Connection::set_bond_type (bt);

  return;
}

int
Bond::debug_print (ostream & os) const
{
  assert (os.good ());

  os << "Bond between " << _a1 << " and " << _a2 << " with bond type ";
  if (IS_SINGLE_BOND (_btype))
    os << '1';
  else if (IS_DOUBLE_BOND (_btype))
    os << '2';
  else if (IS_TRIPLE_BOND (_btype))
    os << '3';
  else if (IS_PERMANENT_AROMATIC_BOND (_btype))
    os << ":";
  else if (AROMATIC_BOND == _btype)
    os << ':';
  else if (0 == _btype)
    os << '0';    /// !!!!!
  else
  {
    cerr << "What type " << _btype << endl;
    os << '?';
  }

  if (! OK_BOND_TYPE (_btype))
    os << " BAD BOND " << hex << _btype << dec;

#ifdef BONDS_KNOW_RING_MEMBERSHIP
  if (UNKNOWN_BOND_NRINGS != _nrings)
    os << " in " << _nrings << " rings";
#endif

  if (IS_AROMATIC_BOND (_btype))
    os << " (aromatic)";

  if (0 == _directional)
    ;
  else if (is_directional_up ())
    os << " up";
  else if (is_directional_down ())
    os << " dn";
  else if (part_of_cis_trans_grouping ())
    os << " ct";
  else if (is_wedge_up())
    os << " UPwedge";
  else if (is_wedge_down())
    os << " DNwedge";
  else
    os << " wedge??";

  return 1;
}

ostream &
operator << (ostream & os, const Bond & b)
{
  os << "Bond between " << b.a1 () << " and " << b.a2 () << " of type " << b.btype ();
  if (b.nrings () >= 0)
    os << ' ' << b.nrings () << " rings";
  if (b.is_directional_up ())
    os << " up";
  else if (b.is_directional_down ())
    os << " down";

  return os;
}

/*
  Change numbering to reflect the fact that atoms A1 and A2 have been
  swapped
  For example swap_atoms (2, 4)
    and
  _a1 = 2, _a2 = 3
  becomes
  _a1 = 4, _a2 = 4
*/

int
Bond::swap_atoms (atom_number_t a1, atom_number_t a2)
{
  int rc = 0;
  if (_a1 == a1)
  {
    _a1 = a2;
    rc++;
  }
  else if (_a1 == a2)
  {
    _a1 = a1;
    rc++;
  }

  if (_a2 == a1)
  {
    _a2 = a2;
    rc++;
  }
  else if (_a2 == a2)
  {
    _a2 = a1;
    rc++;
  }

  return rc;
}

/*
  Count the number of bonds represented by a bond.
  Note that double bonds are counted before aromatic bonds.
*/

int
Bond::number_of_bonds () const
{
  assert (OK_BOND_TYPE (_btype));

  if (IS_SINGLE_BOND (_btype))
    return 1;

  if (IS_DOUBLE_BOND (_btype))
    return 2;

  if (IS_TRIPLE_BOND (_btype))
    return 3;

  if (IS_AROMATIC_BOND (_btype))
    return 1;

  cerr << "What kind of bond is this " << _btype << endl;
  debug_print (cerr);
  iwabort ();

  return -1;
}

#ifdef BONDS_KNOW_RING_MEMBERSHIP

int
Bond::nrings (int & nr) const
{
  if (UNKNOWN_BOND_NRINGS == _nrings)
    return 0;

  assert (_nrings >= 0);

  nr = _nrings;
  return 1;
}

int
Bond::set_nrings (int nr)
{
  assert (nr >= 0);

  _nrings = nr;

  return 1;
}

int
Bond::nrings () const
{
  assert (UNKNOWN_BOND_NRINGS != _nrings);

  return _nrings;
}

int
Bond::in_another_ring ()
{
  if (UNKNOWN_BOND_NRINGS == _nrings)
    _nrings = 1;
  else
    _nrings++;

  return _nrings;
}

int
Bond::nrings_known () const
{
  if (UNKNOWN_BOND_NRINGS == _nrings)
    return 0;

  assert (_nrings >= 0);

  return 1;
}

void
Bond::invalidate_nrings ()
{
//_nrings = UNKNOWN_BOND_NRINGS;
  _nrings = 0;    // too hard otherwise

  set_non_aromatic ();    // if rings are unknown so too is aromaticity

//cerr << "After resetting bond between " << _a1 << " and " << _a2 << " now " << _btype << endl;

  return;
}

#endif

int
Bond::set_a1 (atom_number_t newa1)
{
  assert (newa1 >= 0);
  assert (newa1 != _a2);

  _a1 = newa1;

  return 1;
}

int
Bond::set_a2 (atom_number_t newa2)
{
  assert (newa2 >= 0);
  assert (newa2 != _a1);

  _a2 = newa2;

  return 1;
}

int
Bond::set_a1a2 (atom_number_t newa1, atom_number_t newa2)
{
  assert (newa1 >= 0);
  assert (newa2 >= 0);
  assert (newa1 != newa2);

  _a1 = newa1;
  _a2 = newa2;

  return 1;
}

/*
  We have various kinds of behaviour depending on the value of
  include_aromaticity_in_smiles. haphazardly using bits...,
  formalise if this gets used more...
*/

void
Bond::append_bond_type (IWString & smiles,
                        atom_number_t ato,
                        int include_aromaticity_in_smiles) const
{
  if (IS_AROMATIC_BOND(_btype))
  {
    if (write_smiles_aromatic_bonds_as_colons())
      smiles += AROMATIC_BOND_SYMBOL;
    else if (include_aromaticity_in_smiles)
      ;
    else if (IS_DOUBLE_BOND(_btype))
      smiles += DOUBLE_BOND_SYMBOL;
  }
  else if (IS_SINGLE_BOND (_btype))
  {
    if (include_aromaticity_in_smiles & 2)
      smiles += SINGLE_BOND_SYMBOL;
    else if (0 == _directional)
      ;
    else if (! include_cis_trans_in_smiles())
      ;
    else if (_directional & IW_BOND_DIRECTIONAL_UP)
    {
      if (_a1 == ato)
        smiles += '/';
      else
        smiles += '\\';
    }
    else if (_directional & IW_BOND_DIRECTIONAL_DOWN)
    {
      if (_a1 == ato)
        smiles += '\\';
      else
        smiles += '/';
    }
  }
  else if (IS_DOUBLE_BOND (_btype))
    smiles.add (DOUBLE_BOND_SYMBOL);
  else if (IS_TRIPLE_BOND (_btype))
    smiles.add (TRIPLE_BOND_SYMBOL);
  else
  {
    cerr << "Bond::append_bond_type:unrecognised bond type " << _btype << endl;
    smiles.add ('?');
  }

  return;
}

void
Bond::set_directional_up ()
{
  _directional = ((_directional & ~IW_BOND_DIRECTIONAL_DOWN) | IW_BOND_DIRECTIONAL_UP);

  return;
}

void
Bond::set_directional_down ()
{
  _directional = ((_directional & ~ IW_BOND_DIRECTIONAL_UP) | IW_BOND_DIRECTIONAL_DOWN);

  return;
}

/*
  The bond must be UP from s1 to s2
*/

void
Bond::set_directional_up (atom_number_t s1, atom_number_t s2)
{
  if (s1 == _a1)
    set_directional_up ();
  else
    set_directional_down ();

  return;
}

void
Bond::set_directional_down (atom_number_t s1, atom_number_t s2)
{
  if (s1 == _a1)
    set_directional_down ();
  else
    set_directional_up ();

  return;
}

void
Bond::set_wedge_up ()
{
  int mask = ~ (IW_BOND_DIRECTION_WEDGE_DOWN | IW_BOND_DIRECTION_WEDGE_EITHER);

  _directional = (_directional & mask);

  _directional = (_directional | IW_BOND_DIRECTION_WEDGE_UP);

//cerr << "set_wedge_up: atoms " << _a1 << " to " << _a2 << " _directional now " << hex << _directional << dec << endl;

  return;
}

void
Bond::set_wedge_down ()
{
  int mask = ~ (IW_BOND_DIRECTION_WEDGE_UP | IW_BOND_DIRECTION_WEDGE_EITHER);

  _directional = (_directional & mask);

  _directional = (_directional | IW_BOND_DIRECTION_WEDGE_DOWN);

//cerr << "set_wedge_down: atoms " << _a1 << " to " << _a2 << " _directional now " << hex << _directional << dec << endl;

  return;
}

void
Bond::set_wedge_either ()
{
  int mask = ~ (IW_BOND_DIRECTION_WEDGE_UP | IW_BOND_DIRECTION_WEDGE_DOWN);

  _directional = (_directional & mask);

  _directional = (_directional | IW_BOND_DIRECTION_WEDGE_EITHER);

  return;
}

void
Bond::copy_directionality_specifications (const Bond * rhs)
{
  _directional = rhs->_directional;

  return;
}

int
Bond::set_part_of_cis_trans_grouping (int s)
{
  if (0 == s)
    _directional = (_directional ^ IW_BOND_DIRECTIONAL_DOUBLE_BOND);
  else if (0 == (_btype | DOUBLE_BOND))
  {
    cerr << "Bond::set_part_of_cis_trans_grouping:can only apply to double bonds\n";
    return 0;
  }
  else
    _directional = (_directional | IW_BOND_DIRECTIONAL_DOUBLE_BOND);

  return 1;
}

int
Bond::same_bond_type (const Bond & rhs) const
{
  if (IS_AROMATIC_BOND (_btype))
    return (IS_AROMATIC_BOND (rhs._btype));

  return _btype == rhs._btype;
}
