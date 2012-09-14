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
#include "ostream_and_type.h"

#include "molecule.h"

int
ofstream_and_type::_default_values ()
{
  _valid = 0;

  _molecules_written = 0;
  _verbose = 0;

  _output_type = -1;

  return 1;
}

ofstream_and_type::ofstream_and_type ()
{
  _default_values ();
}

ofstream_and_type::ofstream_and_type (int output_type)
{
  if (_default_values ())
    return;

  if (! set_type (output_type))
    return;

  _valid = 1;

  return;
}

ofstream_and_type::ofstream_and_type (int output_type, const char * fname)
{
  if (!_default_values ())
    return;

  if (! set_type (output_type))
    return;

  if (open (fname))
    _valid = 1;
}

ofstream_and_type::ofstream_and_type (int output_type, IWString & fname)
{
  if (!_default_values ())
    return;

  if (! set_type (output_type))
    return;

  if (open (fname.null_terminated_chars ()))
    _valid = 1;
}

ofstream_and_type::~ofstream_and_type ()
{
//if (_valid && good () && BFILE == _output_type)
//  (*this) << "-1\n";

  return;
}

int
ofstream_and_type::set_type (int t)
{
  if (valid_file_type (_output_type))
  {
    cerr << "ofstream_and_type::set_type: file type already set to " << _output_type << endl;
    return 0;
  }

  if (! valid_file_type (t))
  {
    cerr << "ofstream_and_type::set_type: type " << t << " is invalid\n";
    return 0;
  }

  _output_type = t;

  _valid = 1;

  return 1;
}

int
ofstream_and_type::open (const char * fname)
{
  if ('>' == *fname && strlen (fname) > 1 && '>' == fname[1])
  {
    fname += 2;
    ofstream::open (fname, ios::app);
  }
  else
    ofstream::open (fname, ios::out);

  if (good ())
    _fname = fname;

//cerr << "ofstream_and_type: opened '" << fname << "' good = " << good () << endl;

  return good ();
}

int
ofstream_and_type::open (IWString & fname)
{
  if (fname.starts_with (">>"))
  {
    fname.remove_leading_chars (2);
    ofstream::open (fname.chars (), ios::app);
  }
  else
    ofstream::open (fname.chars (), ios::out);

  if (good ())
    _fname = fname;

  return good ();
}

int
ofstream_and_type::write_molecule (Molecule * m)
{
  if (! _valid)
    return 0;

  if (! good ())
  {
    _valid = 0;
    return 0;
  }

  if (_verbose)
    cerr << _molecules_written + 1 << " writing '" << m->name () << "'\n";

  if (! m->write_molecule (*this, _output_type))
    return 0;

  _molecules_written++;
  return 1;
}

int
ofstream_and_type::write_molecules (const resizable_array_p<Molecule> & molecules)
{
  if (! _valid)
    return 0;

  if (! good ())
  {
    _valid = 0;
    return 0;
  }

  int rc = 0;
  for (int i = 0; i < molecules.number_elements (); i++)
  {
    Molecule * m = molecules[i];
    assert (OK_MOLECULE (m));

    if (! m->write_molecule (*this, _output_type))
      return rc;

    rc++;
  }

  return rc;
}
