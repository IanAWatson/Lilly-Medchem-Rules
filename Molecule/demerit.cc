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
#include "assert.h"

#include "demerit.h"

static int demerit_reason_contains_individual_demerits = 0;

void
set_demerit_reason_contains_individual_demerits (int s)
{
  demerit_reason_contains_individual_demerits = s;

  return;
}

static int store_demerit_reasons_like_tsubstructure = 0;

void
set_store_demerit_reasons_like_tsubstructure (int s)
{
  store_demerit_reasons_like_tsubstructure = s;

  return;
}

static int _rejection_threshold = DEFAULT_REJECTION_THRESHOLD;

void
set_rejection_threshold (int s)
{
  _rejection_threshold = s;
}

int
rejection_threshold ()
{
  return _rejection_threshold;
}

Demerit::Demerit ()
{
  _score = 0;
  _number_different_demerits_applied = 0;

  return;
}

int
Demerit::debug_print (ostream & os) const
{
  os << "Demerit: " << _number_different_demerits_applied << " demerits, ";
  if (_score >= _rejection_threshold)
    os << "REJECTED";
  else
    os << "total " << _score;

  os << ", origin '" << _types << "'\n";

  return 1;
}

int
Demerit::rejected () const
{
  return _score >= _rejection_threshold;
}

void
Demerit::_increment (int increment)
{
  assert (increment > 0);
  if (increment <= 0)
  {
    cerr << "Demerit::_increment:zero increment ignored\n";
    return;
  }

  _score += increment;
  _number_different_demerits_applied++;

  return;
}

int 
Demerit::reject (const const_IWSubstring reason)
{
  _increment (_rejection_threshold);

  _add_hit_type (_rejection_threshold, reason);

  return 1;
}

int
Demerit::extra (int increment, const const_IWSubstring reason)
{
  _increment (increment);

  _add_hit_type (increment, reason);

  return 1;
}

void
Demerit::_add_hit_type (int increment,
                        const const_IWSubstring & reason)
{
  if (store_demerit_reasons_like_tsubstructure)
  {
    if (_types.length ())
      _types.add (' ');

    _types << "(1 matches to 'D" << increment << ' ' << reason << "')";
//  _types << "(1 matches to '" << reason << "')";
  }
  else if (demerit_reason_contains_individual_demerits)
  {
    if (_types.length ())
      _types.add (':');

    _types << 'D' << increment << ' ' << reason;
  }
  else
    _types.append_with_spacer(reason, ':');

  return;
}

int
Demerit::write_in_tdt_form (ostream & os) const
{
  if (_score >= _rejection_threshold)
    os << "REJ<1>\n";

  os << "DMRT<" << _score << ">\n";
  os << "NDMRT<" << _number_different_demerits_applied << ">\n";
  os << "DMRTYP<" << _types << ">\n";

  return 1;
}
