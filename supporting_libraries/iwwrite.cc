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
#ifdef _WIN32
#else
#include <unistd.h>
#endif
#include <sys/types.h>
#include <iostream>
using namespace std;

#include "iwstring.h"

static int
common_write (int fd, const char * s, int nchars)
{
  assert (fd >= 0);

  if (0 == nchars)
    return 1;

  int rc = IW_FD_WRITE (fd, s, nchars);

  if (rc == nchars)
    return 1;

  cerr << "iwstring::common_write: cannot write " << nchars << " bytes to fd " << fd << endl;

  return 0;
}

int
const_IWSubstring::write (int fd) const
{
  return common_write (fd, _data, _nchars);
}

int
IWString::write (int fd) const
{
  return common_write (fd, _things, _number_elements);
}
