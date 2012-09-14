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
#include <sys/types.h>
#ifdef _WIN32
#else
#include <unistd.h>
#endif

#include "iwconfig.h"
#include "iwstring.h"

#define IWWRITE_BLKSIZE 4096

int
IWString::write_whole_blocks_shift_unwritten (int fd)
{
  if (_number_elements < IWWRITE_BLKSIZE)
    return 1;


  int blocks_to_write = _number_elements / IWWRITE_BLKSIZE;

  int chars_written = IW_FD_WRITE (fd, _things, blocks_to_write * IWWRITE_BLKSIZE);

  if (chars_written != blocks_to_write * IWWRITE_BLKSIZE)
  {
    cerr << "IWString::write_whole_blocks_shift_unwritten:cannot write " << blocks_to_write << " blocks to " << fd << endl;
    return 0;
  }

  ::memcpy (_things, _things + chars_written, _number_elements - chars_written);

  _number_elements -= chars_written;

  return chars_written;
}
