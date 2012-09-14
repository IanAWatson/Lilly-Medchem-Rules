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
#ifndef IW_HASH_H
#define IW_HASH_H

#include "iwconfig.h"

class IWStringHash
{
  private:
  public:

#if defined (IW_INTEL_COMPILER)

    static const size_t bucket_size = 4;
    static const size_t min_buckets = 8;
    bool  operator () (const IWString &, const IWString &) const;

#endif

    size_t operator () (const IWString &) const;
};

#endif
