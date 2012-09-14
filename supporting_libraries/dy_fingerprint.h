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
#ifndef IW_DAYLIGHT_FINGERPRINT_H
#define IW_DAYLIGHT_FINGERPRINT_H

extern int du_ascii2bin(const char *ascii, int nchars,
             unsigned char * binary, unsigned int & nbytes);

extern char * du_bin2ascii (int *palen, int blen, char *b);

#endif       
