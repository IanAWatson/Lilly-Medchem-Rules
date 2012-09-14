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
#ifndef IWCONFIG_H
#define IWCONFIG_H

#ifdef _WIN32

#define		_CRT_SECURE_NO_DEPRECATE 1		//Ignore Microsoft security deprecation
	#include	<process.h>
	#include	<io.h>
	typedef int	ssize_t;
	typedef		_int64 INT64;
	#define		String2I64(cTmp) _atoi64(cTmp);
	#define		ios_in_nocreate ios::in
	#define		IW_INTEL_COMPILER
#define IW_STRCPY ::strcpy
#define IW_STRNCPY ::strncpy
#define IW_FD_OPEN ::_open
#define IW_FD_READ ::_read
#define IW_FD_LSEEK ::_lseek
#define IW_FD_WRITE ::_write
#define IW_FD_CLOSE ::_close
#define IW_SPRINTF ::sprintf_s
#define IW_SSCANF ::sscanf_s
#define IW_GETPID ::_getpid
#define IW_IMPLEMENTATIONS_EXPOSED 1

#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <limits.h>
#include <cassert>
#include <assert.h>

#else

#define IW_IMPLEMENTATIONS_EXPOSED 1
#define IW_STRCPY ::strcpy
#define IW_STRNCPY ::strncpy

#define IW_FD_OPEN ::open
#define IW_FD_READ ::read
#define IW_FD_LSEEK ::lseek
#define IW_FD_WRITE ::write
#define IW_FD_CLOSE ::close
#define IW_SPRINTF ::sprintf
#define IW_SSCANF ::sscanf
#define IW_GETPID ::getpid

#define	String2I64(cTmp) atoll(cTmp);
typedef long long INT64;

#endif

#ifdef __INTEL_COMPILER
#define IW_INTEL_COMPILER
#endif

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 \
                     + __GNUC_MINOR__ * 100 \
                     + __GNUC_PATCHLEVEL__)

#endif

#endif
