/***************************************************************************
 *   Copyright (C) 2004 by Oliver Koenig, Marc Wintermantel, Nino Zehnder  *
 *   {koenig,wintermantel,zehnder}@even-ag.ch                              *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef BLASHEADERS_H
#define BLASHEADERS_H BLASHEADERS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_ATLAS
extern "C"{ 
#include <cblas.h>
}
#elif HAVE_MKL
#include <mkl_cblas.h>
#endif

#endif
