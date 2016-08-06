#!/usr/bin/python

# Copyright (c) 2002, Trevor Blackwell. All rights reserved.
# The author can be reached at tlb@tlb.org.
# 
# (This is the MIT License for open source software with capitalization
# fixed.)
# 
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions: 
# 
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# 
# The software is provided "as is", without warranty of any kind,
# express or implied, including but not limited to the warranties of
# merchantability, fitness for a particular purpose and noninfringement.
# In no event shall the authors or copyright holders be liable for any
# claim, damages or other liability, whether in an action of contract,
# tort or otherwise, arising from, out of or in connection with the
# software or the use or other dealings in the software.

# Generate file mtl_specializations.{cc,h} which implement fast versions
# of mult_nonblocked(matrix,matrix,matrix) and mult_nonblocked(matrix,vector,vector).
# In a perfect world, all C++ compilers would be able to collapse the many
# layers of templating in MTL down to a simple inner loop in the matrix
# routines; however g++ as of version 3.1 can't.
# These versions runs about 10x faster than what g++ generates.
# To use, run this script and:
#   - #include "mtl_specializations.h" in any file where you might call these
#   - compile mtl_specialzations.cc into your program


from cgen import *

# Stuff you may want to change:

eltypes=('double','float')
# you might need: eltypes=('double','float','complex<double>','complex<float>')

class matdef:
    def __init__(self,elt,orien,scale,const):
        self.elt=elt
        if (elt=='double'): self.shortelt='d'
        if (elt=='float'): self.shortelt='s'
        if (elt=='complex<double>'): self.shortelt='z'
        if (elt=='complex<float>'): self.shortelt='c'
        self.orien=orien
        self.scale=scale
        self.const=const
    def decl(self):
        r='%s_matrix<gen_dense2D<%s, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<%s_orien, 0, 0, size_t> >' % (
            self.orien,self.elt,self.orien)
        if (self.scale):
            r='%s::scaled_type' % (r)
        return r
    def getnrows(self):
        return '%s.nrows()' % self.name
    def getncols(self):
        return '%s.ncols()' % self.name
    # Get the row/col stride of the matrix. For a row-major matrix,
    # the column stride is 1 and vice-versa
    def getrowstride(self):
        if (self.orien=='column'): return '1'
        if (self.orien=='row'): return '%s.ncols()' % self.name
    def getcolstride(self):
        if (self.orien=='column'): return '%s.nrows()' % self.name
        if (self.orien=='row'): return '1'
    def lead(self):
        if (self.orien=='row'): return '%s.ncols()' % self.name
        if (self.orien=='column'): return '%s.nrows()' % self.name
    def gettrans(self,other):
        if (self.orien==other.orien):
            return 'CblasNoTrans'
        else:
            return 'CblasTrans'
    def order(self):
        if (self.orien=='row'): return 'CblasRowMajor' 
        if (self.orien=='column'): return 'CblasColMajor'
    def data(self):
        if (self.scale):
            return '%s.get_twod().twod.data()' % (self.name)
        else:
            return '%s.data()' % (self.name)
    def pr_setup(self,f):
        f('size_t %s_rs=%s;' % (self.name,self.getrowstride()))
        f('size_t %s_cs=%s;' % (self.name,self.getcolstride()))
        if (self.scale):
            f('%s %s *%s_d=%s.get_twod().twod.data();' % (self.const,self.elt,self.name,self.name))
        else:
            f('%s %s *%s_d=%s.data();' % (self.const,self.elt,self.name,self.name))
        if 0:
            f('fprintf(stderr,"%s: rs=%%d cs=%%d d=%%p\\n",%s_rs,%s_cs,%s_d);' % (
                self.name,self.name,self.name,self.name))

def matdef_poss(const):
    ret=[]
    for elt in eltypes:
        for orien in ('row','column'):
            for scale in const and (0,1) or (0,):
                ret.append(matdef(elt,orien,scale,const))
    return ret
        
class vecdef:
    def __init__(self,elt,scale,const):
        self.elt=elt
        self.scale=scale
        self.const=const
    # resolved type declaration
    def decl(self):
        r='dense1D< %s >' % (self.elt)
        if (self.scale):
            r='%s::scaled_type' % (r)
        return r
    # Get the row/col strides and data pointer
    def pr_setup(self,f):
        if (self.scale):
            f('%s %s *%s_d=%s.rep.data();' % (self.const,self.elt,self.name,self.name))
        else:
            f('%s %s *%s_d=%s.data();' % (self.const,self.elt,self.name,self.name))

def vecdef_poss(const):
    ret=[]
    for elt in eltypes:
        for scale in (0,1):
            ret.append(vecdef(elt,scale,const))
    return ret

def mult_mmm(m1,m2,m3):
    m1.name='m1'
    m2.name='m2'
    m3.name='m3'
    decl='void mult(const %s &m1, const %s &m2, %s &m3)' % (
        m1.decl(),
        m2.decl(),
        m3.decl())
    h('%s;' % decl)
    c('%s' % decl)
    c('{')
    if (0):
        c('}')
        return
    c('size_t nr=%s;' % m1.getnrows())
    c('size_t nc=%s;' % m2.getncols())
    c('size_t k=%s;' % m1.getncols())
    c('assert(k==%s);' % m2.getnrows())
    c('assert(nr==%s);' % m3.getnrows())
    c('assert(nc==%s);' % m3.getncols())

    c('if (nr*nc*k < N_CBLAS) {')

    m1.pr_setup(c)
    m2.pr_setup(c)
    m3.pr_setup(c)

    facop=''
    if (m1.scale or m2.scale):
        facop='*fac'
        c('%s fac=%s*%s;' % (m3.elt,
                             m1.scale and 'm1.get_twod().alpha' or '1.0',
                             m2.scale and 'm2.get_twod().alpha' or '1.0'))

    c('for (size_t r=0; r<nr; r++) {')
    c('%s *m3p=m3_d+r*m3_rs;' % m3.elt)
    c('for (size_t c=0; c<nc; c++) {')
    c('const %s *m1p=m1_d+r*m1_rs;' % m1.elt)
    c('const %s *m2p=m2_d+c*m2_cs;' % m2.elt)
    c('%s tot=0.0;' % m3.elt)
    c('for (size_t i=0; i<k; i++) {')
    c('tot += *m1p * *m2p;')
    c('m1p+=m1_cs;')
    c('m2p+=m2_rs;')
    c('}')
    c('*m3p+=tot%s;' % facop)
    c('m3p+=m3_cs;')
    c('}')
    c('}')
    c('} else {')
    c('/* use cblas */')
    c('const %s ONE=%s(1.0);' % (m3.elt, m3.elt) )
    c('const %s alpha=(%s * %s);' % 
        (m3.elt,(m1.scale and 'm1.get_twod().alpha' or 'ONE'),
        (m2.scale and 'm2.get_twod().alpha' or 'ONE')))
    c('cblas_%sgemm(%s, %s, %s, nr, nc, k, %s, %s, %s, %s, %s, %s, %s, %s);' %
      ( m3.shortelt,
        m3.order(),
        m1.gettrans(m3),
        m2.gettrans(m3),
        (match(r'complex',m3.elt) and '&alpha' or 'alpha'),
        m1.data(),
        m1.lead(),
        m2.data(),
        m2.lead(),
        (match(r'complex',m3.elt) and '&ONE' or 'ONE'),
        m3.data(),
        m3.lead() ))
    c('}')
    c('}')
    c.sep()


def mult_mvv(m1,v2,v3):
    m1.name='m1'
    v2.name='v2'
    v3.name='v3'
    decl='void mult(const %s &m1, const %s &v2, %s &v3)' % (
        m1.decl(),
        v2.decl(),
        v3.decl())
    h('%s;' % decl)
    c('%s' % decl)
    c('{')
    if (0):
        c('}')
        return
    c('typedef %s E;' % m1.elt)
    c('size_t nr=%s;' % m1.getnrows())
    c('size_t nc=%s;' % m1.getncols())
    c('assert(nr==v3.size());')
    c('assert(nc==v2.size());')

    m1.pr_setup(c)
    v2.pr_setup(c)
    v3.pr_setup(c)

    facop=''
    if (m1.scale or v2.scale):
        facop='*fac'
        c('%s fac=%s*%s;' % (v3.elt,
                             m1.scale and 'm1.get_twod().alpha' or '1.0',
                             v2.scale and 'v2.scale' or '1.0'))

    c('%s *v3p=v3_d;' % v3.elt)
    c('for (size_t r=0; r<nr; r++) {')
    c('const %s *m1p=m1_d+r*m1_rs;' % m1.elt)
    c('const %s *v2p=v2_d;' % v2.elt)
    c('%s tot=0.0;' % v3.elt)
    c('for (size_t c=0; c<nc; c++) {')
    c('tot += *m1p * *v2p;')
    c('m1p+=m1_cs;')
    c('v2p++;')
    c('}')
    c('*v3p=tot%s;' % facop)
    c('v3p++;')
    c('}')
    c('}')
    c.sep()

# gw changed this:
def compatible_types(t1,t2,t3):
#    if ((match(r'complex',t1.elt) or match(r'complex',t2.elt)) and
#        not match(r'complex',t3.elt)): return 0
    if ((not t1.elt==t2.elt) or (not t1.elt==t3.elt) or (not t2.elt==t3.elt)):
        return 0
    return 1

cp=cproj('mtl_specializations')
h=cp.h
c=cp.c

c('''
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#define protected public  // heh heh
#include <complex>
#include <mtl/mtl.h>
#include <mtl/matrix.h>
#include <mtl/blais.h>
#include <mtl/dense1D.h>
#include "mtl_specializations.h"

////
// Implemententation of  matrix specializations using atlas/blais functionality
////
#ifdef HAVE_BLAS
#include "BlasHeaders.h"

/* minimal M*N*K to use cblas routines */
#define N_CBLAS 2000

namespace mtl {
''')

h('namespace mtl {')

parms={}
for m1 in matdef_poss('const'):
    for m2 in matdef_poss('const'):
        for m3 in matdef_poss(''):
            if (compatible_types(m1,m2,m3)):
                mult_mmm(m1,m2,m3)

for m1 in matdef_poss('const'):
    for v2 in vecdef_poss('const'):
        for v3 in vecdef_poss(''):
            if (compatible_types(m1,v2,v3)):
                mult_mvv(m1,v2,v3)

h('}')
c('}')
