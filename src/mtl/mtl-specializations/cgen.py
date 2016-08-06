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

import sys
from re import *
import inspect

class cgen:
    def __init__(self,f):
        self.file=f
        self.indent=0
        self.pylines=0
    def line(self,code):
        if (callable(code)):
            code(self)
        else:
            for rawline in code.split('\n'):
                line=sub(r'^\s+','',rawline)
                if match(r'\#',line):
                    self.file.write(line+'\n')
                else:
                    if match(r'\s*}',line):
                        self.indent-=1
                    if self.indent<0:
                        raise ValueError,"Negative indent"
                    self.file.write('  ' * self.indent)
                    self.file.write(line)
                    self.file.write('\n')
                    if (match(r'.*{\s*$',line)):
                        self.indent+=1
    def __call__(self,code):
        if (self.pylines):
            scope=inspect.getouterframes(inspect.currentframe())
            frame,filename,lineno,fname,context,linei=scope[1]
            self.line('# %d \"%s\"' % (lineno,filename))
            # Otherwise, there is a circular reference among stack frames
            del frame
            del scope
        self.line(code)
    def sep(self):
        self.file.write('\n')

class cproj:
    def __init__(self,basename):
        self.basename=basename
        self.c=cgen(open("%s.cc"%basename,"w"))
        self.h=cgen(open("%s.h"%basename,"w"))
        self.c('// COPYRIGHT')
        self.h('// COPYRIGHT')
        self.elems=[]
    def add(self,x):
        self.elems.append(x)
    def emit(self):
        self.h('#ifndef _%s_H' % upper(self.basename))
        self.h('#define _%s_H' % upper(self.basename))
        for x in self.elems:
            x.emit_decl(self.h)
            x.emit_def(self.c)
        self.h('#endif')

class cgen_item:
    def def_arg(self,a):
        return sub(r'=.*$','',a)
    def decl_arg(self,a):
        return a
    def emit_decl(self,f):
        f('// no declaration for %s'%self.name)
    def emit_def(self,f):
        f('// no definition for %s'%self.name)

class cgen_class(cgen_item):
    def __init__(self,name,super=None):
        self.name=name
        self.super=super
        self.decls=[]
        self.data=[]
        self.code=[]
    def add_data(self,type,name):
        self.data.append(Struct(name=name,type=type))
    def add_code(self,name,rtype,args,code):
        self.code.append(Struct(name=name,rtype=rtype,args=args,code=code))
    def add_decl(self,text):
        self.decls.append(text)
    def emit_decl(self,f):
        f('class %s %s {' %(self.name,
                            self.super and ':public %s'%self.super or ''))

        for x in self.decls:
            f(x)
        
        for x in self.code:
            f('%s %s(%s);' % (x.rtype,x.name,
                              ', '.join([self.decl_arg(a) for a in x.args])))
            
        for x in self.data:
            f('%s %s;' % (x.type,x.name))
        
        f('};\n\n')
    def emit_def(self,f):
        for x in self.code:
            f('%s %s::%s(%s) {' % (x.rtype,self.name,x.name,
                              ', '.join([self.def_arg(a) for a in x.args])))
            f(x.code)
            f('}\n\n')

