//-----------------------------------------------------------------------------
// ElementHeaders.h
//
// begin     : Dec 6 2001
// copyright : (c) 2001 by Oliver Koenig, Marc Wintermantel, Nino Zehnder
// email     : {okoenig, wintermantel, nzehnder}@imes.mavt.ethz.ch
// www       : www.structures.ethz.ch
/* 
   This file is part of FELyX (Finite Element Library eXperiment).
   
   FELyX is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   FELyX is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with FELyX; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
//-----------------------------------------------------------------------------
   
#ifndef ElementHeaders_h
#define ElementHeaders_h ElementHeaders_h

////
//// Include all element types implemented
////
#include "Link1.h"
#include "Plane2.h"
#include "Beam3.h"
#include "Beam4.h"
#include "Link8.h"
#include "Plane182.h"
#include "Plane183.h"
#include "Solid185.h"
#include "Solid186.h"
#include "Solid187.h"
#include "Shell93.h"
#include "Shell91.h"
#include "Shell393.h"
#include "Shell391.h"
#include "Darcy2D3.h"
#include "Darcy3D3.h"
#include "Darcy3D4.h"
#include "Darcy3D8.h"

#include "LcmElement.h"



////
//// Register all element types in element factory
////

#include "ElementFactory.h"
using namespace fe_base;

// Create an anonymous namespace to make the functions invisible 
// from other modules
namespace 
{
  // Define create functions for all different elements
  StructElement* Create1()	{ return new Link1; 	}
  StructElement* Create2()	{ return new Plane2;	}
  StructElement* Create3()	{ return new Beam3; 	}
  StructElement* Create4()	{ return new Beam4; 	}
  StructElement* Create8()	{ return new Link8; 	}
  StructElement* Create93()	{ return new Shell93;	}
  StructElement* Create393()	{ return new Shell393;	}
  StructElement* Create91()	{ return new Shell91;	}
  StructElement* Create391()	{ return new Shell391;	}
  StructElement* Create182()	{ return new Plane182;	}
  StructElement* Create183()	{ return new Plane183;	}
  StructElement* Create185()	{ return new Solid185;	}
  StructElement* Create186()	{ return new Solid186;	}
  StructElement* Create187()	{ return new Solid187;	}
  LcmElement* Create55()       { return new Darcy2D3; }
  LcmElement* Create57()       { return new Darcy3D3; }
  LcmElement* Create470()       { return new Darcy3D4; }
  LcmElement* Create70()       { return new Darcy3D8; }

  
  // Register ID's and create functions in StructElementFactory
  const bool reg1    = StructElementFactory::Instance().Register( Link1::Id,    Create1    );
  const bool reg2    = StructElementFactory::Instance().Register( Plane2::Id,   Create2    );
  const bool reg3    = StructElementFactory::Instance().Register( Beam3::Id,    Create3    );
  const bool reg4    = StructElementFactory::Instance().Register( Beam4::Id,    Create4    );
  const bool reg8    = StructElementFactory::Instance().Register( Link8::Id,    Create8    );
  const bool reg93   = StructElementFactory::Instance().Register( Shell93::Id,  Create93   );
  const bool reg393  = StructElementFactory::Instance().Register( Shell393::Id, Create393  );
  const bool reg91   = StructElementFactory::Instance().Register( Shell91::Id,  Create91   );
  const bool reg391  = StructElementFactory::Instance().Register( Shell391::Id, Create391  );
  const bool reg182  = StructElementFactory::Instance().Register( Plane182::Id, Create182  );
  const bool reg183  = StructElementFactory::Instance().Register( Plane183::Id, Create183  );
  const bool reg185  = StructElementFactory::Instance().Register( Solid185::Id, Create185  );
  const bool reg186  = StructElementFactory::Instance().Register( Solid186::Id, Create186  );
  const bool reg187  = StructElementFactory::Instance().Register( Solid187::Id, Create187  );


  // Register ID's and create functions in LcmElementFactory
  const bool reg55  = LcmElementFactory::Instance().Register( Darcy2D3::Id, Create55  );
  const bool reg57  = LcmElementFactory::Instance().Register( Darcy3D3::Id, Create57  );
  const bool reg470  = LcmElementFactory::Instance().Register( Darcy3D4::Id, Create470  );
  const bool reg70  = LcmElementFactory::Instance().Register( Darcy3D8::Id, Create70  );
}

#endif
