#ifndef LIBGEODECOMP_TESTBED_PERFORMANCETESTS_CELL_H
#define LIBGEODECOMP_TESTBED_PERFORMANCETESTS_CELL_H

#include <libgeodecomp/misc/apitraits.h>

#include <libflatarray/macros.hpp>

using namespace LibGeoDecomp;
using namespace LibFlatArray;

class SoACell
{
public:
    class API :
        public APITraits::HasStencil<Stencils::VonNeumann<3, 1> >,
        public APITraits::HasTorusTopology<3>,
        public APITraits::HasSoA
    {};

    double c;
    int a;
    char b;
};

LIBFLATARRAY_REGISTER_SOA(SoACell, ((double)(c))((int)(a))((char)(b)))

class Cell
{
public:
    class API :
        public APITraits::HasStencil<Stencils::VonNeumann<3, 1> >,
        public APITraits::HasTorusTopology<3>
    {};

    double c;
    int a;
    char b;
};

#endif
