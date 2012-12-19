#ifndef _libgeodecomp_parallelization_distributedsimulator_h_
#define _libgeodecomp_parallelization_distributedsimulator_h_

#include <libgeodecomp/io/parallelwriter.h>
#include <libgeodecomp/misc/displacedgrid.h>
#include <libgeodecomp/misc/supervector.h>
#include <libgeodecomp/parallelization/simulator.h>

namespace LibGeoDecomp {

template<typename CELL_TYPE>
class ParallelWriter; 

// fixme: add short doxygen doc for every class
template<typename CELL_TYPE>
class DistributedSimulator : public Simulator<CELL_TYPE>
{
public:
    typedef typename CELL_TYPE::Topology Topology;
    typedef GridBase<CELL_TYPE, Topology::DIMENSIONS> GridType;
    typedef SuperVector<boost::shared_ptr<ParallelWriter<CELL_TYPE> > > WriterVector;

    inline DistributedSimulator(Initializer<CELL_TYPE> *_initializer) : 
        Simulator<CELL_TYPE>(_initializer)
    {}

    // fixme: this method is no longer required. purge it.
    /**
     * Returns the fragment of the current grid, which is available to
     * the current instance.
     */
    virtual void getGridFragment(
        const GridType **grid, 
        const Region<Topology::DIMENSIONS> **validRegion) = 0;

    /**
     * register @a writer which will observe the simulation. The
     * DistributedSimulator will assume that it now owns the
     * ParallelWriter, so it'll delete it upon destruction.
     *
     * fixme: replace @a by \a and @return by \returns ...
     */
    virtual void addWriter(ParallelWriter<CELL_TYPE> *writer)
    {
        writers << boost::shared_ptr<ParallelWriter<CELL_TYPE> >(writer);
    }

protected:
    WriterVector writers;

};

}

#endif
