#ifndef LIBGEODECOMP_IO_MPIIOINITIALIZER_H
#define LIBGEODECOMP_IO_MPIIOINITIALIZER_H

#include <libgeodecomp/config.h>
#ifdef LIBGEODECOMP_WITH_MPI

#include <libgeodecomp/communication/typemaps.h>
#include <libgeodecomp/io/initializer.h>
#include <libgeodecomp/io/mpiio.h>

namespace LibGeoDecomp {

/**
 * This class can set up a simulation from a snapshot generated by
 * MPIIOWriter/ParallelMPIIOWriter. This is especially useful for
 * long-running jobs which might either be shot down because of wall
 * clock limitations or node failures: here checkpoints can save
 * captital amounts of compute time.
 */
template<typename CELL_TYPE>
class MPIIOInitializer : public Initializer<CELL_TYPE>
{
public:
    typedef typename APITraits::SelectTopology<CELL_TYPE>::Value Topology;
    static const int DIM = Topology::DIM;

    explicit MPIIOInitializer(
        const std::string& filename,
        const MPI_Datatype& mpiDatatype = Typemaps::lookup<CELL_TYPE>(),
        const MPI_Comm& comm = MPI_COMM_WORLD) :
        file(filename),
        datatype(mpiDatatype),
        communicator(comm)
    {
        mpiio.readMetadata(
            &dimensions, &currentStep, &maximumSteps, file, communicator);
    }

    virtual void grid(GridBase<CELL_TYPE, DIM> *target)
    {
        Region<DIM> region;
        region << target->boundingBox();
        mpiio.readRegion(target, file, region, communicator, datatype);
    }

    virtual Coord<DIM> gridDimensions() const
    {
        return dimensions;
    }

    virtual unsigned maxSteps() const
    {
        return maximumSteps;
    }

    virtual unsigned startStep() const
    {
        return currentStep;
    }

private:
    std::string file;
    MPI_Datatype datatype;
    MPI_Comm communicator;
    MPIIO<CELL_TYPE> mpiio;
    unsigned currentStep;
    unsigned maximumSteps;
    Coord<DIM> dimensions;
};

}

#endif
#endif
