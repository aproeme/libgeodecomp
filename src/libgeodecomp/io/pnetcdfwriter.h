#ifndef LIBGEODECOMP_IO_PNETCDFWRITER_H
#define LIBGEODECOMP_IO_PNETCDFWRITER_H

#include <libgeodecomp/config.h>
#ifdef LIBGEODECOMP_WITH_MPI

#include <iostream>
#include <mpi.h>

// To do: include pnetcdf once in libgeodecomp.h
//        instead of working with include guards
#ifndef LIBGEODECOMP_INCLUDED_PNETCDF
#define LIBGEODECOMP_INCLUDED_PNETCDF
#include <pnetcdf>
#endif

#include <libgeodecomp/communication/mpilayer.h>
#include <libgeodecomp/io/parallelwriter.h>
#include <libgeodecomp/misc/clonable.h>
#include <libgeodecomp/storage/selector.h>
#include <libgeodecomp/io/mpiio.h>

#include <filesystem>

using namespace PnetCDF;
using namespace PnetCDF::exceptions;

namespace LibGeoDecomp{

/**
 *
 * Writes simulation variables in parallel using pnetcdf (https://github.com/Parallel-NetCDF/PnetCDF)
 * Currently using NetCDF's CDF-2 64-bit offset format (PnetCDF::NcmpiFile::FileFormat::classic2)
 * Aims to adhere to the NetCDF CF convention data model (http://cfconventions.org/)
 * Modelled on BOVWriter
 */


template<typename CELL_TYPE>
class PnetCDFWriter : public Clonable<ParallelWriter<CELL_TYPE>, PnetCDFWriter<CELL_TYPE> >
{
public:
    using ParallelWriter<CELL_TYPE>::period;
    using ParallelWriter<CELL_TYPE>::prefix;
    
    typedef typename ParallelWriter<CELL_TYPE>::Topology Topology;
    
    static const int DIM = Topology::DIM;
    
    PnetCDFWriter(
	const Coord<DIM>& globalDimensions,
	const Selector<CELL_TYPE>& selector,
	const std::string& prefix,
	const unsigned period,
	const unsigned steps,
	const MPI_Comm& communicator = MPI_COMM_WORLD) :
	Clonable<ParallelWriter<CELL_TYPE>, PnetCDFWriter<CELL_TYPE> >(prefix, period),
	globalDimensions(globalDimensions),
	selector(selector),
	steps(steps),
	communicator(communicator)
	{
	    filename = selector.name() + ".nc";
	    restartFilename = selector.name() + "_restart.nc";
	    
	    createFile();
	    writeHeader();

	    createRestartFile();
	    writeRestartHeader();

	    writeCount = 0;
	}
  
  // Only gets called by simulator when step == period for a given PnetCDFWriter's period
    virtual void stepFinished(
        const typename ParallelWriter<CELL_TYPE>::GridType& localGrid,
        const Region<Topology::DIM>& validRegion,
        const Coord<Topology::DIM>& globalDimensions,
        unsigned step,
        WriterEvent event,
        std::size_t rank,
        bool lastCall)
	{
   	    if ((event == WRITER_STEP_FINISHED) && (step % period != 0))
	    {
   	        return;
	    }

	    writeRegion(step, globalDimensions, localGrid, validRegion);
	    writeCount++;
	    
	    if (step == steps)
	    {
		writeRestartRegion(globalDimensions, localGrid, validRegion);
	    }
	}
    
    void createFile()
	{
	    NcmpiFile ncFile(communicator,
			     filename,
			     NcmpiFile::FileMode::replace,
			     NcmpiFile::FileFormat::classic2);

	}
    
    void createRestartFile()
	{
	    NcmpiFile ncFile(communicator,
			     restartFilename,
			     NcmpiFile::FileMode::replace,
			     NcmpiFile::FileFormat::classic2);

	}
    
    
    void writeHeader()
	{
	    NcmpiFile ncFile(communicator,
			     filename,
			     NcmpiFile::FileMode::write,
			     NcmpiFile::FileFormat::classic2);
	    
	    std::vector<std::string> ncmpiDimensionNames(4);
            // ordered according to CF convention
	    ncmpiDimensionNames[0] = "t";
	    ncmpiDimensionNames[1] = "z";
	    ncmpiDimensionNames[2] = "y";
	    ncmpiDimensionNames[3] = "x";

	    std::vector<NcmpiDim> ncmpiDimensions(DIM+1);

	    // First define time dimension, T, to have unlimited extent
	    ncmpiDimensions[0] = ncFile.addDim(ncmpiDimensionNames[0], NC_UNLIMITED);

	    // Then define the other dimensions to have extents
	    // defined by global LibGeoDecomp grid dimensions
	    for (int d = 1; d <= DIM; d++) {
		ncmpiDimensions[d] = ncFile.addDim(ncmpiDimensionNames[d+3-DIM], globalDimensions[DIM-d]);
	    }
	    
	    NcmpiVar var = ncFile.addVar(selector.name(), ncmpiDouble, ncmpiDimensions);
	    NcmpiVar timeVar = ncFile.addVar("time", ncmpiFloat, ncmpiDimensions[0]);
	    // LibGeoDecomp uses unsigned int for time step, convert implicitly to float (ncmpiFloat) for now. 
	    // If needed, can move to cdf-5 format in future, which supports unsigned ints
	    // but which is not (yet) adopted widely 
	    
	    varId = var.getId();
	    timeVarId = timeVar.getId();
	}
    
    
    void writeRestartHeader()
	{
	    NcmpiFile ncFile(communicator,
			     restartFilename,
			     NcmpiFile::FileMode::write,
			     NcmpiFile::FileFormat::classic2);
	    
	    std::vector<std::string> ncmpiDimensionNames(3);
            // ordered according to CF convention
	    ncmpiDimensionNames[0] = "z";
	    ncmpiDimensionNames[1] = "y";
	    ncmpiDimensionNames[2] = "x";
	    
	    std::vector<NcmpiDim> ncmpiDimensions(DIM);
	    
	    for (int d = 0; d < DIM; d++) {
		ncmpiDimensions[d] = ncFile.addDim(ncmpiDimensionNames[d+3-DIM], globalDimensions[DIM-d-1]);
	    }
	    
	    NcmpiVar restartVar = ncFile.addVar(selector.name(), ncmpiDouble, ncmpiDimensions);
	    restartVarId = restartVar.getId();
	}
    
    

    template<typename GRID_TYPE>
    void writeRegion(
	unsigned step,
        const Coord<DIM>& globalDimensions,
        const GRID_TYPE& localGrid,
        const Region<DIM>& region)
    {
	NcmpiFile ncFile(communicator,
			 filename,
			 NcmpiFile::FileMode::write,
			 NcmpiFile::FileFormat::classic2);

	
	NcmpiVar timeVar = NcmpiVar(ncFile, timeVarId);
	float stepFloat = step;
	
	vector<MPI_Offset> index(1); 
	index[0] = writeCount;
	timeVar.putVar_all(index, stepFloat);
	
	NcmpiVar netCDFVar = NcmpiVar(ncFile, varId);
	
	std::vector<double> buffer;
	vector<MPI_Offset> start(DIM+1);
	vector<MPI_Offset> count(DIM+1);
	
	// In 2D, this iterates over streaks of x-indexed values (each iteration is a different y)
	for (typename Region<DIM>::StreakIterator i = region.beginStreak();
	     i != region.endStreak();
             ++i) {
            // the coords need to be normalized because on torus
            // topologies the coordinates may exceed the bounding box
            // (especially negative coordinates may occurr).
            Coord<DIM> coord = Topology::normalize(i->origin, globalDimensions);
	    int xlength = i->endX - i->origin.x();
	    
	    // To do: generalise
	    start[0] = writeCount;
	    start[1] = coord.y(); // Y start index for DIM = 2, or for Z for DIM = 3
	    start[2] = 0; // X start index in 2D, = 0 for each streak iteration
	    count[0] = 1;
	    count[1] = 1; 
	    count[2] = xlength;
	    
	    std::size_t byteSize = xlength;
	    
	    if (buffer.size() != byteSize) {
                buffer.resize(byteSize);
            }
	    
            Region<DIM> tempRegion;
            tempRegion << *i;
            localGrid.saveMember(&buffer[0], MemoryLocation::HOST, selector, tempRegion);
	    
	    netCDFVar.putVar_all(start, count, &buffer[0]);
	}

	
    }

    
template<typename GRID_TYPE>
    void writeRestartRegion(
        const Coord<DIM>& globalDimensions,
        const GRID_TYPE& localGrid,
        const Region<DIM>& region)
    {
	NcmpiFile ncFile(communicator,
			 restartFilename,
			 NcmpiFile::FileMode::write,
			 NcmpiFile::FileFormat::classic2);
	
	NcmpiVar restartVar = NcmpiVar(ncFile, restartVarId);

	std::vector<double> buffer;
	vector<MPI_Offset> start(DIM);
	vector<MPI_Offset> count(DIM);
	
	// In 2D, this iterates over streaks of x-indexed values (each iteration is a different y)
	for (typename Region<DIM>::StreakIterator i = region.beginStreak();
	     i != region.endStreak();
             ++i) {
            // the coords need to be normalized because on torus
            // topologies the coordinates may exceed the bounding box
            // (especially negative coordinates may occurr).
            Coord<DIM> coord = Topology::normalize(i->origin, globalDimensions);
	    int xlength = i->endX - i->origin.x();
	    
	    // To do: generalise
	    start[0] = coord.y(); // Y start index for DIM = 2, or for Z for DIM = 3
	    start[1] = 0; // X start index in 2D, = 0 for each streak iteration
	    count[0] = 1; 
	    count[1] = xlength;
	    
	    // Uncomment to debug
	    /*std::ostringstream debug1;
	    debug1 << "NetCDF write rank" << MPILayer().rank() << ": start=" << start << ", count=" << count << std::endl;
	    std::cout << debug1.str();*/
	    
	    
	    std::size_t byteSize = xlength;
	    
	    if (buffer.size() != byteSize) {
                buffer.resize(byteSize);
            }
	    
            Region<DIM> tempRegion;
            tempRegion << *i;
            localGrid.saveMember(&buffer[0], MemoryLocation::HOST, selector, tempRegion);

	    // Uncomment to debug
	    /*std::ostringstream debug2;
	    debug2 << "NetCDF write rank" << MPILayer().rank() << " netCDF write buffer: " << buffer << std::endl;
	    std::cout << debug2.str();*/
	    
	    
	    restartVar.putVar_all(start, count, &buffer[0]);

	}
    }


    
private:
    MPIIO<CELL_TYPE, Topology> mpiio;
    Coord<DIM> globalDimensions;
    Selector<CELL_TYPE> selector;
    MPI_Comm communicator;
    std::string filename, restartFilename;
    int varId, timeVarId, restartVarId;
    unsigned steps;
    unsigned writeCount;

  
    
};
  
}

#endif
#endif


