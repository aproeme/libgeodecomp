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

#include <iomanip>
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
    friend class PnetCDFWriterTest;
    
    using ParallelWriter<CELL_TYPE>::period;
    using ParallelWriter<CELL_TYPE>::prefix;
    
    typedef typename ParallelWriter<CELL_TYPE>::Topology Topology;
    
    static const int DIM = Topology::DIM;
    
    PnetCDFWriter(
	const Coord<DIM>& gridDimensions,
	const Selector<CELL_TYPE>& selector,
	const std::string& prefix,
	const unsigned period,
	const MPI_Comm& communicator = MPI_COMM_WORLD) :
	Clonable<ParallelWriter<CELL_TYPE>, PnetCDFWriter<CELL_TYPE> >(prefix, period),
	gridDimensions(gridDimensions),
	selector(selector),
	comm(communicator),
	datatype(selector.mpiDatatype())
	{
	    filename = selector.name() + ".nc";
	    createFile();
	    writeHeader();
	}

    // stepFinished function signature based on that in parallelmpiiwriter.h,
    // which conforms to that specified in parallelwriter.h
    virtual void stepFinished(
        const typename ParallelWriter<CELL_TYPE>::GridType& grid,
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
	    
	    writeRegion(step, globalDimensions, grid, validRegion);
	}
    
private:
    MPIIO<CELL_TYPE, Topology> mpiio;
    Coord<DIM> gridDimensions;
    Selector<CELL_TYPE> selector;
    MPI_Comm comm;
    MPI_Datatype datatype;
    std::string filename;
    int varId, timeVarId;
    
    
    void createFile()
	{
	    NcmpiFile ncFile(comm,
			     filename,
			     NcmpiFile::FileMode::replace,
			     NcmpiFile::FileFormat::classic2);

	}
    
    // Define grid (i.e. Cell Member) time series variable in netCDF file header
    void writeHeader()
	{
	    NcmpiFile ncFile(comm,
			     filename,
			     NcmpiFile::FileMode::write,
			     NcmpiFile::FileFormat::classic2);
	    
	    std::vector<std::string> ncmpiDimensionNames(4);
            // ordered according to CF convention
	    ncmpiDimensionNames[0] = "T";
	    ncmpiDimensionNames[1] = "Z";
	    ncmpiDimensionNames[2] = "Y";
	    ncmpiDimensionNames[3] = "X";

	    std::vector<NcmpiDim> ncmpiDimensions(DIM+1);

	    // First define time dimension, T, to have unlimited extent
	    ncmpiDimensions[0] = ncFile.addDim(ncmpiDimensionNames[0], NC_UNLIMITED);

	    // Then define the other dimensions to have extents
	    // defined by LibGeoDecomp grid dimensions
	    for (int d = 1; d <= DIM; d++) {
		ncmpiDimensions[d] = ncFile.addDim(ncmpiDimensionNames[d+3-DIM], gridDimensions[DIM-d]);
	    }
	    
	    NcmpiVar var = ncFile.addVar(selector.name(), ncmpiDouble, ncmpiDimensions);
	    NcmpiVar timeVar = ncFile.addVar("time", ncmpiFloat, ncmpiDimensions[0]);
	    // LibGeoDecomp uses unsigned int for time step, convert implicitly to float (ncmpiFloat) for now. 
	    // If needed, can move to cdf-5 format in future, which supports unsigned ints
	    // but which is not (yet) adopted widely 
	    
	    varId = var.getId();
	    timeVarId = timeVar.getId();
	}


    template<typename GRID_TYPE>
    void writeRegion(
	unsigned step,
        const Coord<DIM>& dimensions,
        const GRID_TYPE& grid,
        const Region<DIM>& region)
    {
	NcmpiFile ncFile(comm,
			 filename,
			 NcmpiFile::FileMode::write,
			 NcmpiFile::FileFormat::classic2);
	
	NcmpiVar timeVar = NcmpiVar(ncFile, timeVarId);
	float stepFloat = step;
	std::cout << "step = " << step << std::endl;
	vector<MPI_Offset> index(1); 
	index[0] = stepFloat; 
	timeVar.putVar_all(index, stepFloat);
	
	
	NcmpiVar netCDFVar = NcmpiVar(ncFile, varId);
	MPI_Aint varLength = mpiio.getLength(datatype);
	std::vector<double> buffer;
	vector<MPI_Offset> start(DIM+1);
	vector<MPI_Offset> count(DIM+1);
	
	start[0] = step; 
	// GENERALISE TO DIM NOT EQUAL TO 2 (push back onto vector?)
	start[2] = 0; // X start index in 2D, = 0 for each streak iteration
	count[0] = 1; 
	count[1] = 1; 
	
	// In 2D, this iterates over streaks of x-indexed values (each iteration is a different y)
	for (typename Region<DIM>::StreakIterator i = region.beginStreak();
	     i != region.endStreak();
             ++i) {
            // the coords need to be normalized because on torus
            // topologies the coordnates may exceed the bounding box
            // (especially negative coordnates may occurr).
            Coord<DIM> coord = Topology::normalize(i->origin, dimensions);
	    int xlength = i->endX - i->origin.x();
	    
	    // GENERALISE TO DIM NOT EQUAL TO 2 (push back onto vector?)
	    start[1] = coord.y(); // Y start index for DIM = 2, or for Z for DIM = 3
	    count[2] = xlength;
	    
	    std::size_t byteSize = xlength;
	    
	    std::cout << "xlength = " + std::to_string(xlength) + '\n';
	    
            if (buffer.size() != byteSize) {
                buffer.resize(byteSize);
            }
	    
            Region<DIM> tempRegion;
            tempRegion << *i;
            grid.saveMember(&buffer[0], MemoryLocation::HOST, selector, tempRegion);

	    //
	    //  FIGURE OUT WHY WITH HIPAR SIMULATOR THE STREAK IS OF SIZE (XLENGTH) 1
	    //  WHEREAS WITH STRIPING SIMULATOR WORKS FINE (XLENGTH = entire x length of grid)
	    // 
	    
	    
	    std::string streak = "";
	    
	    //std::cout << "rank " + std::to_string(MPILayer().rank())	\
//		+ ": " + '\n' + "buffer.size()  = " + std::to_string(buffer.size()) + '\n'; 
	    
	    
	    
	    // PRINT OUT BUFFER
	    for (unsigned int index=0; index<buffer.size(); index++)
	    {
		streak += std::to_string(buffer[index]) + " ";
	    }
	    std::cout << "rank " + std::to_string(MPILayer().rank())	\
		+ ": streak = " + streak + '\n'; 
	    
	    
	    
	    /*std::cout << "rank " + std::to_string(MPILayer().rank())	\
		+ ": start[0] = step = " + std::to_string(start[0]) \
		+ ", start[1] = coord.y() = " + std::to_string(start[1]) \
		+ ", count[2] = xlength = " + std::to_string(count[2]) \
		+ ", buffer[0] = " + std::to_string(buffer[0]) \
		+ ", i->endX = " + std::to_string(i->endX) \
		+ ", i->origin.x() = " + std::to_string(i->origin.x()) \
		+ '\n';*/

	    netCDFVar.putVar_all(start, count, &buffer[0]);
	}
    }
};
  
}

#endif
#endif


