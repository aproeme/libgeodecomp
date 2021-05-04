#ifndef LIBGEODECOMP_IO_PNETCDFINITIALIZER
#define LIBGEODECOMP_IO_PNETCDFINITIALIZER

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
#include <libgeodecomp/io/initializer.h>
#include <libgeodecomp/storage/selector.h>

using namespace PnetCDF;
using namespace PnetCDF::exceptions;

namespace LibGeoDecomp {

/**
 *
 * Initialises a specific grid variable in parallel for all grid Cells by reading
 * from a netCDF file using pnetcdf (https://github.com/Parallel-NetCDF/PnetCDF)
 * Currently using NetCDF's CDF-2 64-bit offset format (PnetCDF::NcmpiFile::FileFormat::classic2)
 * Aims to adhere to the NetCDF CF convention data model (http://cfconventions.org/)
 * 
 */


template<typename CELL_TYPE>
class PnetCDFInitializer : public Initializer<CELL_TYPE> 
{
public:
    typedef typename Initializer<CELL_TYPE>::Topology Topology;
    static const int DIM = Topology::DIM;
    
    explicit PnetCDFInitializer(
	const std::string& file,
	const Selector<CELL_TYPE>& selector,
	const std::string& ncVariableName,
	const unsigned steps) : 
	file(file),
	selector(selector),
	ncVariableName(ncVariableName),
	steps(steps)
	{
	    readHeader();
	}

    void readHeader()
	{
	    NcmpiFile ncFile(MPI_COMM_WORLD,
			     file,
			     NcmpiFile::FileMode::read,
			     NcmpiFile::FileFormat::classic2);
	    
	    // To do: generalise (dimensions, dimension names, including long/lat)
	    NcmpiDim yDim = ncFile.getDim("y");
	    NcmpiDim xDim = ncFile.getDim("x");
	    
	    MPI_Offset yExtent = yDim.getSize();
	    MPI_Offset xExtent = xDim.getSize();
	    
	   globalDimensions = Coord<DIM>(xExtent, yExtent);
	    
	    if (MPILayer().rank() == 0) {
		std::cout << "Reading " + file << std::endl;
		std::cout << "Dimensions: " + globalDimensions.toString() << std::endl; 
	    }
	}

    
    Coord<DIM> gridDimensions() const
	{
	    return globalDimensions;
	}
    
    unsigned maxSteps() const
	{
	    return steps;
	}
    
    unsigned startStep() const
	{
	    return 0;
	}
    
    virtual void grid(GridBase<CELL_TYPE, DIM> *localGrid)
	{
	    Region<DIM> region;
	    region << localGrid->boundingBox();
	    initialiseRegion(region, localGrid, file);
	} 
    
    template<typename GRID_TYPE, int DIM>
    void initialiseRegion(const Region<DIM>& region,
			  GRID_TYPE* localGrid,
			  const std::string& file)
	{
	    NcmpiFile ncFile(MPI_COMM_WORLD,
			     file,
			     NcmpiFile::FileMode::read,
			     NcmpiFile::FileFormat::classic2);
	    
	    NcmpiVar ncVar = ncFile.getVar(ncVariableName);

	    // get the first cell in this rank's grid,
	    // determine its x and y coordinates,
	    // use these to define start point and count
	    // for later getVar call that extracts data
	    // from appropriate part of netcdf file
	    CoordBox<DIM> box = localGrid->boundingBox();
	    
	    std::vector<MPI_Offset> start(DIM), count(DIM);
	    start[0] = box.origin.y();
	    start[1] = box.origin.x();
	    count[0] = box.dimensions.y();
	    count[1] = box.dimensions.x();

	    // Uncomment to debug
	    //std::ostringstream debug1;
	    //debug1 << "NetCDF read rank" << MPILayer().rank() << ": start=" << start << ", count=" << count << std::endl;
	    //std::cout << debug1.str();
	    
	    std::vector<double> buffer;
	    buffer.resize(count[0]*count[1]);

	    ncVar.getVar_all(start, count, &buffer[0]);

	    // Uncomment to debug
	    //std::ostringstream debug2;
	    //debug2 << "NetCDF read rank" << MPILayer().rank() << " netCDF read buffer: " << buffer << std::endl;
	    //std::cout << debug2.str();
	    
	    localGrid->loadMember(&buffer[0], MemoryLocation::HOST, selector, region);
	}
		    
    
protected:
    unsigned steps;

private:
    std::string file;
    std::string ncVariableName;
    Coord<DIM> globalDimensions;
    Selector<CELL_TYPE> selector;

};

}

#endif 
#endif
