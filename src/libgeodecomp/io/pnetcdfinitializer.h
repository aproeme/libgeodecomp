#ifndef LIBGEODECOMP_IO_PNETCDFINITIALIZER
#define LIBGEODECOMP_IO_PNETCDFINITIALIZER

#include <libgeodecomp/config.h>
#ifdef LIBGEODECOMP_WITH_MPI

#include <mpi.h>

// To do: include pnetcdf once in libgeodecomp.h
//        instead of working with include guards
#ifndef LIBGEODECOMP_INCLUDED_PNETCDF
#define LIBGEODECOMP_INCLUDED_PNETCDF
#include <pnetcdf>
#endif

#include <libgeodecomp/io/initializer.h>


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
	const unsigned steps) : 
	file(file),
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

	    
	    
	    // Assume 2D grid
	    NcmpiDim yDim = ncFile.getDim("y");
	    NcmpiDim xDim = ncFile.getDim("x");
	    
	    MPI_Offset yExtent = yDim.getSize();
	    MPI_Offset xExtent = xDim.getSize();
	    
	    dimensions = Coord<DIM>(xExtent, yExtent);
	    
	    if (MPILayer().rank() == 0) {
		std::cout << "Reading " + file << std::endl;
		std::cout << "Extent in X dimension: " << xExtent << std::endl;
		std::cout << "Extent in Y dimension: " << yExtent << std::endl;
		std::cout << "Dimensions: " + dimensions.toString() << std::endl; 
	    }
	    
	}


    Coord<DIM> gridDimensions() const
	{
	    return dimensions;
	}
    
    unsigned maxSteps() const
	{
	    return steps;
	}

    unsigned startStep() const
	{
	    return 0;
	}

    
    virtual void grid(GridBase<CELL_TYPE, DIM> *target)
	{
	    return;
	} 
    
    


private:
    std::string file;
    unsigned steps;
    Coord<DIM> dimensions;

};

}

#endif 
#endif
