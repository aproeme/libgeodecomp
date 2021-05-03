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
	const std::string& varname,
	const unsigned steps) : 
	file(file),
	varname(varname),
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
	    
	    // To do: generalise 
	    NcmpiDim yDim = ncFile.getDim("y");
	    NcmpiDim xDim = ncFile.getDim("x");
	    
	    MPI_Offset yExtent = yDim.getSize();
	    MPI_Offset xExtent = xDim.getSize();
	    
	    dimensions = Coord<DIM>(xExtent, yExtent);
	    
	    if (MPILayer().rank() == 0) {
		std::cout << "Reading " + file << std::endl;
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
	    Region<DIM> region;
	    region << target->boundingBox();
	    readRegion(target, file, region);
	} 
    
    template<typename GRID_TYPE, int DIM>
    void readRegion(GRID_TYPE *grid,
		    const std::string& filename,
		    const Region<DIM>& region)
	{
	    NcmpiFile ncFile(MPI_COMM_WORLD,
			     file,
			     NcmpiFile::FileMode::read,
			     NcmpiFile::FileFormat::classic2);
	    
	    NcmpiVar ncVar = ncFile.getVar(varname);

	    // get the first cell in this rank's grid,
	    // determine its x and y coordinates,
	    // use these to define start point for getVar call
	    // to fetch data from netcdf file
	    CoordBox<DIM> box = grid->boundingBox();
	    
	    std::vector<MPI_Offset> start(DIM), count(DIM);
	    start[0] = box.origin.y();
	    start[1] = box.origin.x();
	    count[0] = box.dimensions.y();
	    count[1] = box.dimensions.x();
	    
	    std::ostringstream debug;
	    debug << "rank" << MPILayer().rank() << ": start=" << start << ", count=" << count << std::endl;
	    std::cout << debug.str();
	    
	    std::vector<double> buffer;
	    buffer.resize(count[0]*count[1]);
	    ncVar.getVar_all(start, count, &buffer[0]);
	    debug << "rank" << MPILayer().rank() << ": " << buffer << std::endl;
	    std::cout << debug.str();
	    
	    
	    //grid.loadMember();
	    
	    // grid->set(
	    // Work backwards 
	}
		    
    


private:
    std::string file;
    std::string varname;
    unsigned steps;
    Coord<DIM> dimensions;

};

}

#endif 
#endif
