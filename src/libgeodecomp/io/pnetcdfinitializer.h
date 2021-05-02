#ifndef LIBGEODECOMP_IO_PNETCDFINITIALIZER
#define LIBGEODECOMP_IO_PNETCDFINITIALIZER

#include <libgeodecomp/config.h>
#ifdef LIBGEODECOMP_WITH_MPI

#include <mpi.h>
#include <pnetcdf>

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
    typedef typename PnetCDFInitializer<CELL_TYPE>::Topology Topology;
    static const int DIM = Topology::DIM;
    
    explicit PnetCDFInitializer(
	const std::string& filename,
	const Selector<CELL_TYPE>& selector,
	const MPI_Comm& communicator = MPI_COMM_WORLD) :
	file(filename),
	selector(selector),
	comm(commmunicator)
	{
	    openFile();
	}

    void readHeader()
	{
	    NcmpiFile ncFile(comm,
			     filename,
			     NcmpiFile::FileMode::read,
			     NcmpiFile::FileFormat::classic2);
	    
	    // Assume 2D grid
	    NcmpiDim yDim = ncFile.getDim("Y");
	    NcmpiDim xDim = ncFile.getDime("X");

	    int yExtent = yDim.getSize();
	    int xExtent = xDim.getSize();
	    	    
	    if (MPILayer().rank() == 0) {
		std::cout << "Reading " + filename << std::endl;
		std::cout << "Extent in X dimension: " + xExtent << std::endl;
		std::cout << "Extent in Y dimension: " + yExtent << std::endl;
	    }
	    
	}
    
    virtual void grid(GridBase<CELL_TYPE, DIM> *target)
    {
        Region<DIM> region;
        region << target->boundingBox();
        //mpiiop.readRegion(target, file, region, communicator, datatype);
    }

    


private:
    std::string file;
    Selector<CELL_TYPE> selector;
    MPI_Comm comm;
	
	

}

#endif 
