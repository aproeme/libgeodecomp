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
 * Helper struct that defines how an application can provide
 * information specifying which netCDF source (filename, netCDF
 * variable name) is to be used to initialise a particular named grid
 * (i.e. Cell Member) variable. The PnetCDFInitializer constructor
 * above expects a vector of an arbitray number of these sources.
 */
template<typename CELL_TYPE>
struct netCDFSource
{
    std::string file;
    std::string variableName;
    Selector<CELL_TYPE> selector; // Selector describing target destination for data
};
    
    
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
	const std::vector<netCDFSource<CELL_TYPE>>& netCDFSources,
	const unsigned steps) : 
	netCDFSources(netCDFSources),
	steps(steps)
	{
	    readHeaders();
	}
    
    void readHeaders()
	{
	    std::vector<Coord<DIM>> sourceDimensions(netCDFSources.size());
	    
	    for(size_t i=0; i<netCDFSources.size(); i++)
	    {
		NcmpiFile ncFile(MPI_COMM_WORLD,
				 netCDFSources[i].file,
				 NcmpiFile::FileMode::read,
				 NcmpiFile::FileFormat::classic2);
	    
		// To do: generalise (dimensions, dimension names, including long/lat)
		NcmpiDim yDim = ncFile.getDim("y");
		NcmpiDim xDim = ncFile.getDim("x");
		
		MPI_Offset yExtent = yDim.getSize();
		MPI_Offset xExtent = xDim.getSize();
		
		sourceDimensions[i] = Coord<DIM>(xExtent, yExtent);
		
		if(i>0)
		{
		    if (sourceDimensions[i] != sourceDimensions[0])
		    {
			if (MPILayer().rank() == 0)
			{
			    std::cout << "Error: dimensions of netCDF files do not all match" << std::endl;
			}
			
			return;
		    }
		}

		// Add copy constructor to Coord class? 
		globalDimensions = Coord<2>(sourceDimensions[0][0], sourceDimensions[0][1]));
	    
		if (MPILayer().rank() == 0)
		{
		    std::cout << "Dimensions: " + globalDimensions.toString() << std::endl; 
		}
		
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

	    for(size_t i=0; i<netCDFSources.size(); i++)
	    {
		initialiseRegion(netCDFSources[i], region, localGrid);
	    }
	} 

    
    template<typename GRID_TYPE, int DIM>
    void initialiseRegion(const netCDFSource<CELL_TYPE>& ncSource,
			  const Region<DIM>& region,
			  GRID_TYPE* localGrid)
	{
	    if (MPILayer().rank() == 0)
	    {
		std::cout << "Initialising grid variable " << ncSource.selector.name()
			  <<  " from netCDF variable " + ncSource.variableName
			  << " in " << ncSource.file << std::endl;
	    }

	    NcmpiFile ncFile(MPI_COMM_WORLD,
			     ncSource.file,
			     NcmpiFile::FileMode::read,
			     NcmpiFile::FileFormat::classic2);
	    
	    NcmpiVar ncVar = ncFile.getVar(ncSource.select.variableName);

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
	    //debug1 << "NetCDF read rank" << MPILayer().rank() << ": localGrid->boundingBox() = " << localGrid->boundingBox() <<  ", start=" << start << ", count=" << count << std::endl;
	    //std::cout << debug1.str();
	    
	    std::vector<double> buffer;
	    buffer.resize(count[0]*count[1]);
	    
	    ncVar.getVar_all(start, count, &buffer[0]);
	    
	    // Uncomment to debug
	    //std::ostringstream debug2;
	    //debug2 << "NetCDF read rank" << MPILayer().rank() << " netCDF read buffer: " << buffer << std::endl;
	    //std::cout << debug2.str();
	    
	    localGrid->loadMember(&buffer[0], MemoryLocation::HOST, ncSource.selector, region);
	}
    
    
protected:
    Coord<DIM> globalDimensions;    
    unsigned steps;
    
private:
    std::vector<netCDFSource<CELL_TYPE>> netCDFSources;
    
};
    
}

#endif 
#endif
