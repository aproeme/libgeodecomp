#ifndef LIBGEODECOMP_IO_PNETCDFWRITER_HA
#define LIBGEODECOMP_IO_PNETCDFWRITER_H

#include <libgeodecomp/config.h>
#ifdef LIBGEODECOMP_WITH_MPI

#include <iostream>
#include <mpi.h>
#include <pnetcdf>

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
 * Writes simulation variables to netCDF's CDF-2 (64-bit offset) format
 * in parallel using pnetcdf. Aims to use CF conventions (http://cfconventions.org/)
 *
 * Needs overall grid size in order to define dimensions of grid variables in
 * netCDF file header prior to simulation. 
 *
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
	comm(communicator)
	{
	    filename = selector.name() + ".nc";
	    createFile();
	    writeHeader();
	}
    
    
    virtual void stepFinished(
        const typename ParallelWriter<CELL_TYPE>::GridType& grid,
        const Region<Topology::DIM>& validRegion,
        const Coord<Topology::DIM>& globalGridDimensions,
        unsigned step,
        WriterEvent event,
        std::size_t rank,
        bool lastCall)
	{}
    
    
    
private:
    Coord<DIM> gridDimensions;
    Selector<CELL_TYPE> selector;
    MPI_Comm comm;
    std::string filename;
    
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
	    try
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
		ncmpiDimensions[0] = ncFile.addDim(ncmpiDimensionNames[0], NC_UNLIMITED);
		
		for (int d = 1; d <= DIM; d++) {
		    ncmpiDimensions[d] = ncFile.addDim(ncmpiDimensionNames[d+3-DIM], gridDimensions[DIM-d]);
		}
	  
		ncFile.addVar(selector.name(), ncmpiDouble, ncmpiDimensions);
	    }
	    catch(NcmpiException& e)
	    {
		std::cout << "PnetCDF" << e.what()  << " error code=" << e.errorCode() << " Error!\n";
	    }
	}
    
};
  
}

#endif
#endif
