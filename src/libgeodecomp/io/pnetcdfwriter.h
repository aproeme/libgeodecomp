#ifndef LIBGEODECOMP_IO_PNETCDFWRITER_H
#define LIBGEODECOMP_IO_PNETCDFWRITER_H

#include <libgeodecomp/config.h>
#ifdef LIBGEODECOMP_WITH_MPI

#include <iostream>
#include <mpi.h>
#include <pnetcdf>

#include <libgeodecomp/communication/mpilayer.h>
#include <libgeodecomp/io/parallelwriter.h>
#include <libgeodecomp/misc/clonable.h>

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
	Clonable<PnetCDFWriter<CELL_TYPE>, PnetCDFWriter<CELL_TYPE> >(prefix, period),
	gridDimensions(gridDimensions),
	selector(selector),
	comm(communicator)
	{
	    filename = selector.name() + ".nc";
	    createFile();
	    writeHeader();
	}

    
private:
    Coord<DIM>& gridDimensions;
    Selector<CELL_TYPE> selector;
    MPI_Comm comm;
    std::string filename;
    
    void createFile()
	{
	    try
	    {
		NcmpiFile ncFile(comm,
				 filename,
				 NcmpiFile::FileMode::newFile,
				 NcmpiFile::FileFormat::classic2);
	    }
	    catch(NcmpiException& e)
	    {
		std::cout << "PnetCDF" << e.what() << " error code=" << e.errorCode() << " Error!\n";
	    }
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
		
		std::vector<NcmpiDim> ncmpiDimensions(DIM);
		ncmpiDimensions[0] = ncFile.addDim(ncmpiDimensionNames[0], NC_UNLIMITED);
		
		for (int d = DIM-1; d > 0; d--) {
		    ncmpiDimensions[d] = ncFile.addDim(ncmpiDimensionNames[d], gridDimensions[(DIM-1)-d]);
		}
	  
		ncFile.addVar(selector.name(), ncmpiDouble, ncmpiDimensions);
	    }
	    catch(NcmpiException& e)
	    {
		std::cout << "PnetCDF" << e.what() << " error code=" << e.errorCode() << " Error!\n";
	    }
	}
    
};
  
}

#endif
#endif
