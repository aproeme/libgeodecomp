#ifndef _libgeodecomp_misc_logger_h_
#define _libgeodecomp_misc_logger_h_

#include "boost/date_time/posix_time/posix_time.hpp"

// fixme: configure this via "configure"
#define LIBGEODECOMP_FEATURE_DEBUG_LEVEL 0
namespace LibGeoDecomp {

/**
 * Logger provides a set of functions and macros to selectively print
 * different amounts/levels of debugging information. Filtering
 * happens at compile time, thus no runtime overhead is incurred for
 * filtered messages.
 */
class Logger {
public:
    static const int FATAL = 0;
    static const int ERROR = 1;
    static const int WARN  = 2;
    static const int INFO  = 3;
    static const int DEBUG = 4;
};

}

#if LIBGEODECOMP_FEATURE_DEBUG_LEVEL < 0
#define LOG(LEVEL, MESSAGE) 
#endif

#if LIBGEODECOMP_FEATURE_DEBUG_LEVEL >= 0
#define LOG(LEVEL, MESSAGE)                                             \
    if (LibGeoDecomp::Logger::LEVEL <= LIBGEODECOMP_FEATURE_DEBUG_LEVEL) { \
        std::cout << #LEVEL[0] << ", ["                                 \
                  << boost::posix_time::to_iso_string(boost::posix_time::second_clock::local_time()) \
                  << "] " << std::right                                 \
                  << std::setw(5) << #LEVEL                             \
                  << " -- " << MESSAGE << "\n";                         \
    }
#endif

#endif
