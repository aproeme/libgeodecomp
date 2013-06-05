#ifndef LIBGEODECOMP_IO_REMOTESTEERER_PIPE_H
#define LIBGEODECOMP_IO_REMOTESTEERER_PIPE_H

#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <libgeodecomp/config.h>
#include <libgeodecomp/misc/stringops.h>
#include <libgeodecomp/mpilayer/mpilayer.h>

namespace LibGeoDecomp {

namespace RemoteSteererHelpers {

/**
 * The purpose of the Pipe is to forward steering data between the
 * connection node and the execution nodes. This forwarding is done
 * asynchronously as Actions from the connection node run in a
 * seperate thread while Handlers are run by the RemoteSteerer, which
 * in turn is triggered by the Simulator. The Pipe also manages
 * scattering/gathering of the data via MPI, which is essentially a
 * 1-n operation. The external interface is thread-safe.
 */
class Pipe
{
public:
    friend class PipeTest;
    typedef StringOps::StringVec StringVec;

#ifdef LIBGEODECOMP_FEATURE_MPI
    Pipe(
        int root = 0,
        MPI::Comm *communicator = &MPI::COMM_WORLD) :
        mpiLayer(communicator),
        root(root)
    {}
#endif

    void addSteeringRequest(const std::string& request)
    {
        boost::lock_guard<boost::mutex> lock(mutex);
        steeringRequests << request;
    }

    void addSteeringFeedback(const std::string& feedback)
    {
        boost::lock_guard<boost::mutex> lock(mutex);
        steeringFeedback << feedback;
    }

    StringVec retrieveSteeringRequests()
    {
        boost::lock_guard<boost::mutex> lock(mutex);
        StringVec requests;
        std::swap(requests, steeringRequests);
        return requests;
    }

    StringVec retrieveSteeringFeedback()
    {
        boost::lock_guard<boost::mutex> lock(mutex);
        StringVec feedback;
        std::swap(feedback, steeringFeedback);
        return feedback;
    }

    void sync()
    {
        boost::lock_guard<boost::mutex> lock(mutex);

#ifdef LIBGEODECOMP_FEATURE_MPI
        broadcastSteeringRequests();
        moveSteeringFeedbackToRoot();
#endif
    }

private:
    boost::mutex mutex;
    StringVec steeringRequests;
    StringVec steeringFeedback;

#ifdef LIBGEODECOMP_FEATURE_MPI
    MPILayer mpiLayer;
    int root;

    void broadcastSteeringRequests()
    {
        int numRequests = mpiLayer.broadcast(steeringRequests.size(), 0);
        if (mpiLayer.rank() != root) {
            steeringRequests.resize(numRequests);
        }

        SuperVector<int> requestSizes(numRequests);

        if (mpiLayer.rank() == root) {
            for (int i = 0; i < numRequests; ++i) {
                requestSizes[i] = steeringRequests[i].size();
            }
        }

        mpiLayer.broadcastVector(&requestSizes, root);
        if (mpiLayer.rank() != root) {
            for (int i = 0; i < numRequests; ++i) {
                steeringRequests[i].resize(requestSizes[i]);
            }
        }

        for (int i = 0; i < numRequests; ++i) {
            mpiLayer.broadcast(&steeringRequests[i][0], requestSizes[i], root);
        }
    }

    /**
     * Will move steering requests from the compute nodes to the root
     * where it can then be forwarded to the user.
     */
    void moveSteeringFeedbackToRoot()
    {
        // marshall all feedback in order to use scalable gather afterwards
        SuperVector<char> localBuffer;
        SuperVector<int> localLengths;
        for (StringVec::iterator i = steeringFeedback.begin();
             i != steeringFeedback.end();
             ++i) {
            localBuffer.insert(localBuffer.end(), i->begin(), i->end());
            localLengths << i->size();
        }

        // how many strings are sent per node?
        SuperVector<int> numFeedback = mpiLayer.gather((int)localLengths.size(), root);

        // all lengths of all strings:
        SuperVector<int> allFeedbackLengths(numFeedback.sum());
        mpiLayer.gatherV(&localLengths[0], localLengths.size(), numFeedback, root, &allFeedbackLengths[0]);

        // gather all messages in a single, giant buffer:
        SuperVector<int> charsPerNode = mpiLayer.gather((int)localBuffer.size(), root);
        SuperVector<char> globalBuffer(charsPerNode.sum());
        mpiLayer.gatherV(&localBuffer[0], localBuffer.size(), charsPerNode, root, &globalBuffer[0]);


        // reconstruct strings:
        int cursor = 0;
        steeringFeedback.resize(0);

        for (SuperVector<int>::iterator i = allFeedbackLengths.begin();
             i != allFeedbackLengths.end();
             ++i) {
            int nextCursor = cursor + *i;
            steeringFeedback << std::string(&globalBuffer[cursor],
                                            &globalBuffer[nextCursor]);
            cursor = nextCursor;
        }
    }
#endif
};

}

}

#endif
