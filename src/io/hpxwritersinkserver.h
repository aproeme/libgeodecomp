#ifndef LIBGEODECOMP_IO_HPXWRITERSINKSERVER_H
#define LIBGEODECOMP_IO_HPXWRITERSINKSERVER_H

#include <libgeodecomp/config.h>
#ifdef LIBGEODECOMP_WITH_HPX_XXX

#include <libgeodecomp/io/parallelwriter.h>
#include <libgeodecomp/io/writer.h>
#include <libgeodecomp/storage/gridvecconv.h>

#include <hpx/hpx_fwd.hpp>
#include <hpx/include/components.hpp>
#include <hpx/lcos/local/spinlock.hpp>

namespace LibGeoDecomp {

/**
 * Forward declaration
 */
template<typename CELL_TYPE>
class DistributedSimulator;

/**
 * HPX component to complement HpxWriterSink
 */
template<typename CELL_TYPE>
class HpxWriterSinkServer : public hpx::components::managed_component_base<HpxWriterSinkServer<CELL_TYPE> >
{
private:
    class RegionInfo;
public:
    typedef CELL_TYPE CellType;
    typedef typename APITraits::SelectTopology<CellType>::Value Topology;
    static const int DIM = Topology::DIM;
    typedef DisplacedGrid<CellType, Topology> GridType;
    typedef Region<Topology::DIM> RegionType;
    typedef Coord<Topology::DIM> CoordType;
    typedef std::vector<CellType> BufferType;
    typedef std::map<unsigned, std::vector<RegionInfo> > RegionInfoMapType;
    typedef std::map<unsigned, std::size_t> StepCountMapType;
    typedef std::map<unsigned, GridType> GridMapType;
    typedef
        std::map<std::size_t, boost::shared_ptr<ParallelWriter<CellType> > >
        ParallelWritersMap;
    typedef
        std::map<std::size_t, boost::shared_ptr<Writer<CellType> > >
        SerialWritersMap;

    typedef hpx::lcos::local::spinlock MutexType;

    HpxWriterSinkServer(
        std::size_t numUpdateGroups = 0) :
        numUpdateGroups(numUpdateGroups),
        nextId(0)
    {}

    HpxWriterSinkServer(
        boost::shared_ptr<ParallelWriter<CellType> > parallelWriter,
        std::size_t numUpdateGroups) :
        numUpdateGroups(numUpdateGroups),
        nextId(0)
    {
        connectParallelWriter(parallelWriter);
    }

    HpxWriterSinkServer(
        boost::shared_ptr<Writer<CellType> > serialWriter,
        std::size_t numUpdateGroups) :
        numUpdateGroups(numUpdateGroups),
        nextId(0)
    {
        connectSerialWriter(serialWriter);
    }

    void stepFinished(
        boost::shared_ptr<BufferType> buffer,
        const RegionType& validRegion,
        const CoordType& globalDimensions,
        unsigned step,
        WriterEvent event,
        std::size_t rank,
        bool lastCall)
    {
        // typedef typename RegionInfoMapType::iterator RegionMapIterator;
        // typedef typename StepCountMapType::iterator StepCountMapIterator;
        // typedef typename GridMapType::iterator GridIterator;

        // MutexType::scoped_lock l(mutex);
        // GridIterator kt = gridMap.find(step);
        // if (kt == gridMap.end())
        // {
        //     CoordBox<DIM> coordBox(CoordType(), globalDimensions);
        //     kt = gridMap.insert(
        //         kt,
        //         std::make_pair(
        //             step,
        //             GridType(coordBox)));
        // }

        // GridVecConv::vectorToGrid(*buffer, &kt->second, validRegion);

        // RegionMapIterator it = regionInfoMap.find(step);
        // if (it == regionInfoMap.end()) {
        //     it = regionInfoMap.insert(
        //             it,
        //             std::make_pair(step, std::vector<RegionInfo>())
        //         );
        // }

        // it->second.push_back(
        //     RegionInfo(validRegion, globalDimensions, event, rank, lastCall)
        // );

        // if (lastCall) {
        //     StepCountMapIterator jt = stepCountMap.find(step);
        //     if (jt == stepCountMap.end())
        //     {
        //         jt = stepCountMap.insert(jt, std::make_pair(step, 1));
        //     }
        //     else
        //     {
        //         ++jt->second;
        //     }

        //     if (jt->second == numUpdateGroups)
        //     {
        //         {
        //             hpx::util::scoped_unlock<MutexType::scoped_lock> ull(l);
        //             notifyWriters(kt->second, step, event);
        //         }
        //         regionInfoMap.erase(it);
        //         stepCountMap.erase(jt);
        //         gridMap.erase(kt);
        //     }
        // }
    }
    HPX_DEFINE_COMPONENT_ACTION(HpxWriterSinkServer, stepFinished, StepFinishedAction);

    std::size_t connectParallelWriter(
        boost::shared_ptr<ParallelWriter<CellType> > parallelWriter)
    {
        MutexType::scoped_lock l(mutex);
        std::size_t id = getNextId();
        parallelWriters.insert(std::make_pair(id, parallelWriter));
        return id;
    }
    HPX_DEFINE_COMPONENT_ACTION(HpxWriterSinkServer, connectParallelWriter, ConnectParallelWriterAction);

    std::size_t connectSerialWriter(
        boost::shared_ptr<Writer<CellType> > serialWriter)
    {
        MutexType::scoped_lock l(mutex);
        std::size_t id = getNextId();
        serialWriters.insert(std::make_pair(id, serialWriter));
        return id;
    }
    HPX_DEFINE_COMPONENT_ACTION(HpxWriterSinkServer, connectSerialWriter, ConnectSerialWriterAction);

    void disconnectWriter(std::size_t id)
    {
        MutexType::scoped_lock l(mutex);
        {
            typename ParallelWritersMap::iterator it = parallelWriters.find(id);
            if (it != parallelWriters.end()) {
                parallelWriters.erase(it);
                freeIds.push_back(id);
                return;
            }
        }
        {
            typename SerialWritersMap::iterator it = serialWriters.find(id);
            if (it != serialWriters.end()) {
                serialWriters.erase(it);
                freeIds.push_back(id);
                return;
            }
        }
    }
    HPX_DEFINE_COMPONENT_ACTION(HpxWriterSinkServer, disconnectWriter, DisconnectWriterAction);

    std::size_t getNumUpdateGroups()
    {
        return numUpdateGroups;
    }
    HPX_DEFINE_COMPONENT_ACTION(HpxWriterSinkServer, getNumUpdateGroups, NumUpdateGroupsAction);

private:
    GridMapType gridMap;
    ParallelWritersMap parallelWriters;
    SerialWritersMap serialWriters;
    std::size_t numUpdateGroups;

    StepCountMapType stepCountMap;
    RegionInfoMapType regionInfoMap;

    std::size_t nextId;
    std::vector<std::size_t> freeIds;

    MutexType mutex;

    std::size_t getNextId()
    {
        std::size_t id = 0;
        if (!freeIds.empty()) {
            id = freeIds.back();
            freeIds.pop_back();
        }
        else {
            id = nextId++;
        }
        return id;
    }

    void notifyWriters(GridType const & grid, unsigned step, WriterEvent event)
    {
        BOOST_FOREACH(typename ParallelWritersMap::value_type& writer,
                      parallelWriters) {
            MutexType::scoped_lock l(mutex);
            typedef typename RegionInfoMapType::iterator RegionInfoIterator;

            RegionInfoIterator it = regionInfoMap.find(step);
            BOOST_ASSERT(it != regionInfoMap.end());
            BOOST_FOREACH(RegionInfo const & regionInfo, it->second) {
                writer.second->stepFinished(
                    grid,
                    regionInfo.validRegion,
                    regionInfo.globalDimensions,
                    step,
                    regionInfo.event, regionInfo.rank,
                    regionInfo.lastCall
                );
            }
        }
        BOOST_FOREACH(typename SerialWritersMap::value_type& writer,
                      serialWriters) {
            writer.second->stepFinished(grid, step, event);
        }
    }

    class RegionInfo
    {
    public:
        RegionInfo(
            const RegionType& validRegion,
            const CoordType& globalDimensions,
            WriterEvent event,
            std::size_t rank,
            bool lastCall) :
            validRegion(validRegion),
            globalDimensions(globalDimensions),
            event(event),
            rank(rank),
            lastCall(lastCall)
        {}

        RegionType validRegion;
        CoordType globalDimensions;
        WriterEvent event;
        std::size_t rank;
        bool lastCall;
    };

};

}

#endif
#endif
