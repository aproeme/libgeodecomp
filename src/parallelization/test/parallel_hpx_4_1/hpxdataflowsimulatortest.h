#include <cxxtest/TestSuite.h>
#include <hpx/hpx.hpp>

#include <libgeodecomp/io/initializer.h>
#include <libgeodecomp/misc/apitraits.h>
#include <libgeodecomp/storage/unstructuredgrid.h>
#include <libgeodecomp/parallelization/hpxdataflowsimulator.h>

using namespace LibGeoDecomp;

namespace LibGeoDecomp {

class DummyMessage
{
public:
    DummyMessage(int senderId = -1,
                 int receiverId = -1,
                 int timestep = -1,
                 int data = -1) :
        senderId(senderId),
        receiverId(receiverId),
        timestep(timestep),
        data(data)
    {}

    template<typename ARCHIVE>
    void serialize(ARCHIVE& archive, int)
    {
        archive & senderId;
        archive & receiverId;
        archive & timestep;
        archive & data;
    }

    int senderId;
    int receiverId;
    int timestep;
    int data;
};

}

LIBGEODECOMP_REGISTER_HPX_COMM_TYPE(DummyMessage)

namespace LibGeoDecomp {

class DummyModel
{
public:
    class API :
        public APITraits::HasUnstructuredTopology
    {};

    DummyModel(int id = -1) :
        id(id)
    {}

    // fixme: use move semantics here
    int update(
        std::vector<hpx::shared_future<DummyMessage> > inputFutures,
        const hpx::shared_future<int>& /* unused */,
        int step,
        std::map<int, hpx::id_type> remoteIDs)
    {
	std::cout << "updating Dummy " << id << " and my neighbors are: [";
        std::vector<DummyMessage> input = hpx::util::unwrapped(inputFutures);

	for (int i = 0; i != input.size(); ++i) {
	    std::cout << input[i].senderId << " -> " << input[i].receiverId;
	}
	std::cout << "]\n";

        for (auto&& iter: remoteIDs) {
            DummyMessage dummyMessage;
            std::cout << "storing ID " << iter.second << ", step: " << step << "\n";
            hpx::apply(typename HPXReceiver<DummyMessage>::receiveAction(), iter.second,  step, dummyMessage);
        }

        return 0;
    }

private:
    int id;
};

 class DummyInitializer : public Initializer<DummyModel>
 {
 public:
     DummyInitializer(int gridSize, int myMaxSteps) :
         gridSize(gridSize),
	 myMaxSteps(myMaxSteps)
     {}

     void grid(GridBase<DummyModel, 1> *grid)
     {
	 CoordBox<1> box = grid->boundingBox();

	 for (CoordBox<1>::Iterator i = box.begin(); i != box.end(); ++i) {
	     grid->set(*i, DummyModel(i->x()));
	 }
     }

    virtual Coord<1> gridDimensions() const
    {
	return Coord<1>(gridSize);
    }

    unsigned startStep() const
    {
	return 0;
    }

    unsigned maxSteps() const
    {
	return myMaxSteps;
    }

    boost::shared_ptr<Adjacency> getAdjacency(const Region<1>& region) const
    {
	boost::shared_ptr<Adjacency> adjacency(new RegionBasedAdjacency());

	for (Region<1>::Iterator i = region.begin(); i != region.end(); ++i) {
            // fixme: have more connections
	    if (i->x() != 0) {
		adjacency->insert(i->x(), i->x() - 1);
	    }

	    if (i->x() != (gridSize - 1)) {
		adjacency->insert(i->x(), i->x() + 1);
	    }
	}

	return adjacency;
    }

 private:
     int gridSize;
     int myMaxSteps;
};

template<typename CELL>
class TestComponent
{
public:
    explicit TestComponent(CELL *cell = 0) :
        cell(cell)
    {}

    CELL *cell;
    std::map<int, std::shared_ptr<HPXReceiver<DummyMessage> > > receivers;
    std::map<int, hpx::id_type> remoteIDs;
};

class HpxDataflowSimulatorTest : public CxxTest::TestSuite
{
public:
    void setUp()
    {
    }

    void testBasic()
    {
        std::cout << "starting dataflow test\n";

        DummyInitializer initializer(50, 13);
        UnstructuredGrid<DummyModel> grid(initializer.gridBox());
        initializer.grid(&grid);
        std::cout << "grid size: " << grid.boundingBox() << "\n";

        typedef hpx::shared_future<int> UpdateResultFuture;
        typedef std::map<int, UpdateResultFuture> TimeStepFutures;

        using hpx::dataflow;
        using hpx::util::unwrapped;
        TimeStepFutures lastTimeStepFutures;
        TimeStepFutures thisTimeStepFutures;
        std::cout << "blah0\n";

        Region<1> localRegion;
        CoordBox<1> box = initializer.gridBox();
        int rank = hpx::get_locality_id();
        for (int i = ((rank + 0) * box.dimensions.x() / 4);
             i <     ((rank + 1) * box.dimensions.x() / 4);
             ++i) {
            localRegion << Coord<1>(i);
        }

        boost::shared_ptr<Adjacency> adjacency = initializer.getAdjacency(localRegion);
        std::cout << "blah1\n";

        // fixme: instantiate components in agas and only hold ids of those
        std::map<int, TestComponent<DummyModel> > components;
        for (Region<1>::Iterator i = localRegion.begin(); i != localRegion.end(); ++i) {
            std::cout << "   blubbA " << *i << "\n";
            TestComponent<DummyModel> component(&grid[*i]);

            std::cout << "   blubbB " << *i << "\n";
            std::vector<int> neighbors;
            adjacency->getNeighbors(i->x(), &neighbors);
            std::cout << "   blubbC " << neighbors << "\n";

            for (auto j = neighbors.begin(); j != neighbors.end(); ++j) {
                component.receivers[*j] = HPXReceiver<DummyMessage>::make(
                    "hpx_receiver_" +
                    StringOps::itoa(*j) +
                    "_to_" +
                    StringOps::itoa(i->x())).get();
            }

            std::cout << "   blubbD " << *i << "\n";
            components[i->x()] = component;
            std::cout << "   blubbE " << *i << "\n";
        }

        for (Region<1>::Iterator i = localRegion.begin(); i != localRegion.end(); ++i) {
            TestComponent<DummyModel>& component = components[i->x()];

            std::vector<int> neighbors;
            adjacency->getNeighbors(i->x(), &neighbors);

            for (auto j = neighbors.begin(); j != neighbors.end(); ++j) {
                std::string linkName = "hpx_receiver_" +
                    StringOps::itoa(i->x()) +
                    "_to_" +
                    StringOps::itoa(*j);

                component.remoteIDs[*j] = hpx::id_type(HPXReceiver<DummyMessage>::find(linkName).get());
            }
        }

        std::cout << "setting up dataflow\n";
        for (Region<1>::Iterator i = localRegion.begin(); i != localRegion.end(); ++i) {
            lastTimeStepFutures[i->x()] = hpx::make_ready_future(UpdateResultFuture());
        }

        std::cout << "blah2\n";
        int maxTimeSteps = initializer.maxSteps();
        for (int t = 0; t < maxTimeSteps; ++t) {

            for (Region<1>::Iterator i = localRegion.begin(); i != localRegion.end(); ++i) {

                std::vector<hpx::shared_future<DummyMessage> > receiveMessagesFutures;
                std::vector<int> neighbors;
                adjacency->getNeighbors(i->x(), &neighbors);
                for (auto j = neighbors.begin(); j != neighbors.end(); ++j) {
                    receiveMessagesFutures << hpx::make_ready_future(DummyMessage(-1, -1, -1, -1));
                    // receiveMessagesFutures << components[i->x()].receivers[*j].get(t);
                }

                auto Operation = boost::bind(&DummyModel::update, *components[i->x()].cell, _1, _2, _3, _4);

                thisTimeStepFutures[i->x()] = dataflow(
                    hpx::launch::async,
                    Operation,
                    receiveMessagesFutures,
                    lastTimeStepFutures[i->x()],
                    t,
                    components[i->x()].remoteIDs);

                // for (Region<1>::Iterator i = localRegion.begin(); i != localRegion.end(); ++i) {

                //     // fixme: use hpxreceiver::receive to get futures
                //     auto Op = unwrapped(boost::bind(&DummyModel::update, &grid[i->x()], _1));
                //     std::vector<ResultsFuture> localDependencies;

                //     std::vector<int> neighbors;
                //     adjacency->getNeighbors(i->x(), &neighbors);

                //     for (std::vector<int>::iterator n = neighbors.begin(); n != neighbors.end(); ++n) {
                //         localDependencies.push_back(lastTimeStepFutures[*n]);
                //     }

                //     thisTimeStepFutures[i->x()] = dataflow(hpx::launch::async, Op, localDependencies);
            }

            using std::swap;
            swap(thisTimeStepFutures, lastTimeStepFutures);
        }

        std::cout << "waiting on futures\n";

        // std::vector<ResultsFuture> finalStep;
        // for (Region<1>::Iterator i = localRegion.begin(); i != localRegion.end(); ++i) {
        //     finalStep.push_back(futures[maxTimeSteps % 2][i->x()]);
        // }

        // fixme: this is ugly
        for (auto&& i: lastTimeStepFutures) {
            i.second.get();
        }
        // hpx::when_all(lastTimeStepFutures).get();

        std::cout << "dataflow test done\n";
    }
};

}
