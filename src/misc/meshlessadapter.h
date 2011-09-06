#ifndef _libgeodecomp_misc_meshlessadapter_h_
#define _libgeodecomp_misc_meshlessadapter_h_

#include <list>
#include <set>
#include <libgeodecomp/misc/floatcoord.h>
#include <libgeodecomp/misc/grid.h>
#include <libgeodecomp/misc/supermap.h>
#include <libgeodecomp/misc/topologies.h>

namespace LibGeoDecomp {

/**
 * A utility class which supports users in porting meshless codes to
 * LibGeoDecomp by superimposing a stencil structure, even though the
 * actual cells may be connected by an irregular graph.
 */
template<class TOPOLOGY=Topologies::Torus<2>::Topology>
class MeshlessAdapter 
{
public:
    friend class MeshlessAdapterTest;
    static const int DIM = TOPOLOGY::DIMENSIONS;
    static const int MAX_SIZE = 100000;

    typedef std::list<FloatCoord<DIM> > CoordList;
    typedef Grid<CoordList, TOPOLOGY> CoordListGrid;
    typedef SuperVector<FloatCoord<DIM> > CoordVec;
    typedef SuperVector<SuperVector<int> > Graph;

    /**
     * creates an MeshlessAdapter which assumes that the coordinates
     * of the cells are elementwise smaller than _dimensions. The
     * stencil cells will be of size boxSize^DIMENSIONS.
     */
    inline MeshlessAdapter(
        const FloatCoord<DIM>& _dimensions=FloatCoord<DIM>(), 
        const double& _boxSize=1) :
        dimensions(_dimensions)
    {
        resetBoxSize(_boxSize);
    }
    
    inline CoordListGrid grid() const
    {
        return CoordListGrid(discreteDim);
    }
    
    inline Coord<DIM> posToCoord(const FloatCoord<DIM>& pos) const
    {
        Coord<DIM> c;
        for (int i = 0; i < DIM; ++i) {
            // cut all overhanging coords
            c.c[i] = std::min(discreteDim.c[i] - 1,
                              (int)(pos.c[i] * scale));
        }
        return c;
    }

    inline void insert(CoordListGrid *grid, const FloatCoord<DIM>& pos) const
    {
        Coord<2> c = posToCoord(pos);
        (*grid)[c].push_back(pos);
    }

    /**
     * checks if the grid cell containing pos or any of its neighbors
     * in its Moore neighborhood contains a vertex which is closer to
     * pos than the boxSize. May return a list of all found vertex IDs
     * if coords is set.
     */
    bool search(
        const CoordListGrid& positions,
        const FloatCoord<DIM>& pos,
        std::set<int> *coords = 0) const
    {
        bool found = false;
        Coord<DIM> center = posToCoord(pos);
   
        CoordBox<DIM> box(CoordDiagonal<DIM>()(-1), CoordDiagonal<DIM>()(3));
        CoordBoxSequence<DIM> s = box.sequence();
        while (s.hasNext()) {
            Coord<DIM> newCenter = center + s.next();
            bool res = searchList(positions[newCenter], pos, coords);
            found |= res;
        }

        return found;
    }

    inline CoordVec findAllPositions(const CoordListGrid& positions) const
    {
        CoordVec ret;

        CoordBox<DIM> box = positions.boundingBox();
        CoordBoxSequence<DIM> s = box.sequence();
        while (s.hasNext()) {
            const std::list<FloatCoord<DIM> >& list = positions[s.next()];
            for (typename std::list<FloatCoord<DIM> >::const_iterator iter = list.begin();
                 iter != list.end();
                 ++iter) {
                ret.push_back(*iter);
            }
        }

        return ret;
    }

    double findOptimumBoxSize(
        const CoordVec& positions, 
        const Graph& graph)
    {
        double upperBorder = -1;
        double lowerBorder = -1;

        Coord<DIM> upperBorderDim;
        Coord<DIM> lowerBorderDim;
        
        // the exponential back-off algorithm allows us to find upper
        // and lower bound for the optimal box size in log time
        if (checkBoxSize(positions, graph)) {
            upperBorder = boxSize;
            upperBorderDim = discreteDim;
            resetBoxSize(boxSize * 0.5);
            while (checkBoxSize(positions, graph))
                resetBoxSize(boxSize * 0.5);
            lowerBorder = boxSize;
            lowerBorderDim = discreteDim;
        } else {
            lowerBorder = boxSize;
            lowerBorderDim = discreteDim;
            resetBoxSize(boxSize * 2);
            while (!checkBoxSize(positions, graph))
                resetBoxSize(boxSize * 2);
            upperBorder = boxSize;
            upperBorderDim = discreteDim;
        }

        // The loop condition is not very tight. I'd rather test for
        // equality here. But the loose test is necessary to keep the
        // number of iterations low in the general case, and finite in
        // special cases where deadlocks would occurr otherwise.
        while (intervalTooLarge(lowerBorderDim, upperBorderDim)) {
            double middle = (upperBorder + lowerBorder) * 0.5;
            resetBoxSize(middle);
            if (checkBoxSize(positions, graph)) {
                upperBorder = middle;
                upperBorderDim = discreteDim;
            } else {
                lowerBorder = middle;
                lowerBorderDim = discreteDim;
            }
        }
      
        double maxBoxSize = 0;
        double nextLower = 0;
        
        // the resulting box size may lead to overhangs on the far
        // boundaries which would get added to the nearest container
        // cells. this step enlarges the box size slightly to ensure a
        // smooth distribution.
        for (int i = 0; i < DIM; ++i) {
            double current = dimensions.c[i] / upperBorderDim.c[i];
            if (current > maxBoxSize) {
                maxBoxSize = current;
                // because the loop condition above is not tight, the
                // next smaller box size might represet a valid
                // solution, too.
                nextLower = dimensions.c[i] / (upperBorderDim.c[i] + 1);
            }
        }
            
        resetBoxSize(nextLower);
        if (checkBoxSize(positions, graph)) 
            return nextLower;

        resetBoxSize(maxBoxSize);
        // this should never happen, but might in arcane cases where
        // the packaging density is much higher on the far boundaries
        // than elsewhere in the the simulation space.
        if (!checkBoxSize(positions, graph))
            throw std::logic_error("failed to determine a valid box size");
        return maxBoxSize;
    }

    bool checkBoxSize(const CoordVec& positions, const Graph& graph)
    {
        for (int i = 0; i < graph.size(); ++i) 
            for (SuperVector<int>::const_iterator n = graph[i].begin(); 
                 n != graph[i].end(); ++n) 
                if (manhattanDistance(positions[i], positions[*n]) > 1)
                    return false;

        return true;
    }

    const Coord<DIM>& getDiscreteDim() const
    {
        return discreteDim;
    }

    SuperMap<std::string, double> reportFillLevels(const CoordVec& positions) const
    {
        SuperMap<Coord<DIM>, int> cache;
        for (typename CoordVec::const_iterator i = positions.begin(); i != positions.end(); ++i) {
            Coord<DIM> c = posToCoord(*i);
            cache[c]++;
        }

        long sum = 0;
        long emptyCells = 0;
        int lowestFill = cache[Coord<DIM>()];
        int highestFill = cache[Coord<DIM>()];

        CoordBox<DIM> box(Coord<DIM>(), discreteDim);
        CoordBoxSequence<DIM> s = box.sequence();
        while (s.hasNext()) {
            Coord<DIM> c = s.next();
            lowestFill  = std::min(cache[c], lowestFill);
            highestFill = std::max(cache[c], highestFill);
            sum += cache[c];

            if (cache[c] == 0)
                ++emptyCells;
        }

        SuperMap<std::string, double> ret;
        ret["emptyCells"]  = emptyCells;
        ret["averageFill"] = 1.0 * sum / discreteDim.prod();
        ret["lowestFill"]  = lowestFill;
        ret["highestFill"] = highestFill;
        return ret;
    }

    const double& getBoxSize() const
    {
        return boxSize;
    }

private:
    FloatCoord<DIM> dimensions;
    Coord<DIM> discreteDim;
    double scale;
    double radius2;
    double boxSize;

    void resetBoxSize(const double& newBoxSize)
    {
        scale = 1 / newBoxSize;
        radius2 = newBoxSize * newBoxSize;
        boxSize = newBoxSize;

        // cut the edges via floor to avoid too thin boundaries, which
        // may interfere with a Torus topology where some cells, for
        // instance on the left border would try to access their left
        // neighbors (actually on the right side of the domain), but
        // would be unable to find them since the resulting grid thin
        // boundary container cells on the right would not contain all
        // required neighbors.
        for (int i = 0; i < DIM; ++i) {
            discreteDim.c[i] = std::max(1.0, std::floor(dimensions.c[i] * scale));
        }

        if (discreteDim.prod() > MAX_SIZE)
            throw std::logic_error("too many container cells are required");
    }

    bool searchList(
        const std::list<FloatCoord<DIM> >& list,
        const FloatCoord<DIM>& pos,
        std::set<int> *coords = 0) const
    {
        bool found = false;

        for (typename std::list<FloatCoord<DIM> >::const_iterator iter = list.begin();
             iter != list.end();
             ++iter) {
            if (distance2(pos, *iter) < radius2) {
                found = true;
                if (coords)
                    coords->insert(iter->id);
            }
        }

        return found;
    }

    /**
     * returns the square of the euclidean distance of a and b.
     */
    double distance2(const FloatCoord<DIM>& a, const FloatCoord<DIM>& b) const
    {
        double dist2 = 0;

        for (int i = 0; i < DIM; ++i) {
            double delta = std::abs(a.c[i] - b.c[i]);
            if (TOPOLOGY::wrapsAxis(i))
                delta = std::min(delta, dimensions.c[i] - delta);
            dist2 += delta * delta;
        }

        return dist2;
    }

    int manhattanDistance(const FloatCoord<DIM>& a, const FloatCoord<DIM>& b) const
    {
        Coord<DIM> coordA = posToCoord(a);
        Coord<DIM> coordB = posToCoord(b);
        Coord<DIM> delta = coordA - coordB;
        int maxDist = 0;
        
        for (int i = 0; i < DIM; ++i) {
            int dist = std::abs(delta.c[i]);
            if (TOPOLOGY::wrapsAxis(i))
                dist = std::min(dist, discreteDim.c[i] - dist);
            maxDist = std::max(dist, maxDist);
        }

        return maxDist;
    }

    bool intervalTooLarge(const Coord<DIM>& a, const Coord<DIM>& b)
    {
        int maxDist = 0;
        for (int i = 0; i < DIM; ++i) {
            maxDist = std::max(maxDist, std::abs(a.c[i] - b.c[i]));
        }
        return maxDist > 1;

    }
};

}

#endif
