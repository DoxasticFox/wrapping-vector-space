#ifndef WVS_H
#define WVS_H

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream> /* Only so that `solve` can echo energy rises */
#include <limits>
#include <sstream>
#include <string>

class Wvs {
/* IMPORTANT: In the interests of speed, this code has been written missing many
 * basic safety features like rangechecks and exceptions. Be careful using it!
 */
private:
	const unsigned MAX_POINTS;
	unsigned dimensionality;
	unsigned numPoints;
	float**  points;
public:
	            Wvs(unsigned maxPoints, unsigned dimensionality=3);
	float       computeDistance(float* p, float* q);
	float       computeDistance(unsigned i, unsigned j);
	float       computeSquaredDistance(float* p, float* q);
	float       computeSquaredDistance(unsigned i, unsigned j);
	float       computeEnergy();
	float       computeEnergyDueToPoint(unsigned  i);
	float       computeEnergyInPoints(unsigned  i, unsigned j);
	float       computeLength(float* p);
	float       computeLength(unsigned i);
	float*      getPoint(unsigned i);
	float*      getPointClone(unsigned i);
	std::string pointToString(float* p);
	std::string pointToString(unsigned i);
	std::string pointsToString();
	std::string pointsToPyDict();
	unsigned    computeNearestPoint(float* p);
	unsigned    computeNearestPoint(unsigned    i);
	unsigned*   computeNearestPoints(unsigned i, unsigned n);
	unsigned    computeFurthestPoint(float* p);
	unsigned    computeFurthestPoint(unsigned    i);
	unsigned    getNumPoints();
	unsigned    getDimensionality();
	unsigned    getMaxPoints();
	unsigned    addRandomPoint();
	unsigned    addPoint(float* p);
	void        computeNearestPoints(unsigned &i, unsigned &j);
	void        computeFurthestPoints(unsigned &i, unsigned &j);
	void        copyPoint(float* source, float* dest);
	void        perturbPoint(float* p, float d);
	void        updatePoint(float* p, unsigned i);
	void        solve();
	float       solveStep(float inc);
	void        madd(float* a, float* b, float* out);
	void        msub(float* a, float* b, float* otu);
	float       madd(float a, float b);
	float       msub(float a, float b);
};

#endif
