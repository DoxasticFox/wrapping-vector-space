#include "Wvs.h"

Wvs::Wvs(unsigned maxPoints, unsigned dimensionality)
	: MAX_POINTS(maxPoints)
{
	this->dimensionality = dimensionality;
	
	this->numPoints = 0;
	
	this->points = new float*[MAX_POINTS];
	#pragma omp parallel for schedule(static,256) num_threads(4)
	for (unsigned i = 0; i < this->MAX_POINTS; i++) {
		this->points[i] = new float[dimensionality];
	}
}
float Wvs::computeDistance(unsigned i, unsigned j)
{
	return this->computeDistance(this->points[i], this->points[j]);
}
float Wvs::computeDistance(float* p, float* q)
{
	float distance = this->computeSquaredDistance(p, q);
	distance = sqrt(distance);
	return distance;
}
float Wvs::computeSquaredDistance(unsigned i, unsigned j)
{
	return this->computeSquaredDistance(this->points[i], this->points[j]);
}
float Wvs::computeSquaredDistance(float* p, float* q)
{
	float distance = 0.0;
	for (unsigned i = 0; i < this->dimensionality; i++) {
		float distanceI = this->msub(p[i], q[i]);
		distance += distanceI * distanceI;
	}
	return distance;
}
float Wvs::computeEnergy()
{
	float energy = 0.0;
	for (unsigned i = 0; i < this->numPoints; i++) {
		for (unsigned j = i + 1; j < this->numPoints; j++) {
			energy += this->computeEnergyInPoints(i, j);
		}
	}
	return energy;
}
float Wvs::computeEnergyDueToPoint(unsigned i)
{
	float energy = 0.0;
	#pragma omp parallel for reduction(+:energy) schedule(static,256) num_threads(4)
	for (unsigned i_ = 0; i_ < this->numPoints; i_++) {
		energy += this->computeEnergyInPoints(i, i_);
	}
	return energy;
}
float Wvs::computeEnergyInPoints(unsigned i, unsigned j)
{
	float distance = this->computeDistance(i, j);
	if (distance <= std::numeric_limits<float>::epsilon()) {
		return 0.0;
	}
	return - this->getWeight(i, j) / distance;
}
float Wvs::computeLength(unsigned i)
{
	return this->computeLength(this->points[i]);
}
float Wvs::computeLength(float* p)
{
	float length = 0.0;
	for (unsigned i = 0; i < this->dimensionality; i++) {
		length += p[i] * p[i];
	}
	length = sqrt(length);
	return length;
}
float Wvs::randUniform(float range)
{
	int rand = std::rand() - RAND_MAX / 2; // interval [-RAND_MAX/2, RAND_MAX/2]
	float rand_ = rand / (RAND_MAX / 2.0); // interval [-1, 1]
	return rand_ * range;
}
float Wvs::getDefaultWeight()
{
	return sqrt(2.0);
}
float* Wvs::getPoint(unsigned i)
{
	return this->points[i];
}
float* Wvs::getPointClone(unsigned i)
{
	float* point = new float[this->dimensionality];
	this->copyPoint(this->points[i], point);
	return point;
}
std::string Wvs::pointToString(float* p)
{
	if (this->dimensionality == 0) {
		return "[]";
	}

	typedef std::numeric_limits<float> fl;
	std::ostringstream ss;
	ss.precision(fl::digits10);

	ss << "[ ";
	for (unsigned i = 0; i < this->dimensionality - 1; i++) {
		ss << std::fixed << p[i];
		ss << ",";
	}
	ss << std::fixed << p[this->dimensionality - 1];
	ss << "]";
	return ss.str();
}
std::string Wvs::pointToString(unsigned i)
{
	return this->pointToString(this->points[i]);
}
std::string Wvs::pointsToString()
{
	std::ostringstream ss;
	for (unsigned i = 0; i < this->numPoints; i++) {
		ss << this->pointToString(i) << std::endl;
	}
	return ss.str();
}
std::string Wvs::pointsToPyDict()
{
	if (this->numPoints == 0) {
		return "";
	}

	std::ostringstream ss;
	ss << "{" << std::endl;
	for (unsigned i = 0; i < this->numPoints - 1; i++) {
		ss << "\t";
		ss << "\"" << std::to_string(i) << "\": ";
		ss << this->pointToString(i);
		ss << "," << std::endl;
	}
	ss << "\t";
	ss << "\"" << std::to_string(this->numPoints - 1) << "\": ";
	ss << this->pointToString(this->numPoints - 1);
	ss << std::endl;
	ss << "}";
	return ss.str();
}
unsigned Wvs::computeNearestPoint(unsigned i)
{
	return this->computeNearestPoint(this->points[i]);
}
unsigned Wvs::computeNearestPoint(float* p)
{
	float minDistance = std::numeric_limits<float>::max();
	unsigned minPoint;
	for (unsigned i = 0; i < this->numPoints; i++) {
		if (p == this->points[i]) {
			continue;
		}
		float thisDistance = this->computeDistance(p, this->points[i]);
		if (thisDistance <= minDistance) {
			minDistance = thisDistance;
			minPoint = i;
		}
	}
	return minPoint;
}
/*
unsigned* Wvs::computeNearestPoints(unsigned i, unsigned n)
{
	
}
*/
unsigned Wvs::computeFurthestPoint(unsigned i)
{
	return this->computeFurthestPoint(this->points[i]);
}
unsigned Wvs::computeFurthestPoint(float* p)
{
	float maxDistance = std::numeric_limits<float>::min();
	unsigned maxPoint;
	for (unsigned i = 0; i < this->numPoints; i++) {
		float thisDistance = this->computeDistance(p, this->points[i]);
		if (thisDistance > maxDistance) {
			maxDistance = thisDistance;
			maxPoint = i;
		}
	}
	return maxPoint;
}
unsigned Wvs::getNumPoints()
{
	return this->numPoints;
}
unsigned Wvs::getDimensionality()
{
	return this->dimensionality;
}
unsigned Wvs::getMaxPoints()
{
	return this->MAX_POINTS;
}
unsigned Wvs::addRandomPoint()
{
	float* point = new float[this->dimensionality];
	for (unsigned i = 0; i < this->dimensionality; i++) {
		float rand;
		rand = (float) std::rand();
		rand /= (float) RAND_MAX;	// interval [0, 1]
		rand *= 2.0f;				// interval [0, 2]

		point[i] = rand;
	}

	unsigned pointNumber = this->addPoint(point);
	delete point;
	return pointNumber;
}
unsigned Wvs::addPoint(float* p)
{
	this->numPoints++;
	this->updatePoint(p, this->numPoints);
	return this->numPoints;
}
void Wvs::computeNearestPoints(unsigned &i, unsigned &j)
{
	float minDistance = std::numeric_limits<float>::max();
	for (unsigned i_ = 0; i_ < this->numPoints; i_++) {
		unsigned j_ = this->computeNearestPoint(i_);
		float thisDistance = this->computeDistance(i_, j_);
		if (thisDistance < minDistance) {
			minDistance = thisDistance;
			i = i_;
			j = j_;
		}
	}
}
void Wvs::computeFurthestPoints(unsigned &i, unsigned &j)
{
	float maxDistance = std::numeric_limits<float>::min();
	for (unsigned i_ = 0; i_ < this->numPoints; i_++) {
		unsigned j_ = this->computeFurthestPoint(i_);
		float thisDistance = this->computeDistance(i_, j_);
		if (thisDistance > maxDistance) {
			maxDistance = thisDistance;
			i = i_;
			j = j_;
		}
	}
}
void Wvs::copyPoint(float* source, float* dest)
{
	std::memcpy(dest, source, this->dimensionality * sizeof(float));
}
void Wvs::perturbPoint(float* p, float d)
{
	for (unsigned i = 0; i < this->dimensionality; i++) {
		p[i] += ((rand() % 3) - 1) * d;
	}
}

void Wvs::solve()
{
	float totalEnergyRise = 0.0;

	for (int i = 0; i < 10; i++) {
		for (int j = 10; j >= 0; j--) {
			float energyRise = this->solveStep(2.0 * pow(j+3.5, -3.5));
			//std::cout << energyRise << std::endl;
			totalEnergyRise += energyRise;
		}
		//std::cout << std::endl;
	}

	std::cout << totalEnergyRise << " *" << std::endl;
}
float Wvs::solveStep(float inc)
{
	float totalEnergyRise = 0.0;

	for (unsigned i = 0; i < this->numPoints; i++) {
		float energyRise = 0.0;
		
		float* savedPoint = this->getPointClone(i);

		energyRise -= this->computeEnergyDueToPoint(i);
		this->perturbPoint(this->points[i], inc);
		energyRise += this->computeEnergyDueToPoint(i);

		if (energyRise >= 0) {
			// Keep change
			totalEnergyRise += energyRise;
		} else {
			// Undo change
			this->updatePoint(savedPoint, i);
		}
		delete savedPoint;
	}

	return totalEnergyRise;
}
void Wvs::updatePoint(float* p, unsigned i)
{
	this->copyPoint(p, this->points[i]);
}
float Wvs::madd(float* p, float* q, float* out)
{
	for (unsigned i = 0; i < this->dimensionality; i++) {
		out[i] = this->madd(p[i], q[i]);
	}
}
float Wvs::msub(float* p, float* q, float* out)
{
	for (unsigned i = 0; i < this->dimensionality; i++) {
		out[i] = this->msub(p[i], q[i]);
	}
}
float Wvs::madd(float a, float b)
{
	a += b;
	return a + 2 * (a < -1) - 2 * (a > 1);
}
float Wvs::msub(float a, float b)
{
	a += b;
	return 2 * (a < 0) + (a % 2);
}
