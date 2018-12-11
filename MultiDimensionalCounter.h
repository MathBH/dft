#pragma once
#include <vector>
/*
	MultiDimensionalCounter

	Object that serves the purpose of automating multi-dimensional iteration.
	Specify number of dimensions and size of each and calls to nextPos() will output a
	vector of indicies for the multidimensional space in linear order.
*/
class MultiDimensionalCounter
{
private:
	int numDimensions;
	bool capped;
	std::vector<int> indicies;
	std::vector<int> idxMax;
	void increment();
public:
	MultiDimensionalCounter();
	MultiDimensionalCounter(std::vector<int> dimensions);
	~MultiDimensionalCounter();

	std::vector<int> nextPos();
	bool hasNext();
};

