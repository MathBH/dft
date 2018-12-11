#include "stdafx.h"
#include "MultiDimensionalCounter.h"

MultiDimensionalCounter::MultiDimensionalCounter()
{
}

MultiDimensionalCounter::MultiDimensionalCounter(std::vector<int> dimensions)
{
	numDimensions = dimensions.size();
	indicies = std::vector<int>(numDimensions,0);
	idxMax = dimensions;
	capped = false;
}


MultiDimensionalCounter::~MultiDimensionalCounter()
{
}

void MultiDimensionalCounter::increment()
{
	int d = 0;
	int capBuff;
	do {
		capBuff = idxMax[d];
		indicies[d] ++;
		indicies[d] %= capBuff;
	} while (indicies[d] == 0 && ++d < numDimensions);
	if (d == numDimensions && indicies[numDimensions - 1] == 0)
		capped = true;
}

std::vector<int> MultiDimensionalCounter::nextPos()
{
	std::vector<int> output = std::vector<int>(indicies);
	increment();
	return output;
}

bool MultiDimensionalCounter::hasNext()
{
	return !capped;
}
