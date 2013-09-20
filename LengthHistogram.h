#pragma once
#include <direct.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
using namespace std;
using namespace System;
using namespace System::Runtime::InteropServices;

namespace AutoSlicing 
{

	//struct LengthHistogram holds the plain summatiom, weightedaverage and smoothedweightedaverage values
	public ref struct LengthHistogram
	{
		double _XCoordinate;
		double _densitySummation;		
		double _weightedDensityAverage;	
		double _smoothedWeightedAverage;
		char   _PKFlag;
	};

}//end of namespace AutoSlicing

