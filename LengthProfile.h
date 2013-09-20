#pragma once
#include <direct.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include "LengthHistogram.h"
#include "TraceMap.h"
using namespace std;
using namespace System;
using namespace System::Runtime::InteropServices;


#define LEFTFRAMESIZE (2)
#define RIGHTFRAMESIZE (2)
#define ORDEROFDERIVATIVE (0)
#define POLYNOMIALORDER (2)

namespace AutoSlicing {

	//struct LengthHistogram holds the plain summatiom, weightedaverage and smoothedweightedaverage values
	public ref class CLengthProfile
	{
		public:
			CLengthProfile(void);
			bool GenerateLengthProfile(array<LengthHistogram ^> ^lengthHisto, array<CRawData ^> ^vRawData, const int xBins, int yBins, String ^sDirectory);
			bool Generate2ndorderSmoothedWeightedAverages (array<LengthHistogram ^> ^lengthHisto, String ^sDirectory, int xBins);		
			bool SmoothingInterface(array<LengthHistogram ^> ^lengthHisto, const int xBins);
			bool WriteRScript(String ^sDirectory, String ^sErr);

		public:
			array<LengthHistogram ^> ^lengthHisto;  //array holding the length histogram
	};
}//end of namespace AutoSlicing

