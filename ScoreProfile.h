#pragma once
#include <direct.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include "TraceMap.h"
#include "AutoSlice.h"
#include "AutoSliceOptions.h"

using namespace std;
using namespace System;
using namespace System::Runtime::InteropServices;

namespace AutoSlicing {

	public ref class CScoreProfile
	{
		public:
			CScoreProfile(void);
			CScoreProfile(const CAutoSliceOptions ^cAutoSliceOptions);
			bool GenerateScoreProfile(array<CRawData ^> ^vRawData, int xBins, int yBins, String ^sDirectory, double xBoundary, double yBoundary, int sliceNumber);

			double yBorders, _YMin, _YMax;

		private:
			ref struct ScoreHistogram
			{
				double _YCoordinate;
				double _densitySummation;
				//double _densityAverage;	
			};

			array<ScoreHistogram ^> ^scoreHisto;
			double _totalMolecules;
			double FindYThreshold(int yBins, double dThresholdPercent);
			CAutoSliceOptions sXMLOptions;
	};
}//end of namespace AutoSlicing
