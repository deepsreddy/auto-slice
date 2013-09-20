#pragma once
#include <direct.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include "LengthHistogram.h"
#include "PreselectionItem.h"

using namespace std;
using namespace System;
using namespace System::Runtime::InteropServices;
using namespace System::Collections::Generic;

namespace AutoSlicing {

	//CRawData for holding the raw length ('x'), score ('y'), molecule density ('z') values
	public ref class CRawData
	{
		public:
			double	_moleculelength;
			double	_distancescore;
			double	_moleculedensity;
	};

	public ref class CTraceMap
	{
		public:
			CTraceMap(void);
			CTraceMap(String ^sDirectory);

			bool AssembleRawData(String ^sError, array<double> ^xInterp, array<double> ^yInterp, array<double, 2> ^zInterp);
			
			List <PreselectionItem^>^  EvaluateProfiles(String ^sErr);

            void WriteTrace2TraceMap(const int xBins, const int yBins);

		private:
			int xBins;
			int yBins;
			String ^_sDirectory;
			array<CRawData ^> ^rawData;
			array<LengthHistogram ^> ^lengthHisto;  //array holding the length histogram

			//FSSTODO remove obsolete
			//bool ReadFirstThreeLines(ifstream &rawTraceFile);
			//bool WriteRScript(String ^sDirectory, String ^sErr);
	};
} //end of namespace AutoSlicing

