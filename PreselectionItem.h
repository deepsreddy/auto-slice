#pragma once
#include <direct.h>
#include <vector>
#include <string>

using namespace std;
using namespace System;
using namespace System::Runtime::InteropServices;

namespace AutoSlicing 
{
    public ref class PreselectionItem
    {
		public:

		PreselectionItem(void);
		PreselectionItem(double lowerLenLimit, double upperLenLimit, double cutoffDistanceValue, double cutoffDistanceMin, bool _real);

		public:
        double LowerLenLimit;
        double UpperLenLimit;
        double CutoffDistanceValue;
        double CutoffDistanceMin;
		bool _Real;

	};

} //end of namespace AutoSlicing

