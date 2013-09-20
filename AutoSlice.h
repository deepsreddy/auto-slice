#pragma once
#include "lengthprofile.h"
#include "PreselectionItem.h"

using namespace System::Collections::Generic;

#define VERY_SMALL       (1e-12) 

//Length independent thresholding parameters
//#define MIN_NUMBEROFMOLECULES	(100.0)
//#define TROUGHTOPEAKRATIO_LIMIT	(0.7)

//#define ADD 1
//#define DISCARD 0
//Length dependent thresholding parameters
#define TROUGHTOPEAKRATIO_LOWERLENGTH	(0.75)
#define TROUGHTOPEAKRATIO_UPPERLENGTH (0.55)
#define NUMBEROFMOLECULES_LOWERLENGTH	(150)
#define NUMBEROFMOLECULES_UPPERLENGTH (100)
#define BOUNDARY_LENGTH	(100.0)
#define MINIMUM_LOWERPERCENTMOLECULES (3.0)
#define MINIMUM_UPPERPERCENTMOLECULES (1.0)

//#define PERCENTMOLECULES_LOWERLENGTH (0.10)
//#define PERCENTMOLECULES_UPPERLENGTH (0.05)
//#define PERCENTMOLECULES_LOWERSLICETHRESHOLD (0.10)
//#define PERCENTMOLECULES_UPPERSLICETHRESHOLD (0.05)
//#define PERCENTLOWERLENGTHMOLECULES_JUNKTHRESHOLD (0.20)
//#define PERCENTHIGHERLENGTHMOLECULES_JUNKTHRESHOLD (0.150)
//#define BOUNDARY_LENGTH	(80.0)

namespace AutoSlicing {

	public ref class CAutoSlice : public CLengthProfile
	{
		public:
			CAutoSlice(void);
			CAutoSlice(String ^sDirectory);
			bool PerformAutoSlicing(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, String ^sErr);

			//struct for properties of final slices
			ref struct finalSliceProperties
			{
				bool _real;
				int	_sliceWidth;
				double _leftBoundary;
				double _rightBoundary;
				double _totalMolecules;
				double _minMolecules;
				double _percentMolecules;
			};
			List<finalSliceProperties ^>	^finalSlices;
			List<PreselectionItem ^>		^preselectionItems;
		private:
			ref struct profileProperties
			{
				int		_iPosition;
				double	_dCurrentSignal;
				double	_dPreviousSignal;
				double	_dNextSignal;
			};

			ref struct sliceProperties
			{
				int _iStart;
				int _iEnd;
				int _iPeak;
				double	_peakSignal; 
				double _totalMolecules;
				double _minMolecules;
				double _percentMolecules;
				double _trough2peakratio1;
				double _trough2peakratio2;
				double _trough1Position;
				double _trough2Position;
				double _peakPosition;
			};

			String ^_sDirectory;
			List<profileProperties ^> ^peakProfiles;
			List<profileProperties ^> ^troughProfiles;
			List<sliceProperties ^> ^sliceProfiles;
			double _UpperSliceMolecules;
			double _LowerSliceMolecules;
			int _iLowerSliceCount;
			int _iUpperSliceCount;
			double _AverageSliceMolecules;

			bool EstimatePeakTroughs(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, String ^sErr);
			bool EstimateTrough2PeakRatios(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, String ^sErr);
			bool EvaluatePercentMolecules(const double, const double, array<LengthHistogram ^> ^lengthHisto1);
			bool AutoBoundaryDetection(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, String ^sErr);
			
//			bool evaluateLengthProfile(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, String ^sErr);
//			bool cleanProfiles(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, String ^sErr);
//			bool autoBoundaryDetection(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, String ^sErr);
//			bool check4Slices(int &iCurrentSlicePosition, int &iEndSlicePosition);
//			bool check4CurrentJunkSlice(int &iCurrentSlicePosition, int iSlicePosition, double dThreshold);
//			bool check4NextJunkSlice(int &iCurrentSlicePosition, int iNumberMoleculesComparison);
//			bool advanceSlices(int &iCurrentSlicePosition);
//			//bool check4UpperLengthSlice(int &iCurrentSlicePosition, int &iEndSlicePosition);
	};
}

