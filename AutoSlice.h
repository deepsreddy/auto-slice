#pragma once
#include "lengthprofile.h"
#include "PreselectionItem.h"
#include "AutoSliceOptions.h"
using namespace System::Collections::Generic;

#define VERY_SMALL       (1e-12)
#define EQUALITY_TOLERANCE (1e-6)

namespace AutoSlicing {

	public ref class CAutoSlice : public CLengthProfile
	{
		public:
			CAutoSlice(void);
			CAutoSlice(String ^sDirectory, const CAutoSliceOptions ^cAutoSliceOptions);
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
			int	_iFinalSliceCount;
			double _AverageSliceMolecules;
			CAutoSliceOptions sXMLOptions;

			bool EstimatePeakTroughs(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, String ^sErr);
			bool EstimateTrough2PeakRatios(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, String ^sErr);
			bool EvaluatePercentMolecules(const double, const double, array<LengthHistogram ^> ^lengthHisto1);
			bool AutoBoundaryDetection(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, String ^sErr);
			int FindNextPosition(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, int iCurrentPosition);
			int FindPreviousPosition(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, int iCurrentPosition);
			bool AppendSuperSlice(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1);
			int GetIndex(int iCount, int iNthSlice);
			double GetSumMoleculeDensity(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, int iStartIndex, int iEndIndex);
			void OutputBoundaries();
			void CopySlices();
	};
}

