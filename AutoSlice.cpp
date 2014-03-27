#include "stdafx.h"
#include "AutoSlice.h"
#include "BasicFunctions.h"
using namespace System;
using namespace System::Collections::Generic;

namespace AutoSlicing {

	//Default constructor
	CAutoSlice::CAutoSlice(void)
	{

	}

	//Constructor for initializing output directory
	CAutoSlice::CAutoSlice(String ^sDirectory, const CAutoSliceOptions ^cAutoSliceOptions)
	{
		_sDirectory = sDirectory;
		_UpperSliceMolecules = 0.0;
		_LowerSliceMolecules = 0.0;
		_iLowerSliceCount = 0;
		_iUpperSliceCount = 0;

		sXMLOptions.operator=(cAutoSliceOptions);
	}

	//Auto Slicing - external call for performing autoslicing, 
	bool CAutoSlice::PerformAutoSlicing(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, String ^sErr)
	{
		if(!EstimatePeakTroughs(iLengthHistoSize, lengthHisto1, sErr))
		{
			sErr->Format("Something wrong with evaluating profile");
			return false;
		}

		if(!EstimateTrough2PeakRatios(iLengthHistoSize, lengthHisto1, sErr))
		{
			sErr->Format("Something wrong with evaluating profile");
			return false;
		}

		//Check for super slice
		if(sXMLOptions._bSuperSliceOption > 0 && _iFinalSliceCount > 2)
		{
			if(!AppendSuperSlice(iLengthHistoSize, lengthHisto1))
				Console::WriteLine(L"Cannot append the super slice");
		}
		
		//Output boundaries to Boundaries.txt
		OutputBoundaries();

		//Copy slice boundary information to FSS Preselection list
		CopySlices();

		return true;
	}

	//evaluate Length Profile for peaks and troughs
	bool CAutoSlice::EstimatePeakTroughs(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, String ^sErr)
	{
		peakProfiles = gcnew List< profileProperties ^ >();
		troughProfiles = gcnew List< profileProperties ^ >();

		//Make sure to sweep the length histo for cleaning up
		
		double diff1 = 0.0, diff2 = 0.0;
		double dCurrentSignal;
		double dNextSignal;
		double dPreviousSignal;
		int iCurrentPosition = 0, iNextPosition, iPreviousPosition;

		profileProperties^ initialProfile = gcnew profileProperties();

		//  FSSTODO JEP This next loop was added back after we removed something similar to it during the initial port
		while((dCurrentSignal = lengthHisto1[iCurrentPosition]->_smoothedWeightedAverage) <= VERY_SMALL)
			iCurrentPosition++;

		//Go to the previous '0' molecules position
		if(iCurrentPosition > 0)
			iCurrentPosition--;
		initialProfile->_dCurrentSignal = lengthHisto1[iCurrentPosition]->_smoothedWeightedAverage;
		initialProfile->_dNextSignal = lengthHisto1[iCurrentPosition+1]->_smoothedWeightedAverage;

		if(iCurrentPosition > 0)
		{
			initialProfile->_dPreviousSignal = lengthHisto1[iCurrentPosition-1]->_smoothedWeightedAverage;
			initialProfile->_iPosition = iCurrentPosition - 1;
		}
		else
		{
			initialProfile->_dPreviousSignal = -1;
			initialProfile->_iPosition = iCurrentPosition;
		}

		lengthHisto1[iCurrentPosition]->_PKFlag = 'T';
		troughProfiles->Add(initialProfile);
		iCurrentPosition++;

		//Estimate 'peaks' or 'troughs' based on the local derivative
		while (iCurrentPosition < (iLengthHistoSize-1))
		{
			dCurrentSignal = lengthHisto1[iCurrentPosition]->_smoothedWeightedAverage;

			//dNextSignal = lengthHisto1[iCurrentPosition+1]->_smoothedWeightedAverage;
			iNextPosition = FindNextPosition(iLengthHistoSize, lengthHisto1, iCurrentPosition);
			dNextSignal = lengthHisto1[iNextPosition]->_smoothedWeightedAverage;

			//dPreviousSignal = lengthHisto1[iCurrentPosition-1]->_smoothedWeightedAverage;
			iPreviousPosition = FindPreviousPosition(iLengthHistoSize, lengthHisto1, iCurrentPosition);
			dPreviousSignal = lengthHisto1[iPreviousPosition]->_smoothedWeightedAverage;

			diff1 = dPreviousSignal - dCurrentSignal;
			diff2 = dNextSignal - dCurrentSignal;

			if((diff1 > 0 && diff2 > 0) || ((diff1 < 0 && diff2 < 0)))
			{
				profileProperties^ tempProfile = gcnew profileProperties();

				tempProfile->_dCurrentSignal = dCurrentSignal;
				tempProfile->_dNextSignal = dNextSignal;
				tempProfile->_dPreviousSignal = dPreviousSignal;
				tempProfile->_iPosition = iCurrentPosition;

				if(diff1 > 0 && diff2 > 0)
				{
					lengthHisto1[iCurrentPosition]->_PKFlag = 'T';
					troughProfiles->Add(tempProfile);
				}
				else if(diff1 < 0 && diff2 < 0)
				{
					lengthHisto1[iCurrentPosition]->_PKFlag = 'P';
					peakProfiles->Add(tempProfile);
				}
			}

			//Advance current position
			iCurrentPosition = FindNextPosition(iLengthHistoSize, lengthHisto1, iCurrentPosition);
		}

		//Go back to the last position
		if(iCurrentPosition > (iLengthHistoSize-1))
			iCurrentPosition = (iLengthHistoSize-1);

		//Add last position as trough
		profileProperties^ tempProfile = gcnew profileProperties();
		dCurrentSignal = lengthHisto1[iCurrentPosition]->_smoothedWeightedAverage;
		dPreviousSignal = lengthHisto1[iCurrentPosition-1]->_smoothedWeightedAverage;
		tempProfile->_dCurrentSignal = dCurrentSignal;
		tempProfile->_dNextSignal = -1;
		tempProfile->_dPreviousSignal = dPreviousSignal;
		tempProfile->_iPosition = iCurrentPosition;
		troughProfiles->Add(tempProfile);
		lengthHisto1[iCurrentPosition]->_PKFlag = 'T';


		//Console::WriteLine(L"Number of peak points: {0}", peakProfiles->Count);
		//Console::WriteLine(L"Number of trough points: {0}", troughProfiles->Count);

		//Outputting the peak/trough positions in a file
		//FSS: No need to output this file
		String ^fileName = gcnew String("PeakTroughProfiles.txt");
		fileName = _sDirectory + fileName;

		//Convert System::String to char array to read in ofstream
		ofstream profileFile(AutoSlicing::CBasicFunctions::SystemStringToCharArray(fileName));

		profileFile << "moleculelength\tweightedDensity\tsmoothedDensity\tPKflag\n";

		iCurrentPosition = 0;
		while(iCurrentPosition < iLengthHistoSize)
		{
			profileFile << lengthHisto1[iCurrentPosition]->_XCoordinate << "\t" << lengthHisto1[iCurrentPosition]->_weightedDensityAverage << "\t" << lengthHisto1[iCurrentPosition]->_smoothedWeightedAverage << "\t";
			if(lengthHisto1[iCurrentPosition]->_PKFlag == 'P' || lengthHisto1[iCurrentPosition]->_PKFlag == 'T') 
				profileFile << lengthHisto1[iCurrentPosition]->_PKFlag << "\n";
			else
				profileFile << "NA" << "\n";
			iCurrentPosition++;
		}

		profileFile.close();
		profileFile.clear();
		
		return true;
	}

	//Estimate next position
	int CAutoSlice::FindNextPosition(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, int iCurrentPosition)
	{
		int iNextPosition = iCurrentPosition+1;
		double dDiff;

		if(iNextPosition <= (iLengthHistoSize-1))
			dDiff = abs(lengthHisto1[iNextPosition]->_smoothedWeightedAverage - lengthHisto1[iCurrentPosition]->_smoothedWeightedAverage);
		else
			return iCurrentPosition;

		while(dDiff < EQUALITY_TOLERANCE && iNextPosition < (iLengthHistoSize-1))
		{
			iNextPosition++;
			dDiff = abs(lengthHisto1[iNextPosition]->_smoothedWeightedAverage - lengthHisto1[iCurrentPosition]->_smoothedWeightedAverage);
		}

		return iNextPosition;
	}

	//Estimate previous position
	int CAutoSlice::FindPreviousPosition(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, int iCurrentPosition)
	{
		int iPreviousPosition = iCurrentPosition-1;
		double dDiff;

		if(iCurrentPosition > 0)
			dDiff = abs(lengthHisto1[iPreviousPosition]->_smoothedWeightedAverage - lengthHisto1[iCurrentPosition]->_smoothedWeightedAverage);

		while(dDiff < EQUALITY_TOLERANCE && iPreviousPosition > 0)
		{
			iPreviousPosition--;
			dDiff = abs(lengthHisto1[iPreviousPosition]->_smoothedWeightedAverage - lengthHisto1[iCurrentPosition]->_smoothedWeightedAverage);
		}

		return iPreviousPosition;
	}

	//Clean profiles for adding slices
	//The name of this function might be misleading. This function basically creates preliminary slices by taking trough-peak-trough positions
	//and calculates the trough1/peak & trough2/peak ratios. Each of these preliminary slices is populated in the list <sliceProfiles>
	bool CAutoSlice::EstimateTrough2PeakRatios(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, String ^sErr)
	{
		//sliceProfiles holds the pseudo-slices with ratios and other properties
		sliceProfiles = gcnew List<sliceProperties ^>();
		int iCount, iDiff;

		//Making sure the index is not going out of bound for samples with troughs <= peaks
		if(peakProfiles->Count < troughProfiles->Count)
			iCount = peakProfiles->Count;
		else
		{
			iDiff = (peakProfiles->Count - troughProfiles->Count);
			if(iDiff == 0)
				iCount = (peakProfiles->Count-1);
			else
				iCount = (peakProfiles->Count-iDiff-1);
		}

		int i = 0, j, iStart, iEnd, iPeak;
		double tMolecules = 0.0, ratio1, ratio2, totalMolecules = 0.0, upperMolecules = 0.0, lowerMolecules = 0.0;
		while(i < (iCount))
		{
			iStart = troughProfiles[i]->_iPosition;
			iEnd = troughProfiles[i+1]->_iPosition;
			iPeak = peakProfiles[i]->_iPosition;
			tMolecules = 0.0;
			j = iStart;
			while(j < iEnd)
			{
				tMolecules += lengthHisto1[j]->_densitySummation;
				j++;
			}

			sliceProperties^ tempSlice = gcnew sliceProperties();
			tempSlice->_iStart = iStart;
			tempSlice->_iPeak = iPeak;
			tempSlice->_iEnd = iEnd;
			tempSlice->_totalMolecules = tMolecules;
			tempSlice->_peakSignal = peakProfiles[i]->_dCurrentSignal;

			ratio1 = (troughProfiles[i]->_dCurrentSignal/peakProfiles[i]->_dCurrentSignal);
			ratio2 = (troughProfiles[i+1]->_dCurrentSignal/peakProfiles[i]->_dCurrentSignal);

			tempSlice->_trough2peakratio1 = ratio1;
			tempSlice->_trough2peakratio2 = ratio2;

			tempSlice->_trough1Position = lengthHisto1[iStart]->_XCoordinate;
			tempSlice->_trough2Position = lengthHisto1[iEnd]->_XCoordinate;
			tempSlice->_peakPosition = lengthHisto1[i]->_XCoordinate;

			if(tempSlice->_trough2Position >= sXMLOptions._dBoundaryLength)
				upperMolecules += tMolecules;
			else
				lowerMolecules += tMolecules;

			sliceProfiles->Add(tempSlice);
			
			i++;
		}

		_LowerSliceMolecules = lowerMolecules;
		_UpperSliceMolecules = upperMolecules;

		//Evaluate percent molecules
		EvaluatePercentMolecules(lowerMolecules, upperMolecules, lengthHisto1);
		AutoBoundaryDetection(iLengthHistoSize, lengthHisto1, sErr);

		return true;
	}

	//Evaluating percent of molecules within the preliminary slices
	bool CAutoSlice::EvaluatePercentMolecules(const double lowerMolecules, const double upperMolecules, array<LengthHistogram ^> ^lengthHisto1)
	{
		int iCount = sliceProfiles->Count;
		int i = 0, iStart, iEnd;
		while(i < iCount)
		{
			if(sliceProfiles[i]->_trough2Position >= sXMLOptions._dBoundaryLength)
				sliceProfiles[i]->_percentMolecules = ((sliceProfiles[i]->_totalMolecules/upperMolecules)*100.0);
			else
				sliceProfiles[i]->_percentMolecules = ((sliceProfiles[i]->_totalMolecules/lowerMolecules)*100.0);
			i++;
		}
		
		String ^fileName = gcnew String("peaktroughRatios.txt");
		fileName = _sDirectory + fileName;

		//Convert System::String to char array to read in ofstream
		ofstream profileFile((AutoSlicing::CBasicFunctions::SystemStringToCharArray(fileName)), std::iostream::out);
		profileFile << "moleculelength_start\tmoleculelength_end\ttotalMolecules\tPercentMolecules\ttroughDensity1\ttroughDensity2\tpeakDensity\ttrough1peakRatio\ttrough2peakRatio\n";

		iCount = sliceProfiles->Count;
		i = 0;
		while(i < iCount)
		{
			iStart = sliceProfiles[i]->_iStart;
			iEnd = sliceProfiles[i]->_iEnd;

			profileFile << lengthHisto1[iStart]->_XCoordinate << "\t" << lengthHisto1[iEnd]->_XCoordinate << "\t" << sliceProfiles[i]->_totalMolecules << "\t";
			profileFile << sliceProfiles[i]->_percentMolecules << "\t";
			profileFile << troughProfiles[i]->_dCurrentSignal << "\t" << troughProfiles[i+1]->_dCurrentSignal << "\t" << peakProfiles[i]->_dCurrentSignal << "\t";
			profileFile << (((double) troughProfiles[i]->_dCurrentSignal)/((double) peakProfiles[i]->_dCurrentSignal)) << "\t";
			profileFile << (((double) troughProfiles[i+1]->_dCurrentSignal)/((double) peakProfiles[i]->_dCurrentSignal)) << "\n";

			i++;
		}

		profileFile.close();
		profileFile.clear();

		return true;
	}

	//Detecting boundaries automatically from the preliminary slices calculated in the previous function and stored in sliceProfiles.
	//The List 'finalSlices' holds the final boundaries.
	bool CAutoSlice::AutoBoundaryDetection(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, String ^sErr)
	{
		int iNumberofSlices = sliceProfiles->Count;
		int iCurrentSlice = 0, iSliceCount = 0, iCurrIndex;
		double leftBoundary, minMolecules = 100000000, totalMolecules;
		double ratioThreshold, dMinPercentMolecules, dNumberMoleculesThreshold;
		finalSliceProperties ^tempSlice = gcnew  finalSliceProperties();
		finalSlices = gcnew List<finalSliceProperties ^>();

		while(iCurrentSlice < iNumberofSlices)
		{
			finalSliceProperties ^tempSlice = gcnew finalSliceProperties();
			leftBoundary = sliceProfiles[iCurrentSlice]->_trough1Position;

			if(sliceProfiles[iCurrentSlice]->_trough2Position > sXMLOptions._dBoundaryLength)
			{
				dNumberMoleculesThreshold = sXMLOptions._dMinimumNumberofMoleculesUpperLength;
				dMinPercentMolecules = sXMLOptions._dMinimumUpperPercentMolecules;
				ratioThreshold = sXMLOptions._dTroughtoPeakRatioUpperLengthLimit;
			}
			else
			{
				dNumberMoleculesThreshold = sXMLOptions._dMinimumNumberofMoleculesLowerLength;
				dMinPercentMolecules = sXMLOptions._dMinimumLowerPercentMolecules;
				ratioThreshold = sXMLOptions._dTroughtoPeakRatioLowerLengthLimit;
			}

			//Discarding junk slices
			if(sliceProfiles[iCurrentSlice]->_trough2peakratio1 > ratioThreshold && sliceProfiles[iCurrentSlice]->_trough2peakratio2 > ratioThreshold)
			{
				if(sliceProfiles[iCurrentSlice]->_percentMolecules < dMinPercentMolecules || sliceProfiles[iCurrentSlice]->_totalMolecules < dNumberMoleculesThreshold)
				{
					iCurrentSlice++;
					continue;
				}
			}
			
			iCurrIndex = iCurrentSlice;
			minMolecules = 100000000;
			totalMolecules = 0;

			
			//Check for the previous slice length and check for slice based on LOWER/UPPER length
			while((iCurrentSlice+1) < iNumberofSlices && (sliceProfiles[iCurrentSlice]->_trough2peakratio2 > ratioThreshold))
			{
				if(sliceProfiles[iCurrentSlice]->_trough2Position > sXMLOptions._dBoundaryLength)
					ratioThreshold = sXMLOptions._dTroughtoPeakRatioUpperLengthLimit;
				else
					ratioThreshold = sXMLOptions._dTroughtoPeakRatioLowerLengthLimit;

				
				if(sliceProfiles[iCurrentSlice+1]->_trough2peakratio1 <= ratioThreshold && sliceProfiles[iCurrentSlice]->_trough2peakratio2 <= ratioThreshold)
					break;

				iCurrentSlice++;
			}
			
			if(iCurrentSlice == iNumberofSlices)
				iCurrentSlice--;

			tempSlice->_leftBoundary = leftBoundary;
			tempSlice->_rightBoundary = sliceProfiles[iCurrentSlice]->_trough2Position;
			
			//Get total number of molecules in the slice
			while(iCurrIndex <= iCurrentSlice)
			{
				if(sliceProfiles[iCurrIndex]->_totalMolecules < minMolecules)
					minMolecules = sliceProfiles[iCurrIndex]->_totalMolecules;
				totalMolecules += sliceProfiles[iCurrIndex]->_totalMolecules;
				iCurrIndex++;
			}
			
			tempSlice->_minMolecules = minMolecules;
			tempSlice->_totalMolecules = totalMolecules;
			
			if(tempSlice->_rightBoundary > sXMLOptions._dBoundaryLength)
			{
				dNumberMoleculesThreshold = sXMLOptions._dMinimumNumberofMoleculesUpperLength;
				dMinPercentMolecules = sXMLOptions._dMinimumUpperPercentMolecules;
				tempSlice->_percentMolecules = (totalMolecules/_UpperSliceMolecules)*100;
			}
			else
			{
				dNumberMoleculesThreshold = sXMLOptions._dMinimumNumberofMoleculesLowerLength;
				dMinPercentMolecules = sXMLOptions._dMinimumLowerPercentMolecules;
				tempSlice->_percentMolecules = (totalMolecules/_LowerSliceMolecules)*100;
			}

			if(tempSlice->_totalMolecules <= dNumberMoleculesThreshold || tempSlice->_percentMolecules <= dMinPercentMolecules)
				tempSlice->_real = false;
			else
			{
				tempSlice->_real = true;
				iSliceCount++;
			}

			finalSlices->Add(tempSlice);
			iCurrentSlice++;
		}

		//Merge thin slices with slice width < 5.0 to the previous slice
		int iSlices = (finalSlices->Count);
		int index = (finalSlices->Count-1);
		double sliceWidth;
		while(index >= 0)
		{
			if(finalSlices[index]->_real)
			{
				sliceWidth = (finalSlices[index]->_rightBoundary - finalSlices[index]->_leftBoundary);
				if(sliceWidth < 5.0 && index > 0)
				{
					if((finalSlices[index]->_leftBoundary-finalSlices[index-1]->_rightBoundary) <= VERY_SMALL)
					{
						finalSlices[index-1]->_rightBoundary = finalSlices[index]->_rightBoundary;
						finalSlices[index]->_real = false;
						index--;
						iSlices--;
					}
				}
			}
			else
				iSlices--;
			index--;
		}

		_iFinalSliceCount = iSlices;

		return true;
	}

	//Check for if a superslice need to be appended
	bool CAutoSlice::AppendSuperSlice(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1)
	{
		int iSlices = (finalSlices->Count);
		int iLastSliceIndex = GetIndex(iSlices, _iFinalSliceCount);
		
		if(iLastSliceIndex < 0)
		{
			Console::WriteLine("Last slice index < 0");
			return false;
		}

		int iPenultimateSliceIndex = GetIndex(iSlices, (_iFinalSliceCount-1));

		if(iPenultimateSliceIndex < 0)
		{
			Console::WriteLine("Penultimate slice index < 0");
			return false;
		}

		int iStartIndex = (int) (finalSlices[iPenultimateSliceIndex]->_rightBoundary - lengthHisto1[0]->_XCoordinate);
		int iEndIndex = iLengthHistoSize;
		double dSuperSliceMolecules1 = GetSumMoleculeDensity(iLengthHistoSize, lengthHisto1, iStartIndex, iEndIndex);

		iStartIndex = (int) (finalSlices[iLastSliceIndex]->_rightBoundary - lengthHisto1[0]->_XCoordinate);
		double dSuperSliceMolecules2 = GetSumMoleculeDensity(iLengthHistoSize, lengthHisto1, iStartIndex, iEndIndex);

		finalSliceProperties ^tempSlice = gcnew  finalSliceProperties();
		if(dSuperSliceMolecules1 >= sXMLOptions._dSuperSliceMoleculesThreshold || dSuperSliceMolecules2 >= sXMLOptions._dSuperSliceMoleculesThreshold)
		{
			//finalSliceProperties ^tempSlice = gcnew  finalSliceProperties();
			tempSlice->_leftBoundary = finalSlices[iPenultimateSliceIndex]->_rightBoundary;
			tempSlice->_rightBoundary = lengthHisto1[iLengthHistoSize-1]->_XCoordinate;
			tempSlice->_totalMolecules = dSuperSliceMolecules1;
			tempSlice->_real = true;
			finalSlices->Add(tempSlice);
			_iFinalSliceCount++;
		}
		else
		{
			Console::WriteLine("#Molecules < Threshold: {0}, {1}, {2}", dSuperSliceMolecules1, tempSlice->_leftBoundary, iLengthHistoSize);
			return false;
		}

		return true;
	}

	//Get index of the nth slice
	int CAutoSlice::GetIndex(int iCount, int iNthSlice)
	{
		int iIndex = 0, iIndexNthSlice = 0;
		while(iIndex < iCount)
		{
			if(finalSlices[iIndex]->_real)
			{
				iIndexNthSlice++;

				if(iIndexNthSlice == iNthSlice)
					return iIndex;
			}

			iIndex++;
		}

		return -1;
	}

	//Get molecule density
	double CAutoSlice::GetSumMoleculeDensity(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, int iStartIndex, int iEndIndex)
	{
		double dSumDensity = 0.0;

		while(iStartIndex < iEndIndex && iStartIndex < iLengthHistoSize)
		{
			dSumDensity += lengthHisto1[iStartIndex]->_densitySummation;
			iStartIndex++;
		}

		return dSumDensity;
	}

	//Output boundaries file
	void CAutoSlice::OutputBoundaries()
	{
		//Outputting the boundaries to the Boundaries.txt
		//FSS: Please note that the number of real slices is not "finalSlices->Count" but "iSlices" as outputted in this file. 
		//FSS: These boundaries are used to output in the 'R' script file
		int index = 0;
		String ^fileName = gcnew String("Boundaries.txt");
		fileName = _sDirectory + fileName;

		//Convert System::String to char array to read in ofstream
		ofstream boundariesFile((AutoSlicing::CBasicFunctions::SystemStringToCharArray(fileName)), std::iostream::out);
		boundariesFile << "left_boundary\tright_boundary\n";
		boundariesFile << _iFinalSliceCount << "\n";
		while(index < finalSlices->Count)
		{
			if(finalSlices[index]->_real)
			{
				boundariesFile << finalSlices[index]->_leftBoundary << "\t" << finalSlices[index]->_rightBoundary << "\t"; 
				boundariesFile << finalSlices[index]->_minMolecules << "\t" << finalSlices[index]->_totalMolecules << "\n";
			}
			index++;
		}

		boundariesFile.close();
		boundariesFile.clear();
	}

	//Copy final slices to FSS preselection list
	void CAutoSlice::CopySlices()
	{
        //  FSS - this while loop was added...
		// Build data to return to Preselection
		int index = 0;
		preselectionItems = gcnew List<PreselectionItem ^>();
		while(index < finalSlices->Count /* && finalSlices[index]->_real */)
		{
			// other two variables are filled in later
			PreselectionItem ^ psItem = gcnew PreselectionItem(finalSlices[index]->_leftBoundary, finalSlices[index]->_rightBoundary, 0.0, 0.0, finalSlices[index]->_real);
			if (psItem) 
			{
				preselectionItems->Add(psItem );
			}
			index++;
		}
	}
} // end AutoSlicing namespace