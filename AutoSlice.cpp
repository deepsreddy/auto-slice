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
	CAutoSlice::CAutoSlice(String ^sDirectory)
	{
		_sDirectory = sDirectory;
		_UpperSliceMolecules = 0.0;
		_LowerSliceMolecules = 0.0;
		_iLowerSliceCount = 0;
		_iUpperSliceCount = 0;
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
		int iCurrentPosition = 0;

		profileProperties^ initialProfile = gcnew profileProperties();

		//  FSSTODO JEP This next loop was added back after we removed something similar to it during the initial port
		while((dCurrentSignal = lengthHisto1[iCurrentPosition]->_smoothedWeightedAverage) <= VERY_SMALL)
			iCurrentPosition++;

		//Go to the previous '0' molecules position
		if(iCurrentPosition > 0)
			iCurrentPosition--;
		initialProfile->_dCurrentSignal = lengthHisto1[iCurrentPosition]->_smoothedWeightedAverage;
		initialProfile->_dNextSignal = lengthHisto1[iCurrentPosition+1]->_smoothedWeightedAverage;
//		initialProfile->_dPreviousSignal = -1;
//		initialProfile->_iPosition = iCurrentPosition;
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
			dNextSignal = lengthHisto1[iCurrentPosition+1]->_smoothedWeightedAverage;
			dPreviousSignal = lengthHisto1[iCurrentPosition-1]->_smoothedWeightedAverage;

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
			iCurrentPosition++;
		}

		//Go back to the last position
		iCurrentPosition--;

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


		Console::WriteLine(L"Number of peak points: {0}", peakProfiles->Count);
		Console::WriteLine(L"Number of trough points: {0}", troughProfiles->Count);

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

	//Clean profiles for adding slices
	//The name of this function might be misleading. This function basically creates preliminary slices by taking trough-peak-trough positions
	//and calculates the trough1/peak & trough2/peak ratios. Each of these preliminary slices is populated in the list <sliceProfiles>
	bool CAutoSlice::EstimateTrough2PeakRatios(int iLengthHistoSize, array<LengthHistogram ^> ^lengthHisto1, String ^sErr)
	{
		//sliceProfiles holds the pseudo-slices with ratios and other properties
		sliceProfiles = gcnew List<sliceProperties ^>();
		int iCount = peakProfiles->Count;
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

			if(tempSlice->_trough2Position >= BOUNDARY_LENGTH)
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
			if(sliceProfiles[i]->_trough2Position >= BOUNDARY_LENGTH)
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

//		autoBoundaryDetection(iLengthHistoSize, lengthHisto1, sErr);

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

		while(iCurrentSlice < (iNumberofSlices))
		{
			finalSliceProperties ^tempSlice = gcnew finalSliceProperties();
			leftBoundary = sliceProfiles[iCurrentSlice]->_trough1Position;

			if(sliceProfiles[iCurrentSlice]->_trough2Position > BOUNDARY_LENGTH)
			{
				dNumberMoleculesThreshold = NUMBEROFMOLECULES_UPPERLENGTH;
				dMinPercentMolecules = MINIMUM_UPPERPERCENTMOLECULES;
				ratioThreshold = TROUGHTOPEAKRATIO_UPPERLENGTH;
			}
			else
			{
				dNumberMoleculesThreshold = NUMBEROFMOLECULES_LOWERLENGTH;
				dMinPercentMolecules = MINIMUM_LOWERPERCENTMOLECULES;
				ratioThreshold = TROUGHTOPEAKRATIO_LOWERLENGTH;
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
//			sliceMode = check4Slices(iCurrentSlice, iSliceCount);
//			if(sliceMode == ADD)
//			{
			while((iCurrentSlice+1) < iNumberofSlices && (sliceProfiles[iCurrentSlice]->_trough2peakratio2 > ratioThreshold))
			{
				if(sliceProfiles[iCurrentSlice]->_trough2Position > BOUNDARY_LENGTH)
					ratioThreshold = TROUGHTOPEAKRATIO_UPPERLENGTH;
				else
					ratioThreshold = TROUGHTOPEAKRATIO_LOWERLENGTH;

				
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
			if(tempSlice->_rightBoundary > BOUNDARY_LENGTH)
			{
				dNumberMoleculesThreshold = NUMBEROFMOLECULES_UPPERLENGTH;
				dMinPercentMolecules = MINIMUM_UPPERPERCENTMOLECULES;
				tempSlice->_percentMolecules = (totalMolecules/_UpperSliceMolecules)*100;
			}
			else
			{
				dNumberMoleculesThreshold = NUMBEROFMOLECULES_LOWERLENGTH;
				dMinPercentMolecules = MINIMUM_LOWERPERCENTMOLECULES;
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

//			tempSlice->_real = true;
//			finalSlices->Add(tempSlice);
//			iSliceCount++;
////			}
//			iCurrentSlice++;
//		}


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


		index = 0;

		String ^fileName = gcnew String("Boundaries.txt");
		fileName = _sDirectory + fileName;

		//Convert System::String to char array to read in ofstream
		ofstream boundariesFile((AutoSlicing::CBasicFunctions::SystemStringToCharArray(fileName)), std::iostream::out);
		boundariesFile << "left_boundary\tright_boundary\n";
		boundariesFile << iSlices << "\n";
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

        //  FSS - this while loop was added...
		// Build data to return to Preselection
		index = 0;
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
		return true;
	}

#if 0
	//Check for each slide to add or discard, returns '1' for 'Add', '0' for 'Discard' 
	//Please note the code needs to be cleaned a bit
	bool CAutoSlice::check4Slices(int &iCurrentSlicePosition, int &iSliceCount)
	{
		//default set to DISCARD the slice
		double ratio1 = sliceProfiles[iCurrentSlicePosition]->_trough2peakratio1;
		double ratio2 = sliceProfiles[iCurrentSlicePosition]->_trough2peakratio2;
		double dTrough2PeakRatioThreshold;
		double xPosition = sliceProfiles[iCurrentSlicePosition]->_peakPosition;
		bool bLastSlice = false, bMode;
		
		if((iCurrentSlicePosition+1) == (sliceProfiles->Count))
			bLastSlice = true;

		if(iSliceCount > 0)
		{
			if(finalSlices[iSliceCount-1]->_rightBoundary > BOUNDARY_LENGTH)
				dTrough2PeakRatioThreshold = TROUGHTOPEAKRATIO_UPPERLENGTH;
			else
				dTrough2PeakRatioThreshold = TROUGHTOPEAKRATIO_LOWERLENGTH;
		}
		else
			dTrough2PeakRatioThreshold = TROUGHTOPEAKRATIO_LOWERLENGTH;

		if(ratio1 <= dTrough2PeakRatioThreshold)
		{
			if(ratio2 <= dTrough2PeakRatioThreshold)
			{
				if(!bLastSlice)
				{
					if(sliceProfiles[iCurrentSlicePosition+1]->_trough1Position > BOUNDARY_LENGTH)
						dTrough2PeakRatioThreshold = TROUGHTOPEAKRATIO_UPPERLENGTH;
					else
						dTrough2PeakRatioThreshold = TROUGHTOPEAKRATIO_LOWERLENGTH;

					if(sliceProfiles[iCurrentSlicePosition+1]->_trough2peakratio1 > dTrough2PeakRatioThreshold)
					{
						if(sliceProfiles[iCurrentSlicePosition+1]->_trough2peakratio2 > dTrough2PeakRatioThreshold)
						{
							//Check for junk slice
							if(!check4CurrentJunkSlice(iCurrentSlicePosition, (iSliceCount-1), dTrough2PeakRatioThreshold))
								bMode = check4NextJunkSlice(iCurrentSlicePosition, sliceProfiles[iCurrentSlicePosition]->_totalMolecules);
							else
								return DISCARD;
							
							if(bMode)
								return ADD;
							else
							{
								//Advance to slices until the next slice is part of the current one
								bMode = advanceSlices(iCurrentSlicePosition);
								if(bMode)
									return ADD;
							}
						}
						else
						{
							//Confirm that its not junk slice
							if(!check4CurrentJunkSlice(iCurrentSlicePosition, (iSliceCount-1), dTrough2PeakRatioThreshold))
								bMode = check4NextJunkSlice(iCurrentSlicePosition, sliceProfiles[iCurrentSlicePosition]->_totalMolecules);
							else
								return DISCARD;

							if(bMode == true)
								return ADD;
							else
							{
								//Advance to next slice which is part of the current slice
								iCurrentSlicePosition++;
								return ADD;
							}
						}
					}
					else
					{
						if(!check4CurrentJunkSlice(iCurrentSlicePosition, (iSliceCount-1), dTrough2PeakRatioThreshold))
							return ADD;
					}
				}
				else
				{
					//Check if the last slice is junk
					if(!check4CurrentJunkSlice(iCurrentSlicePosition, (iSliceCount-1), dTrough2PeakRatioThreshold))
						bMode = check4NextJunkSlice(iCurrentSlicePosition, sliceProfiles[iCurrentSlicePosition]->_totalMolecules);
					else
						return DISCARD;

					if(bMode == true)
						return DISCARD;
					else
						return ADD;
				}
			}
			else	
			{
				if(!bLastSlice)
				{
					//check for junk slice
					if(!check4CurrentJunkSlice(iCurrentSlicePosition, (iSliceCount-1), dTrough2PeakRatioThreshold))
						bMode = check4NextJunkSlice(iCurrentSlicePosition, sliceProfiles[iCurrentSlicePosition]->_totalMolecules);
					else
						return DISCARD;

					if(bMode == true)
						return DISCARD;

					//sliceWidth = (sliceProfiles[iCurrentSlicePosition]->_trough2Position-sliceProfiles[iCurrentSlicePosition]->_trough1Position);
					//if( (sliceWidth < 5) || sliceProfiles[iCurrentSlicePosition]->_totalMolecules < MIN_NUMBEROFMOLECULES)
						//return DISCARD;
			
					//Advance to where the next slice is part of the current slice
					bMode = advanceSlices(iCurrentSlicePosition);

					if(bMode)
						return ADD;
				}
				else
				{
					//Check if the last slice is junk
					if(!check4CurrentJunkSlice(iCurrentSlicePosition, (iSliceCount-1), dTrough2PeakRatioThreshold))
						bMode = check4NextJunkSlice(iCurrentSlicePosition, sliceProfiles[iCurrentSlicePosition]->_totalMolecules);
					else
						return DISCARD;

					if(bMode == true)
						return DISCARD;
					else
						return ADD;
				}
			}
		}
		else
		{
			//Advance until we get <LIMIT ratio1 & possibly discard very thin slices with less molecules
			//if(ratio2 <= TROUGHTOPEAKRATIO_LOWERLENGTH)
			{
				//sliceWidth = (sliceProfiles[iCurrentSlicePosition]->_trough2Position-sliceProfiles[iCurrentSlicePosition]->_trough1Position);
				if(!bLastSlice)
				{
					//check for junk slice
					if(!check4CurrentJunkSlice(iCurrentSlicePosition, (iSliceCount-1), dTrough2PeakRatioThreshold))
						bMode = check4NextJunkSlice(iCurrentSlicePosition, sliceProfiles[iCurrentSlicePosition]->_totalMolecules);
					else
						return DISCARD;

					if(bMode == true)
						return DISCARD;

					//Advance to where the next slice is part of the current slice
					bMode = advanceSlices(iCurrentSlicePosition);

					if(bMode)
						return ADD;
				}
				else
				{
					//Check if the last slice is junk
					if(!check4CurrentJunkSlice(iCurrentSlicePosition, (iSliceCount-1),  dTrough2PeakRatioThreshold))
						bMode = check4NextJunkSlice(iCurrentSlicePosition, sliceProfiles[iCurrentSlicePosition]->_totalMolecules);
					else
						return DISCARD;

					if(bMode == true)
						return DISCARD;
					else
						return ADD;
				}
			}
		}
	
		return DISCARD;
	}

	//Check for current slice if its junk or not
	bool CAutoSlice::check4CurrentJunkSlice(int &iCurrentSlicePosition, int iSlicePosition, double dRatioThreshold)
	{
		double dPercentMoleculesThreshold;
		double dNumberMoleculesThreshold, ratio, sliceWidth;

		sliceWidth = (sliceProfiles[iCurrentSlicePosition]->_trough2Position-sliceProfiles[iCurrentSlicePosition]->_trough1Position);

		if(iSlicePosition <= 0 || iCurrentSlicePosition <=0)
			return false;

		///TODO check this - looks wrong...

		//if(sliceProfiles[iCurrentSlicePosition]->_trough1Position > BOUNDARY_LENGTH)
		if(finalSlices[iSlicePosition]->_rightBoundary > BOUNDARY_LENGTH)
		{
			dPercentMoleculesThreshold = PERCENTMOLECULES_UPPERLENGTH;
			dNumberMoleculesThreshold = NUMBEROFMOLECULES_UPPERLENGTH;
		}
		else
		{
			dPercentMoleculesThreshold = PERCENTMOLECULES_LOWERLENGTH;
			dNumberMoleculesThreshold = NUMBEROFMOLECULES_LOWERLENGTH;
		}

		if(finalSlices[iSlicePosition]->_rightBoundary > BOUNDARY_LENGTH)
			dPercentMoleculesThreshold = PERCENTMOLECULES_UPPERSLICETHRESHOLD;	// this overwrites what was just done above
		else
			dPercentMoleculesThreshold = PERCENTMOLECULES_LOWERSLICETHRESHOLD;
	
		ratio = (sliceProfiles[iCurrentSlicePosition]->_totalMolecules/finalSlices[iSlicePosition]->_minMolecules);

		if(sliceProfiles[iCurrentSlicePosition]->_trough2peakratio1 >= dRatioThreshold && sliceProfiles[iCurrentSlicePosition]->_trough2peakratio2 >= dRatioThreshold)
        {
            //  FSS - the line commented out was added by FSS to replace the next line - which is correct?
            //if(ratio < dPercentMoleculesThreshold)
			if(ratio < dPercentMoleculesThreshold || sliceWidth < 5.0)
				return true;
			else if(sliceProfiles[iCurrentSlicePosition]->_totalMolecules <= dNumberMoleculesThreshold)
				return true;
		}
		else if(sliceProfiles[iCurrentSlicePosition]->_trough2peakratio1 >= dRatioThreshold || sliceProfiles[iCurrentSlicePosition]->_trough2peakratio2 >= dRatioThreshold)
		{
			if(ratio < dPercentMoleculesThreshold)
				return true;
			else if(sliceProfiles[iCurrentSlicePosition]->_totalMolecules <= dNumberMoleculesThreshold)
				return true;
		}

		//ratio = (sliceProfiles[iCurrentSlicePosition]->_totalMolecules/sliceProfiles[iCurrentSlicePosition-1]->_totalMolecules);

		///TODO check this logic!

		if(ratio < dPercentMoleculesThreshold)
				return true;
		else if(sliceProfiles[iCurrentSlicePosition]->_totalMolecules < dNumberMoleculesThreshold)
			return true;
		else if(sliceWidth < 5.0 && ratio < dPercentMoleculesThreshold)	// can't get here? (already returned true above...)
			return true;
		else 
			return false;
	}

	//Check for if the next slice is junk or not
	bool CAutoSlice::check4NextJunkSlice(int &iCurrentSlicePosition, int iNumberMoleculesComparison)
	{
		double dPercentMoleculesThreshold;
        double dNumberMoleculesThreshold, ratio;
        //  FSS - added the following variable and all code related to it
        double dPercentMoleculesJunkThreshold;

		if(iCurrentSlicePosition < (sliceProfiles->Count-1))
		{
			if(sliceProfiles[iCurrentSlicePosition+1]->_trough1Position > BOUNDARY_LENGTH)
			{
				dPercentMoleculesThreshold = PERCENTMOLECULES_UPPERLENGTH;
				dPercentMoleculesJunkThreshold = PERCENTHIGHERLENGTHMOLECULES_JUNKTHRESHOLD;
			}
			else
			{
				dPercentMoleculesThreshold = PERCENTMOLECULES_LOWERLENGTH;
				dPercentMoleculesJunkThreshold = PERCENTLOWERLENGTHMOLECULES_JUNKTHRESHOLD;
			}
		}

		if(sliceProfiles[iCurrentSlicePosition]->_trough1Position > BOUNDARY_LENGTH)
			dNumberMoleculesThreshold = NUMBEROFMOLECULES_UPPERLENGTH;
		else
			dNumberMoleculesThreshold = NUMBEROFMOLECULES_LOWERLENGTH;
	
		//Check for junk slices
		if(iCurrentSlicePosition < (sliceProfiles->Count-1))
		{
			//ratio = (sliceProfiles[iCurrentSlicePosition+1]->_totalMolecules/sliceProfiles[iCurrentSlicePosition]->_totalMolecules);
			ratio = (sliceProfiles[iCurrentSlicePosition+1]->_totalMolecules/iNumberMoleculesComparison);

			if(ratio < dPercentMoleculesThreshold)
                return true;
            //  FSS - commented out the following two lines
			//else if(sliceProfiles[iCurrentSlicePosition+1]->_totalMolecules < dNumberMoleculesThreshold)
			//	return true;
			else if(sliceProfiles[iCurrentSlicePosition+1]->_trough2peakratio1 >= TROUGHTOPEAKRATIO_LOWERLENGTH && sliceProfiles[iCurrentSlicePosition+1]->_trough2peakratio2 >= TROUGHTOPEAKRATIO_LOWERLENGTH)
			{
				//if(ratio <= PERCENTMOLECULES_LOWERLENGTH)
				if(ratio <= dPercentMoleculesJunkThreshold)
					return true;
				else
					return false;
			}
			else
				return false;
		}
        //  FSS - commented out the following two lines
		//else if(sliceProfiles[iCurrentSlicePosition]->_totalMolecules < dNumberMoleculesThreshold)
		//	return true;
		else
			return false;
	}

	//Advance slices until they are part of the current one
	bool CAutoSlice::advanceSlices(int &iCurrentSlicePosition)
	{
		double dTrough2PeakRatioThreshold;
		bool bIncrement = false, bMode = false;
		int iCurrSlice = iCurrentSlicePosition;
		Console::WriteLine(L"Entering advanceSlices() - iCurrentSlicePosition: {0}", iCurrentSlicePosition);

		if(sliceProfiles[iCurrentSlicePosition+1]->_trough1Position > BOUNDARY_LENGTH)
			dTrough2PeakRatioThreshold = TROUGHTOPEAKRATIO_UPPERLENGTH;
		else
			dTrough2PeakRatioThreshold = TROUGHTOPEAKRATIO_LOWERLENGTH;

		while(iCurrentSlicePosition < (sliceProfiles->Count - 1) && sliceProfiles[iCurrentSlicePosition+1]->_trough2peakratio2 > dTrough2PeakRatioThreshold)
		{
			bMode = check4NextJunkSlice(iCurrentSlicePosition, sliceProfiles[iCurrSlice]->_totalMolecules);
			if(!bMode)
				iCurrSlice++;
			iCurrentSlicePosition++;
		}
		///TODO May want to remove one of these? Increments above in While loop then again here - skips a slice?
		//if(bIncrement)
		//  FSSTODO - JEP 7/8/2013 - added the following conditional to avoid an exception after returning to the calling function
		if(iCurrentSlicePosition < (sliceProfiles->Count - 1))
		{
			iCurrentSlicePosition++;
		}
		// Console::WriteLine(L"Skipped incrementing iCurrentSlicePosition");
		Console::WriteLine(L" leaving advanceSlices() - iCurrentSlicePosition: {0}", iCurrentSlicePosition);

		return true;
	}
#endif
}