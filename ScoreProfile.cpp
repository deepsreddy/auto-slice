#include "stdafx.h"
#include "ScoreProfile.h"
#include "BasicFunctions.h"

namespace AutoSlicing 
{

	//Default constructor
	CScoreProfile::CScoreProfile(void)
	{
	}

	//Estimate the score Profile: calculates the plain summation of score values in each slice defined by the left and right boundary in the length axis
	bool CScoreProfile::GenerateScoreProfile(array<CRawData ^> ^vRawData, int xBins, int yBins, String ^sDirectory, double leftBoundary, double rightBoundary, int sliceNumber)
	{
		long idataSize = xBins*yBins;
		scoreHisto = gcnew array<ScoreHistogram  ^>(yBins);

		long rawIndex = 0;
		for(long i=0; i < yBins; i++)
		{
			scoreHisto[i] = gcnew ScoreHistogram();
			scoreHisto[i]->_YCoordinate = vRawData[rawIndex]->_distancescore;
			scoreHisto[i]->_densitySummation = 0.0;	

			rawIndex += xBins;
		}

		_totalMolecules = 0.0;
		//_YMin, _YMax holds the minimum and maximum y-score values
		_YMin = 1000.00;
		_YMax = -10.0;

		//Summation of the score values defined by the left and right boundary values
		for(long i=0; i < yBins; i++)
		{
			rawIndex = i*xBins;


			for(long j=0; j < xBins; j++)
			{
				if(vRawData[rawIndex]->_moleculelength >= leftBoundary && vRawData[rawIndex]->_moleculelength <= rightBoundary)
				{
					scoreHisto[i]->_densitySummation += vRawData[rawIndex]->_moleculedensity;
					_totalMolecules += vRawData[rawIndex]->_moleculedensity;
					if(vRawData[rawIndex]->_moleculedensity > 0.0 && scoreHisto[i]->_YCoordinate < _YMin)
						_YMin = scoreHisto[i]->_YCoordinate;
					if(vRawData[rawIndex]->_moleculedensity > 0.0 && scoreHisto[i]->_YCoordinate > _YMax)
						_YMax = scoreHisto[i]->_YCoordinate;
				}
				rawIndex++;
			}
		}
		// Removed commented code here

		//Output
		//Outputting the score summation (histogram) values for all the slices
		//FSS: Output of these histograms are not needed
		double avgDensity = 0;
		///TODO Take out file IO, or #DEFINE it for testing only
		string fileName = "scoreProfile-Slice";
		string ssliceNumber = std::to_string((long long)sliceNumber);
		fileName = AutoSlicing::CBasicFunctions::SystemStringToCharArray(sDirectory) + fileName + ssliceNumber + ".txt";

		ofstream scoreProfileFile(fileName);
		scoreProfileFile << "Similarity Score\tDensity Summation\tDensity Average\n";
		for(long i=0; i < yBins; i++)
		{
			avgDensity = (scoreHisto[i]->_densitySummation/(double) xBins);
			scoreProfileFile << scoreHisto[i]->_YCoordinate << "\t" << abs(scoreHisto[i]->_densitySummation) << "\t" << avgDensity << "\n";
		}
		scoreProfileFile.close();
		scoreProfileFile.clear();

//		if (sliceNumber-1 >= 0 )
//		{
//			yBorders[sliceNumber-1] = FindYThreshold(yBins, 0.15); 
//		}
		//Find y-score threshold value corresponding to 6% or 15% (YSCORE_THRESHOLD) of the area within the y-axis density summation
		//which is dependent on the slice length range
		if(rightBoundary >= BOUNDARY_LENGTH)
			yBorders = FindYThreshold(yBins, YSCORE_HIGHERTHRESHOLD);
		else
			yBorders = FindYThreshold(yBins, YSCORE_LOWERTHRESHOLD);

		return true;
	}


	//Find the thresholding borders for each slice
	double CScoreProfile::FindYThreshold(int yBins, double dThresholdPercent)
	{
		double yBorder, sumMolecules = 0.0, prevSum = 0.0;

		// The slice is examined from the top.
		// When 15% of the total number of molecules have been accumulated, that Y coordinate is used as the preliminary upper bound of the slice 
		double _ThresholdMolecules = _totalMolecules*dThresholdPercent;

		double prevPeakPosition, currPeakPosition;
		
		int i = (yBins-1);
		sumMolecules += scoreHisto[i]->_densitySummation;
		prevSum = sumMolecules;
		i--;
		while(sumMolecules <= _ThresholdMolecules)
		{
			prevSum = sumMolecules;
			sumMolecules += scoreHisto[i]->_densitySummation;
			i--;
		}

		i++;
		prevPeakPosition = scoreHisto[i+1]->_YCoordinate;
		currPeakPosition = scoreHisto[i]->_YCoordinate;

		//Calculate the minimum y-score at which the percent of molecules fall below certain threshold
		// of the total area in the region defined by left and right boundaries.
		//The threshold is passed on as a parameter (dThresholdPercent) to this function
		yBorder = prevPeakPosition + (((currPeakPosition-prevPeakPosition)*(_ThresholdMolecules-prevSum))/(sumMolecules-prevSum));

		//Rounding to two decimal points
		return ceil(yBorder * 100) / 100;
	}

} // end of namespace AutoSlicing
