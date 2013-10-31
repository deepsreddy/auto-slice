#include "stdafx.h"
#include "TraceMap.h"
#include "BasicFunctions.h"
#include "LengthProfile.h"
#include "ScoreProfile.h"
#include "AutoSlice.h"

namespace AutoSlicing {

	//Default constructor
	CTraceMap::CTraceMap(void)
	{
		_sDirectory = "";
		xBins = 0;
		yBins = 0;

		//Creating reference for AutoSliceOptions
		cAutoSliceOptions = gcnew CAutoSliceOptions();
	}

	CTraceMap::CTraceMap(String ^ folder)
	{
		_sDirectory = folder;
		xBins = 0;
		yBins = 0;

		//Creating reference for AutoSliceOptions
		cAutoSliceOptions = gcnew CAutoSliceOptions();
	}

	//Read raw data file - Trace2TraceMap_Molecules.txt
	//********************* Comments for FSS ********************* 
	//If you manage to directly pass the data from Matlab to .NET, re-write this function to 
	//load the values from the arrays passed through Matlab
	//********************* End Comments for FSS *****************
    bool CTraceMap::AssembleRawData(String ^sError, array<double> ^xInterp, array<double> ^yInterp, array<double, 2> ^zInterp)
	{
		xBins = xInterp->GetLength(0);
		yBins = yInterp->GetLength(0);
		
		//Create arrays of rawData struct
		long irawArraySize = xBins*yBins;  
		rawData = gcnew array<CRawData ^>(irawArraySize);

		int irawArrayIndex;
		int yIndex;
		int xIndex;
		for (yIndex = 0; yIndex < yBins; yIndex++)
		{
			for (xIndex = 0; xIndex < xBins; xIndex++)
			{
				irawArrayIndex = (yIndex * xBins) + xIndex;

				if (irawArrayIndex >= irawArraySize)
					irawArrayIndex = irawArraySize - 1;
				rawData[irawArrayIndex] = gcnew CRawData();
				rawData[irawArrayIndex]->_moleculelength = (double)xInterp[xIndex];
				rawData[irawArrayIndex]->_distancescore = -1.0 * (double)yInterp[yIndex];
				rawData[irawArrayIndex]->_moleculedensity = zInterp[yIndex, xIndex];
			}
		}

        WriteTrace2TraceMap( xBins, yBins);

		//Create an instance of CLengthProfile for estimation of weightedAverages
		CLengthProfile ^lengthProfile = gcnew CLengthProfile();
		lengthHisto = gcnew array<LengthHistogram  ^>(xBins);
		if(!lengthProfile->GenerateLengthProfile(lengthHisto, rawData, xBins, yBins, _sDirectory))
		{
			sError->Format("Something wrong with estimating length profiles");
			return false;
		}

		return true;
	}

    // This function is for test only. It writes out trace2Trace data for comparison with PGX's matlab output - input to their experimental program.
    // The text file outputs in this code should all be #defined for tracing or debugging only 
    void CTraceMap::WriteTrace2TraceMap(const int xBins, const int yBins)
    {
        // Write out the Trace2TraceMap data for testing and algorithm refinement

		String ^fileName = gcnew String("Trace2TraceMap.txt");
		fileName = _sDirectory + fileName;	
        fileName = "Trace2TraceMap.txt";
		fileName = _sDirectory + fileName;

		//Convert System::String to char array to read in iostream
		ofstream Trace2TraceMapFile(AutoSlicing::CBasicFunctions::SystemStringToCharArray(fileName));

		if(!Trace2TraceMapFile.is_open())
		{
		//	sError->Format(L"Cannot open raw file: {0}", srawFile);
			return;
		}

		Trace2TraceMapFile << "#numxbins:" << xBins << "\n";
        Trace2TraceMapFile << "#numybins:" << yBins << "\n";
		Trace2TraceMapFile << "Length\tscore\tmoelcules\n";

        Trace2TraceMapFile.precision(6);
        // There is code in the slicing algorithm that discards values < 0, so we save _distancescore (Yvals) as opposite sign
        // they get flipped back before charting, and here, since the user expects them to be negative.
		int xIndex = 0;
		int yIndex = 0;
		int iIndex = 0;
		while(yIndex < xBins)
		{
		    while(xIndex < yBins && iIndex < (xBins*yBins))
		    {
			    Trace2TraceMapFile << rawData[iIndex]->_moleculelength << "\t";
			    Trace2TraceMapFile << (-1.0 * rawData[iIndex]->_distancescore) << "\t";   
			    Trace2TraceMapFile << rawData[iIndex]->_moleculedensity << "\n";

			    xIndex++;
                iIndex++;
            }
			yIndex++;
			xIndex = 0;
		}

		Trace2TraceMapFile.close();
		Trace2TraceMapFile.clear();

    }


	//Evaluate Profiles for length, & score profiles and auto-slicing
	List<PreselectionItem ^>^ CTraceMap::EvaluateProfiles( String ^sErr)
	{
		//Read options
		cAutoSliceOptions->ReadOptions();

		CLengthProfile ^lengthProfile = gcnew CLengthProfile();
		lengthProfile->Generate2ndorderSmoothedWeightedAverages(lengthHisto, _sDirectory, xBins);

		//Create an instance of CAutoSlice for performing auto-slicing
		CAutoSlice ^autoSlice = gcnew CAutoSlice(_sDirectory, cAutoSliceOptions);
		if(!autoSlice->PerformAutoSlicing(xBins, lengthHisto, sErr))
		{
			sErr->Format("Something wrong with auto slicing");
			return nullptr; 
		}
		
		//Create an instance of CScoreProfile for score histograms for each slice and estimation
		//the thresholds on the y-axis
		CScoreProfile ^scoreProfile = gcnew CScoreProfile(cAutoSliceOptions);

		//y-Borders - store the borders for each slice for the 15% threshold values
//		scoreProfile->yBorders = gcnew array<double>(autoSlice->finalSlices->Count);

		//For each slice, estimate the threshold along the y-axis
		int iIndex = 0, sliceNumber = 0;
		while(iIndex < autoSlice->finalSlices->Count)
		{
			if(autoSlice->finalSlices[iIndex]->_real)
			{
				double dLeftBoundary = (int) autoSlice->finalSlices[iIndex]->_leftBoundary;
				double dRightBoundary = (int) autoSlice->finalSlices[iIndex]->_rightBoundary;

				//Invoke GenerateScoreProfile for calculating the score profiles
				if(!scoreProfile->GenerateScoreProfile(rawData, xBins, yBins, _sDirectory, dLeftBoundary, dRightBoundary, (sliceNumber+1)))
				{
					sErr->Format("Something wrong with estimating score profiles");
					return nullptr;	
				}

				// For now, CutoffDistanceMin is defaulting to -2
				autoSlice->preselectionItems[iIndex]->CutoffDistanceValue = -1 * (scoreProfile->yBorders);
				sliceNumber++;
			}
			iIndex++;
		}

        // Write Boundaries file
		String ^fileName = gcnew String("Boundaries.txt");
		fileName = _sDirectory + fileName;

        iIndex = 0;
		//Convert System::String to char array to read in ofstream
		ofstream xBoundariesFile((AutoSlicing::CBasicFunctions::SystemStringToCharArray(fileName)), std::iostream::app);
		xBoundariesFile << autoSlice->finalSlices->Count << "\n";

        int realCount = 0;
		while(iIndex < autoSlice->finalSlices->Count)
		{
			if(autoSlice->finalSlices[iIndex]->_real)
			{
                xBoundariesFile << scoreProfile->yBorders  << "\n";
                realCount++;
            }
            iIndex++;
        }
        xBoundariesFile << realCount << "\n";

		xBoundariesFile.close();
		xBoundariesFile.clear();

		//Write settings in XML file
		cAutoSliceOptions->WriteOptions();

		return autoSlice->preselectionItems;
	}

	 //List<PreselectionItem ^>^ CTraceMap::PreselectionItems()
	 //{
		// return ^ autoSlice->preselectionItems;
	 //}
} //end AutoSlicing namespace
