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
	}

	CTraceMap::CTraceMap(String ^ folder)
	{
		_sDirectory = folder;
		xBins = 0;
		yBins = 0;
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
		CLengthProfile ^lengthProfile = gcnew CLengthProfile();
		lengthProfile->Generate2ndorderSmoothedWeightedAverages(lengthHisto, _sDirectory, xBins);

		//Create an instance of CAutoSlice for performing auto-slicing
		CAutoSlice ^autoSlice = gcnew CAutoSlice(_sDirectory);
		if(!autoSlice->PerformAutoSlicing(xBins, lengthHisto, sErr))
		{
			sErr->Format("Something wrong with auto slicing");
			return nullptr; 
		}
		
		//Create an instance of CScoreProfile for score histograms for each slice and estimation
		//the thresholds on the y-axis
		CScoreProfile ^scoreProfile = gcnew CScoreProfile();

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

		return autoSlice->preselectionItems;
	}

	////Write r script for genertating maps with borders
	//bool CTraceMap::WriteRScript(String ^sDirectory, String ^sErr)
	//{
	//	//char line[1000];
	//	//int sliceCount, i = 0, leftBoundary, rightBoundary, minMolecules, maxMolecules;
	//	//double tenPYThreshold, fifteenPYThreshold, TwentyPYThreshold, TwentyFivePYThreshold;
	//	//array<int> ^xBorders;

	//	//String ^fileName = gcnew String("Boundaries.txt");
	//	//fileName = sDirectory + fileName;
	//	//String ^rFileName = gcnew String("Trace2TraceMap-XYBorders.r");
	//	//rFileName = sDirectory + rFileName;
	//	////Convert System::String to char array to read in ofstream
	//	//ifstream BoundariesFile((AutoSlicing::CBasicFunctions::SystemStringToCharArray(fileName)), std::iostream::in);
	//	//ofstream RScriptFile((AutoSlicing::CBasicFunctions::SystemStringToCharArray(rFileName)), std::iostream::out);

	//	//BoundariesFile.getline(line, 1000);
	//	//BoundariesFile >> sliceCount;

	//	////y-Borders - store the borders for each slice
	//	//xBorders = gcnew array<int>(sliceCount*2);
	//	//
	//	//RScriptFile << "#gridmap\n";
	//	//RScriptFile << "library(grid)\n";
	//	//RScriptFile << "library(graphics)\n";
	//	//RScriptFile << "library(lattice)\n";
	//	//RScriptFile << "library(ggplot2)\n";
	//	//RScriptFile << "tracemap <- read.table(\"Trace2TraceMap_Molecules.txt\",header=TRUE)\n";
	//	//RScriptFile << "NumMolecules <- (tracemap$molecules)\n";
	//	//RScriptFile << "png(\"Trace2TraceMap_XYBorders.png\",height=3000,width=4000,res=900, bg=\"white\")\n";
	//	//RScriptFile << "q <- qplot(x=length,y=score,data=tracemap,geom=\"tile\",ylab=\"similarity score\", xlab=\"molecule length (microns)\", fill=log(NumMolecules)) + scale_fill_continuous(low=\"yellow\", high=\"dark red\")\n";

	//	//while(i < sliceCount)
	//	//{
	//	//	BoundariesFile >> leftBoundary >> rightBoundary >> minMolecules >> maxMolecules;
	//	//	RScriptFile << "q <- q + geom_vline(xintercept = " << leftBoundary << ", colour=\"black\", linetype=\"longdash\", size=0.2)\n"; 
	//	//	RScriptFile << "q <- q + geom_vline(xintercept = " << rightBoundary << ", colour=\"red\", linetype=\"longdash\", size=0.2)\n"; 
	//	//	xBorders[2*i] = leftBoundary;
	//	//	xBorders[2*i+1] = rightBoundary;
	//	//	i++;
	//	//}

	//	//i = 0;
	//	//BoundariesFile >> sliceCount;
	//	//while(i < sliceCount)
	//	//{
	//	//	BoundariesFile >> 	tenPYThreshold >> fifteenPYThreshold >> TwentyPYThreshold >> TwentyFivePYThreshold;
	//	//	RScriptFile << "q <- q + geom_segment(aes(x=" << xBorders[2*i] <<",y=" << -1*tenPYThreshold << ", xend=" << xBorders[2*i+1] << ", yend=" << -1*tenPYThreshold << "), size=0.2, linetype=\"longdash\")\n";
	//	//	RScriptFile << "q <- q + geom_segment(aes(x=" << xBorders[2*i] <<",y=" << -1*fifteenPYThreshold << ", xend=" << xBorders[2*i+1] << ", yend=" << -1*fifteenPYThreshold << "), size=0.2, linetype=\"longdash\")\n";
	//	//	RScriptFile << "q <- q + geom_segment(aes(x=" << xBorders[2*i] <<",y=" << -1*TwentyPYThreshold << ", xend=" << xBorders[2*i+1] << ", yend=" << -1*TwentyPYThreshold << "), size=0.2, linetype=\"longdash\")\n";
	//	//	RScriptFile << "q <- q + geom_segment(aes(x=" << xBorders[2*i] <<",y=" << -1*TwentyFivePYThreshold << ", xend=" << xBorders[2*i+1] << ", yend=" << -1*TwentyFivePYThreshold << "), size=0.2, linetype=\"longdash\")\n";
	//	//	i++;
	//	//}
	//	//
	//	//RScriptFile << "grid.newpage()\n";
	//	//RScriptFile << "pushViewport(viewport(height=1.0,width=1.0))\n";
	//	//RScriptFile << "q$grid.fill <- \"white\"\n";
	//	//RScriptFile << "q <- q + theme_bw()\n";
	//	//RScriptFile << "print(q,newpage=FALSE)\n";
	//	//RScriptFile << "upViewport()\n";
	//	//RScriptFile << "dev.off()\n";

	//	return true;
	//}

	 //List<PreselectionItem ^>^ CTraceMap::PreselectionItems()
	 //{
		// return ^ autoSlice->preselectionItems;
	 //}
} //end AutoSlicing namespace
