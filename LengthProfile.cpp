#include "stdafx.h"
#include "LengthProfile.h"
#include "BasicFunctions.h"
#include "sgfilter.h"

namespace AutoSlicing {

	//Default constructor
	CLengthProfile::CLengthProfile(void)
	{

	}

	//Estimate Length Profile: calculates the plain summation and weighted average of the length histogram
	bool CLengthProfile::GenerateLengthProfile(array<LengthHistogram ^> ^lengthHisto, array<CRawData ^> ^vRawData, const int xBins, int yBins, String ^sDirectory)
	{
		long idataSize = xBins*yBins;

		//Initialize the lengthHisto[] array
		for(long i=0; i < xBins; i++)
		{
			lengthHisto[i] = gcnew LengthHistogram();
			lengthHisto[i]->_XCoordinate = vRawData[i]->_moleculelength;
			lengthHisto[i]->_densitySummation = 0.0;

			lengthHisto[i]->_weightedDensityAverage = 0.0;
			lengthHisto[i]->_smoothedWeightedAverage = 0.0;
			lengthHisto[i]->_PKFlag = '0';
		}

		double sumDistanceScore;
		//Summation and weighted average - read from the vRawData array passed from CTraceMap
		for(long i=0; i < xBins; i++)
		{
			long j;
			try
			{
				sumDistanceScore = 0.0;
				for(j=i; j < idataSize; j+=xBins)
				{
					try
					{
						lengthHisto[i]->_densitySummation += vRawData[j]->_moleculedensity;
						lengthHisto[i]->_weightedDensityAverage += (vRawData[j]->_distancescore*vRawData[j]->_moleculedensity);
						sumDistanceScore += vRawData[j]->_distancescore;
					}
					catch (Exception^ )
					{
						return false;
					}
				}
			}
			catch (Exception^ )
			{
				return false;
			}
			
			if(sumDistanceScore > 0.0)
				lengthHisto[i]->_weightedDensityAverage = (lengthHisto[i]->_weightedDensityAverage/sumDistanceScore);
			
		}
		///TODO FSS
		//Output of weighted averages for applying smoothing
		//********************* Comments for FSS ********************* 
		//Output of the file weightedAverages.txt is not needed if we directly pass the double array to the function "sgfilter" in SGFilter.cpp
		//********************* End Comments for FSS *****************

		String ^fileName = gcnew String("lengthProfile.txt");
		String ^weightedFileName = gcnew String("weightedAverages.txt");
		fileName = sDirectory + fileName;
		weightedFileName = sDirectory + weightedFileName;
		Console::WriteLine(L"\nOutput filename: {0}", fileName);
		ofstream lengthProfileFile(AutoSlicing::CBasicFunctions::SystemStringToCharArray(fileName));
		ofstream weightedProfileFile(AutoSlicing::CBasicFunctions::SystemStringToCharArray(weightedFileName));
		lengthProfileFile << "Molecule_Length\tDensity_Summation\tWeighted_Density_Average\n";
		
		for(long i=0; i < xBins; i++)
		{
			lengthProfileFile << lengthHisto[i]->_XCoordinate << "\t" << abs(lengthHisto[i]->_densitySummation) << "\t" << lengthHisto[i]->_weightedDensityAverage << "\n";
			weightedProfileFile << lengthHisto[i]->_XCoordinate << "\t" << lengthHisto[i]->_weightedDensityAverage << "\n";
		}
		
		lengthProfileFile.close();
		lengthProfileFile.clear();
		weightedProfileFile.close();
		weightedProfileFile.clear();

        if(!Generate2ndorderSmoothedWeightedAverages (lengthHisto, sDirectory, xBins))
		{
			Console::WriteLine("Smoothed Weighted Averages error");
            return false;
		}

		return true;
	}

	bool CLengthProfile::Generate2ndorderSmoothedWeightedAverages (array<LengthHistogram ^> ^lengthHisto, String ^sDirectory, int xBins)
	{
		//The smoothed weightedAverages is estimated by running it through the Savitzky-Golay function (SGFilter)
		if(!SmoothingInterface(lengthHisto, xBins))
		{
			Console::WriteLine("Smoothing error");
            return false;
		}


        // Write out the smoothed data for testing and algorithm refinement

		String ^fileName = gcnew String("2ndorderSmoothedWeightedAverages.txt");
		fileName = sDirectory + fileName;	

		//Convert System::String to char array to read in ifstream
		ofstream smoothedAvgFile(AutoSlicing::CBasicFunctions::SystemStringToCharArray(fileName));

		if(!smoothedAvgFile.is_open())
		{
			//sError->Format(L"Cannot open raw file: {0}", srawFile);
			return false;
		}

		smoothedAvgFile << "#numxbins:" << xBins << "\n";
		smoothedAvgFile << "MoleculeLength\twAverage\tsFrame5\n";
		int iIndex = 0;


		while(iIndex < xBins)
		{
			smoothedAvgFile << lengthHisto[iIndex]->_XCoordinate << "\t";
			smoothedAvgFile << lengthHisto[iIndex]->_weightedDensityAverage << "\t";
			smoothedAvgFile << lengthHisto[iIndex]->_smoothedWeightedAverage << "\n";
			iIndex++;
		}

		smoothedAvgFile.close();
		smoothedAvgFile.clear();

		String ^sErr;
		WriteRScript(sDirectory, sErr);

		return true;
		
	}

    //Marshalling - interfacing with native code for smoothing in SGFilter.cpp
    // Update LengthHistogram with smothed values
    bool CLengthProfile::SmoothingInterface(array<LengthHistogram ^> ^lengthHisto, const int xBins)
	{
        int iIndex = 0;
        //Managed double array for weighted averages
        array<double> ^weightedAvg = gcnew array<double>(xBins);
        array<double> ^smoothedAvg = gcnew array<double>(xBins);
        double smoothedArray[1000];

		//Copy the weighted Averages in lengthHisto to the double array
        while(iIndex < xBins)
		{
			weightedAvg[iIndex] = lengthHisto[iIndex]->_weightedDensityAverage;
			smoothedAvg[iIndex] = 0.0;
            iIndex++;
		}
		//Marshal to native double pointer
		pin_ptr<double> pin = &weightedAvg[0];

        //Call to native SGFilter - returns '0' for successful call and non-'0' for unsuccessful call
		if(sgfilter_main(pin, smoothedArray, xBins, LEFTFRAMESIZE, RIGHTFRAMESIZE, ORDEROFDERIVATIVE, POLYNOMIALORDER))
		{
			Console::WriteLine(L"Error in smoothing function");
            return false;
		}

        //Marshal native array to managed array 
		Marshal::Copy(static_cast<IntPtr>(smoothedArray),smoothedAvg, 0, (xBins));

        iIndex = 0;
        while(iIndex < xBins)
		{
			lengthHisto[iIndex]->_smoothedWeightedAverage = smoothedAvg[iIndex];
            iIndex++;
		}

		return true;
	}

	//Write r script for genertating lenght histograms - weighted and smoothed
	bool CLengthProfile::WriteRScript(String ^sDirectory, String ^sErr)
	{
		String ^rFileName = gcnew String("SmoothedProfiles.r");
		rFileName = sDirectory + rFileName;
		//Convert System::String to char array to read in ofstream
		ofstream RScriptFile((AutoSlicing::CBasicFunctions::SystemStringToCharArray(rFileName)), std::iostream::out);

		RScriptFile << "#gridmap\n";
		RScriptFile << "library(grid)\n";
		RScriptFile << "library(graphics)\n";
		RScriptFile << "library(lattice)\n";
		RScriptFile << "library(ggplot2)\n";
		RScriptFile << "smoothProfile <- read.table(\"2ndorderSmoothedWeightedAverages_cpp.txt\",header=TRUE)\n";
		RScriptFile << "png(\"2ndOrderSmoothedProfilePlot.png\",height=3000,width=4000,res=900, bg=\"white\")\n";
		RScriptFile << "q <- qplot(x=MoleculeLength,y=wAverage,data=smoothProfile,geom=\"line\",ylab=\"normalized density\", xlab=\"molecule length (microns)\")\n";
		RScriptFile << "q <- q + geom_line(aes(y=smoothProfile$sFrame5, x=MoleculeLength), color=\"red\", size=0.4)\n";
		RScriptFile << "grid.newpage()\n";
		RScriptFile << "pushViewport(viewport(height=1.0,width=1.0))\n";
		RScriptFile << "q$grid.fill <- \"white\"\n";
		RScriptFile << "q <- q + theme_bw()\n";
		RScriptFile << "print(q,newpage=FALSE)\n";
		RScriptFile << "upViewport()\n";
		RScriptFile << "dev.off()\n";

		RScriptFile.close();
		RScriptFile.clear();

		return true;
	}
} // end of namespace AutoSlicing
