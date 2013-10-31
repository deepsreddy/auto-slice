#include "stdafx.h"
#include "AutoSliceOptions.h"
using namespace System::IO;
using namespace System::Text::RegularExpressions;

namespace AutoSlicing {

	CAutoSliceOptions::CAutoSliceOptions(void)
	{
		
	}
	//=============================================================================================================
	//Setting default values if the settings file is not loaded
	void CAutoSliceOptions::operator=(const CAutoSliceOptions ^sXMLOptions)
	{
		_dTroughtoPeakRatioLowerLengthLimit = sXMLOptions->_dTroughtoPeakRatioLowerLengthLimit;
		_dTroughtoPeakRatioUpperLengthLimit = sXMLOptions->_dTroughtoPeakRatioUpperLengthLimit;
		_dMinimumNumberofMoleculesLowerLength = sXMLOptions->_dMinimumNumberofMoleculesLowerLength;
		_dMinimumNumberofMoleculesUpperLength = sXMLOptions->_dMinimumNumberofMoleculesUpperLength;
		_dBoundaryLength = sXMLOptions->_dBoundaryLength;
		_dMinimumLowerPercentMolecules = sXMLOptions->_dMinimumLowerPercentMolecules;
		_dMinimumUpperPercentMolecules = sXMLOptions->_dMinimumUpperPercentMolecules;
		_dYScoreLowerLengthThreshold = sXMLOptions->_dYScoreLowerLengthThreshold;
		_dYScoreUpperLengthThreshold = sXMLOptions->_dYScoreUpperLengthThreshold;
	}
	//=============================================================================================================
	//Setting default values if the settings file is not loaded
	void CAutoSliceOptions::SetDefaultParameters(void)
	{
		_dTroughtoPeakRatioLowerLengthLimit = 0.75;
		_dTroughtoPeakRatioUpperLengthLimit = 0.55;
		_dMinimumNumberofMoleculesLowerLength = 100.0;
		_dMinimumNumberofMoleculesUpperLength = 50.0;
		_dBoundaryLength = 80.0;
		_dMinimumLowerPercentMolecules = 3.0;
		_dMinimumUpperPercentMolecules = 0.5;
		_dYScoreLowerLengthThreshold = 0.06;
		_dYScoreUpperLengthThreshold = 0.15;

		//Just in case: eliminate old format file
		System::String ^ aPath      = System::Reflection::Assembly::GetExecutingAssembly()->Location;
		System::String ^ dirPath    = System::IO::Path::GetDirectoryName(aPath);
		System::String ^ fileName   = System::String::Concat(dirPath, "\\AutoSlicingConfig.xml");
	}
	//=============================================================================================================
	//Reading the parameters from the xml file
	void CAutoSliceOptions::ReadOptions()
	{
		System::String ^aPath = System::Reflection::Assembly::GetExecutingAssembly()->Location;
		System::String ^dirPath = System::IO::Path::GetDirectoryName(aPath);
		System::String ^fileName = System::String::Concat(dirPath, "\\AutoSlicingConfig.xml");

		try
		{
			StreamReader ^fs = gcnew StreamReader(fileName);
			XmlTextReader ^reader = gcnew XmlTextReader(fs);

			try
			{
				reader->WhitespaceHandling = WhitespaceHandling::None;
				reader->MoveToContent();
        
				//read group title tag
				reader->Read(); 

				ReadDoubleElement(reader, "TroughtoPeakRatioLowerLengthLimit", (_dTroughtoPeakRatioLowerLengthLimit));
				ReadDoubleElement(reader, "TroughtoPeakRatioUpperLengthLimit", _dTroughtoPeakRatioUpperLengthLimit);
				ReadDoubleElement(reader, "MinimumNumberofMoleculesLowerLength", _dMinimumNumberofMoleculesLowerLength);
				ReadDoubleElement(reader, "MinimumNumberofMoleculesUpperLength", _dMinimumNumberofMoleculesUpperLength);
				ReadDoubleElement(reader, "BoundaryLength", _dBoundaryLength);
				ReadDoubleElement(reader, "MinimumLowerPercentMolecules", _dMinimumLowerPercentMolecules);
				ReadDoubleElement(reader, "MinimumUpperPercentMolecules", _dMinimumUpperPercentMolecules);
				ReadDoubleElement(reader, "YScoreLowerLengthThreshold", _dYScoreLowerLengthThreshold);
				ReadDoubleElement(reader, "YScoreUpperLengthThreshold", _dYScoreUpperLengthThreshold);

				//read closing group tag
				reader->Read();
			}
			catch( char * str )
			{
				//Console::WriteLine(L"Error in reading xml settings: {0}", gcnew String(str));
				reader->Close();
				fs->Close();
				throw str;
			}
			catch(...)
			{
				//Console::WriteLine(L"Error in reading xml settings");
				reader->Close();
				fs->Close();
				throw;
			}

			reader->Close();
			fs->Close();
		}
		catch(char *str)
		{
			String ^errMsg;
			errMsg = gcnew String(str);
			//Console::WriteLine(L"Error Message: {0}", errMsg);
			SetDefaultParameters();
		}
		catch(...)
		{
			//Console::WriteLine("Could not read the configuration options; continue with defaults.");
			SetDefaultParameters();
		}
	}
	//=============================================================================================
	//Read double elements from XML
	void CAutoSliceOptions::ReadDoubleElement(XmlTextReader ^ reader, String ^sParameterText, double %dValue)
	{
		reader->ReadStartElement(sParameterText);
		dValue = (XmlConvert::ToDouble(reader->Value));
		//Console::WriteLine(L"Value read: {0}", dValue);
		reader->Read();
		reader->ReadEndElement();
	}
	//=============================================================================================================
	//Write AutoSlicing options to XML
	void CAutoSliceOptions::WriteOptions()
	{
		try
		{
			System::String ^aPath = System::Reflection::Assembly::GetExecutingAssembly()->Location;
			System::String ^dirPath = System::IO::Path::GetDirectoryName(aPath);

			StreamWriter ^fs = gcnew StreamWriter(System::String::Concat(dirPath, "\\AutoSlicingConfig.xml"), false);
			XmlTextWriter ^writer = gcnew XmlTextWriter(fs);;

			//Use indentation for readability.
			writer->Formatting = Formatting::Indented;
			writer->Indentation = 4;
         
			//Write the root
			writer->WriteStartDocument();
			writer->WriteStartElement("AutoSlicingOptions");
        
			//Write parameters
			WriteElement(writer, "TroughtoPeakRatioLowerLengthLimit", (_dTroughtoPeakRatioLowerLengthLimit));
			WriteElement(writer, "TroughtoPeakRatioUpperLengthLimit", (_dTroughtoPeakRatioUpperLengthLimit));
			WriteElement(writer, "MinimumNumberofMoleculesLowerLength", (_dMinimumNumberofMoleculesLowerLength));
			WriteElement(writer, "MinimumNumberofMoleculesUpperLength", (_dMinimumNumberofMoleculesUpperLength));
			WriteElement(writer, "BoundaryLength", (_dBoundaryLength));
			WriteElement(writer, "MinimumLowerPercentMolecules", (_dMinimumLowerPercentMolecules));
			WriteElement(writer, "MinimumUpperPercentMolecules", (_dMinimumUpperPercentMolecules));
			WriteElement(writer, "YScoreLowerLengthThreshold", (_dYScoreLowerLengthThreshold));
			WriteElement(writer, "YScoreUpperLengthThreshold", (_dYScoreUpperLengthThreshold));

			//Close tag for the root element.
			writer->WriteEndElement(); //"AutoSlicingOptions"
			writer->WriteEndDocument();
              
			//Write the XML to file and close the writer.
			writer->Close();
			fs->Close();
		}
		catch (...)
		{
			//Console::WriteLine("Writing Options Error");
		}
	}
	//=============================================================================================================
	//Write double elements to XML
	//=============================================================================================================
	void CAutoSliceOptions::WriteElement(XmlTextWriter ^writer, String ^sParameterText, double dValue)
	{
		 writer->WriteStartElement(sParameterText);
		 writer->WriteString(System::Convert::ToString(dValue));
		 writer->WriteEndElement();
	}
}
