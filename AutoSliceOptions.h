#pragma once
#using <mscorlib.dll>
#using <System.Xml.dll>
#using <system.dll>
using namespace System;
using namespace System::Xml; 

namespace AutoSlicing {

	public ref class CAutoSliceOptions
	{
		public:
			CAutoSliceOptions(void);
			void ReadOptions();
			void SetDefaultParameters(void);
			void WriteOptions();
			void operator = (const CAutoSliceOptions ^sXMLoptions);

		private:
			void ReadDoubleElement(XmlTextReader ^ reader, String ^sParameterText, double %dValue);
			void WriteElement(XmlTextWriter ^writer, String ^sParameterText, double dValue);

		public:
			double _dTroughtoPeakRatioLowerLengthLimit;
			double _dTroughtoPeakRatioUpperLengthLimit;
			double _dMinimumNumberofMoleculesLowerLength;
			double _dMinimumNumberofMoleculesUpperLength;
			double _dBoundaryLength;
			double _dMinimumLowerPercentMolecules;
			double _dMinimumUpperPercentMolecules;
			double _dYScoreLowerLengthThreshold;
			double _dYScoreUpperLengthThreshold;
	};
}


