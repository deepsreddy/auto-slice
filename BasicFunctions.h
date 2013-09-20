#pragma once
using namespace System;
using namespace System::Runtime::InteropServices;

namespace AutoSlicing {
	public ref class CBasicFunctions
	{
		public:
			CBasicFunctions(void);
			static char* SystemStringToCharArray(String ^managedString);
	};
} // end of namespace FPConverter
