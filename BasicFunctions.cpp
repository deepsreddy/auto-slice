#include "stdafx.h"
#include "BasicFunctions.h"

namespace AutoSlicing {

	// default constructor
	CBasicFunctions::CBasicFunctions(void)
	{
	}

	//Convert System String to char array
	char * CBasicFunctions::SystemStringToCharArray(String ^managedString)
	{
		IntPtr ptrToNativeString = Marshal::StringToHGlobalAnsi(managedString);
		char* nativeString = static_cast<char*>(ptrToNativeString.ToPointer());
		return nativeString;
	}
} // end of namespace FPConverter
