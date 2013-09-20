#include "stdafx.h"
#include "PreselectionItem.h"

namespace AutoSlicing 
{        

    //double LowerLenLimit { get; set; }
    //double UpperLenLimit { get; set; }
    //double CutoffDistanceValue { get; set; }
    //double CutoffDistanceMin { get; set; }

	PreselectionItem::PreselectionItem()
    {
        LowerLenLimit = 0.0;
        UpperLenLimit = 0.0;
        CutoffDistanceValue = 0.0;
        CutoffDistanceMin = -2.0;
    }


    PreselectionItem::PreselectionItem(double lowerLenLimit, double upperLenLimit, double cutoffDistanceValue, double cutoffDistanceMin, bool _real)
    {
        LowerLenLimit = lowerLenLimit;
        UpperLenLimit = upperLenLimit;
        CutoffDistanceValue = cutoffDistanceValue;
        CutoffDistanceMin = cutoffDistanceMin;
		_Real = _real;
    }

};