//============================================================================
// Name        : main.cpp
// Author      : Ronny Haensch
// Version     : 1.0
// Copyright   : -
// Description : only calls processing and test routines
//============================================================================

#include <iostream>

#include "Pcv5.h"

using namespace std;

// usage: path to image in argv[1]
// main function. loads and saves image
int main(int argc, char** argv) {

	// will contain path to point lists (taken from argv)
	string imgPoints, ctrPoints;

    // check if image paths were defined
    if (argc != 3){
		cerr << "Usage: pcv5 <path to image points> <path to control points>" << endl;
		cerr << "Press enter to continue..." << endl;
		cin.get();
		return -1;
    }else{
	    // if yes, assign it to variable fname
	    imgPoints = argv[1];
	    ctrPoints = argv[2];
	}
	
	// construct processing object
	Pcv5 pcv5;

	// start processing
	pcv5.run(imgPoints, ctrPoints);

	cout << "Press enter to continue..." << endl;
	cin.get();

	return 0;

}
