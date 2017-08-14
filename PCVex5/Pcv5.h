//============================================================================
// Name        : Pcv5.h
// Author      : Ronny Haensch
// Version     : 1.0
// Copyright   : -
// Description : header file for fifth PCV assignment
// Submitted by: Group_G PCV ws[16/17]
//============================================================================

#include <iostream>
#include <fstream>
#include <list>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

class Pcv5{

	public:
		// constructor
		Pcv5(void){};
		// destructor
		~Pcv5(void){};
		
		// processing routine
		void run(string, string);

	private:
		// given functions
		int readMatchingPoints(string fname, Mat& x1, Mat& x2);
		int readControlPoints(string fname, Mat& x1, Mat& x2, Mat& Xm);
		void savePointList(string fname, Mat& points);
		// functions to be implemented
		// --> edit ONLY these functions!		
		Mat getCondition2D(Mat& p);
		Mat getCondition3D(Mat& p);
		Mat applyH(Mat& geomObj, Mat& H, string type);
		Mat solve_dlt(Mat& A);
		Mat getDesignMatrix_fundamental(Mat& n1, Mat& n2);
		Mat getDesignMatrix_homography3D(Mat& fst, Mat& snd);
		Mat getFundamentalMatrix(Mat& p1, Mat& p2);
		void forceSingularity(Mat& F);
		void decondition_fundamental(Mat& T1, Mat& T2, Mat& F);
		void decondition_homography3D(Mat& T1, Mat& T2, Mat& H);
		Mat homography3D(Mat& Xp, Mat& Xm);
		void defineCameras(Mat& F, Mat& P1, Mat& P2);
		void getEpipols(Mat& F, Mat& e1, Mat& e2);
		Mat makeSkewMatrix(Mat& e);

		Mat linearTriangulation(Mat& P1, Mat& P2, Mat& x1, Mat& x2);
};
