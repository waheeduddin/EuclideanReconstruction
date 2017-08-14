//============================================================================
// Name        : Pcv5.cpp
// Author      : Ronny Haensch
// Version     : 1.0
// Copyright   : -
// Description : Triangulation
// submitted by: Group_G PCV ws[16/17]
//============================================================================

#include "Pcv5.h"

// compute fundamental matrix
/*
fst	first set of points
snd	second set of points
return	the estimated fundamental matrix
*/
Mat Pcv5::getFundamentalMatrix(Mat& fst, Mat& snd){
	// TODO !!!
	Mat T_fst = getCondition2D(fst);//get the condition matrix
	Mat T_snd = getCondition2D(snd);//get the condition matrix
	Mat fstPoints = applyH(fst, T_fst, "point");//get the conditioned points
	Mat sndPoints = applyH(snd, T_snd, "point");//get the conditioned points
	Mat A = getDesignMatrix_fundamental(fstPoints, sndPoints);
	Mat intermediateF = solve_dlt(A);
	forceSingularity(intermediateF);
	decondition_fundamental(T_fst, T_snd, intermediateF);

	return intermediateF;
}

// solve homogeneous equation system by usage of SVD
/*
A			the design matrix
return		the estimated fundamental matrix
*/
Mat Pcv5::solve_dlt(Mat& A){
	// this one function will be used to solve fundamental design matrix as well as homography design matrix. The algorithm is of course generic to
	//solve homogenous system in either 9 or 16 variables
	if(A.cols == 9){ //check if its a fundamental matrix 
		cv::SVD svd(A, SVD::FULL_UV);

		Mat F = Mat::eye(3, 3, CV_32FC1);

		Mat tempW = Mat::zeros(9, 1, CV_32FC1);
		for (int i = 0; i < svd.w.rows; i++) {
			tempW.at<float>(i, 0) = svd.w.at<float>(i, 0);  // making sure that we have 9 singular values. 
		}

		float min = tempW.at<float>(0, 0);
		int minIndex = 0;   //find the minimum singular value
		for (int i = 0; i < tempW.rows; i++) {
			if (min > tempW.at<float>(i, 0)) {
				min = tempW.at<float>(i, 0);
				minIndex = i;
			}
		}

		//reshaping the vt row into a 3x3 F matrix
		//implementing change due to feedback "you use arbitrary flips of signs at many many places in your code" 
		F.at<float>(0, 0) =  (svd.vt.at<float>(minIndex, 0));// before this statement was: F.at<float>(0, 0) =  -1 * (svd.vt.at<float>(minIndex, 0));
		F.at<float>(0, 1) =  (svd.vt.at<float>(minIndex, 1));
		F.at<float>(0, 2) =  (svd.vt.at<float>(minIndex, 2));

		F.at<float>(1, 0) =  (svd.vt.at<float>(minIndex, 3));
		F.at<float>(1, 1) =  (svd.vt.at<float>(minIndex, 4));
		F.at<float>(1, 2) =  (svd.vt.at<float>(minIndex, 5));

		F.at<float>(2, 0) = (svd.vt.at<float>(minIndex, 6));
		F.at<float>(2, 1) = (svd.vt.at<float>(minIndex, 7));
		F.at<float>(2, 2) = (svd.vt.at<float>(minIndex, 8));

		return F;
	}
	else{
		if (A.cols == 16) {//if its the design matrix of homography 3D
			cv::SVD svd(A, SVD::FULL_UV);

			Mat H = Mat::eye(4, 4, CV_32FC1);

			Mat tempW = Mat::zeros(16, 1, CV_32FC1);
			for (int i = 0; i < svd.w.rows; i++) {
				tempW.at<float>(i, 0) = svd.w.at<float>(i, 0);  // making sure that we have 16 singular values. 
			}

			float min = tempW.at<float>(0, 0);
			int minIndex = 0;   //find the minimum singular value
			for (int i = 0; i < tempW.rows; i++) {
				if (min >= tempW.at<float>(i, 0)) {
					min = tempW.at<float>(i, 0);
					minIndex = i;
				}
			}

			//reshaping the vt row into a 4x4 H matrix
			//implementing change due to feedback "you use arbitrary flips of signs at many many places in your code" 
			H.at<float>(0, 0) =  (svd.vt.at<float>(minIndex, 0));//before this statement was:H.at<float>(0, 0) =  -1 * (svd.vt.at<float>(minIndex, 0));
			H.at<float>(0, 1) =  (svd.vt.at<float>(minIndex, 1));
			H.at<float>(0, 2) =  (svd.vt.at<float>(minIndex, 2));
			H.at<float>(0, 3) =  (svd.vt.at<float>(minIndex, 3));

			H.at<float>(1, 0) =  (svd.vt.at<float>(minIndex, 4));
			H.at<float>(1, 1) =  (svd.vt.at<float>(minIndex, 5));
			H.at<float>(1, 2) =  (svd.vt.at<float>(minIndex, 6));
			H.at<float>(1, 3) =  (svd.vt.at<float>(minIndex, 7));

			H.at<float>(2, 0) =  (svd.vt.at<float>(minIndex, 8));
			H.at<float>(2, 1) =  (svd.vt.at<float>(minIndex, 9));
			H.at<float>(2, 2) =  (svd.vt.at<float>(minIndex, 10));
			H.at<float>(2, 3) =  (svd.vt.at<float>(minIndex, 11));

			H.at<float>(3, 0) =  (svd.vt.at<float>(minIndex, 12));
			H.at<float>(3, 1) =  (svd.vt.at<float>(minIndex, 13));
			H.at<float>(3, 2) =  (svd.vt.at<float>(minIndex, 14));
			H.at<float>(3, 3) =  (svd.vt.at<float>(minIndex, 15));

			return H;
		}
		else { // as a default case, a 3x3 NULL matrix is returned
			Mat F = Mat::zeros(3, 3, CV_32FC1);
			return F;
		}
	}	
}

// decondition a fundamental matrix that was estimated from conditioned point clouds
/*
T_fst	conditioning matrix of first set of points
T_snd	conditioning matrix of second set of points
F	conditioned fundamental matrix that has to be un-conditioned (in-place)
*/
void Pcv5::decondition_fundamental(Mat& T_fst, Mat& T_snd, Mat& F){
	// TODO !!!
	Mat deconditionedF = Mat::eye(3, 3, CV_32FC1);
	deconditionedF = T_snd.t() * F * T_fst;//learned in PCV lectures
	F = deconditionedF;
}

// define the design matrix as needed to compute fundamental matrix
/*
fst		first set of points
snd		second set of points
return		the design matrix to be computed
*/
Mat Pcv5::getDesignMatrix_fundamental(Mat& fst, Mat& snd){
	// TODO !!!
	int totalNumOfPoints = snd.cols;
	Mat A(totalNumOfPoints, 9, CV_32FC1);
	for (int i = 0; i < totalNumOfPoints; i++) {
		A.at<float>(i, 0) = fst.at<float>(0, i) * snd.at<float>(0, i);
		A.at<float>(i, 1) = fst.at<float>(1, i) * snd.at<float>(0, i);
		A.at<float>(i, 2) = snd.at<float>(0, i);

		A.at<float>(i, 3) = fst.at<float>(0, i) * snd.at<float>(1, i);
		A.at<float>(i, 4) = fst.at<float>(1, i) * snd.at<float>(1, i);
		A.at<float>(i, 5) = snd.at<float>(1, i);

		A.at<float>(i, 6) = fst.at<float>(0, i);
		A.at<float>(i, 7) = fst.at<float>(1, i);

		A.at<float>(i, 8) = 1.0;
	}

	return A;
}

// define the design matrix as needed to compute fundamental matrix
/*
fst	first set of points
snd	second set of points
A	the design matrix to be computed
*/
Mat Pcv5::getDesignMatrix_homography3D(Mat& fst, Mat& snd){
	// TODO !!!
	int size = snd.cols;
	Mat A(size * 3, 16, CV_32FC1);

	//change after feedback "wrong order of point sets in getDesignMatrix_homography3D"
	//The order of point sets has been changed. it is now the same as dictated by slide 14/20 of PCV lecture "pcv7_WS1617_3Dhomo"
	for (int i = 0; i < size; i++) {
		
		A.at<float>(i * 6, 0) = -snd.at<float>(3, i) * fst.at<float>(0, i);
		A.at<float>(i * 6, 1) = -snd.at<float>(3, i) * fst.at<float>(1, i);
		A.at<float>(i * 6, 2) = -snd.at<float>(3, i) * fst.at<float>(2, i);
		A.at<float>(i * 6, 3) = -snd.at<float>(3, i) * fst.at<float>(3, i);
		
		A.at<float>(i * 6, 4) = 0;
		A.at<float>(i * 6, 5) = 0;
		A.at<float>(i * 6, 6) = 0;
		A.at<float>(i * 6, 7) = 0;

		A.at<float>(i * 6, 8) = 0;
		A.at<float>(i * 6, 9) = 0;
		A.at<float>(i * 6, 10) = 0;
		A.at<float>(i * 6, 11) = 0;

		A.at<float>(i * 6, 12) = snd.at<float>(0, i) * fst.at<float>(0, i);
		A.at<float>(i * 6, 13) = snd.at<float>(0, i) * fst.at<float>(1, i);
		A.at<float>(i * 6, 14) = snd.at<float>(0, i) * fst.at<float>(2, i);
		A.at<float>(i * 6, 15) = snd.at<float>(0, i) * fst.at<float>(3, i);

		A.at<float>(i * 6 + 1, 0) = 0;
		A.at<float>(i * 6 + 1, 1) = 0;
		A.at<float>(i * 6 + 1, 2) = 0;
		A.at<float>(i * 6 + 1, 3) = 0;

		A.at<float>(i * 6 + 1, 4) = -snd.at<float>(3, i) * fst.at<float>(0, i);
		A.at<float>(i * 6 + 1, 5) = -snd.at<float>(3, i) * fst.at<float>(1, i);
		A.at<float>(i * 6 + 1, 6) = -snd.at<float>(3, i) * fst.at<float>(2, i);
		A.at<float>(i * 6 + 1, 7) = -snd.at<float>(3, i) * fst.at<float>(3, i);
		
		A.at<float>(i * 6 + 1, 8) = 0;
		A.at<float>(i * 6 + 1, 9) = 0;
		A.at<float>(i * 6 + 1, 10) = 0;
		A.at<float>(i * 6 + 1, 11) = 0;

		A.at<float>(i * 6 + 1, 12) = snd.at<float>(1, i) * fst.at<float>(0, i);
		A.at<float>(i * 6 + 1, 13) = snd.at<float>(1, i) * fst.at<float>(1, i);
		A.at<float>(i * 6 + 1, 14) = snd.at<float>(1, i) * fst.at<float>(2, i);
		A.at<float>(i * 6 + 1, 15) = snd.at<float>(1, i) * fst.at<float>(3, i);

		A.at<float>(i * 6 + 2, 0) = 0;
		A.at<float>(i * 6 + 2, 1) = 0;
		A.at<float>(i * 6 + 2, 2) = 0;
		A.at<float>(i * 6 + 2, 3) = 0;

		A.at<float>(i * 6 + 2, 4) = 0;
		A.at<float>(i * 6 + 2, 5) = 0;
		A.at<float>(i * 6 + 2, 6) = 0;
		A.at<float>(i * 6 + 2, 7) = 0;

		A.at<float>(i * 6 + 2, 8) = -snd.at<float>(3, i) * fst.at<float>(0, i);
		A.at<float>(i * 6 + 2, 9) = -snd.at<float>(3, i) * fst.at<float>(1, i);
		A.at<float>(i * 6 + 2, 10) = -snd.at<float>(3, i) * fst.at<float>(2, i);
		A.at<float>(i * 6 + 2, 11) = -snd.at<float>(3, i) * fst.at<float>(3, i);

		A.at<float>(i * 6 + 2, 12) = snd.at<float>(2, i) * fst.at<float>(0, i);
		A.at<float>(i * 6 + 2, 13) = snd.at<float>(2, i) * fst.at<float>(1, i);
		A.at<float>(i * 6 + 2, 14) = snd.at<float>(2, i) * fst.at<float>(2, i);
		A.at<float>(i * 6 + 2, 15) = snd.at<float>(2, i) * fst.at<float>(3, i);

		A.at<float>(i * 6 + 3, 0) = -snd.at<float>(3, i) * fst.at<float>(0, i);
		A.at<float>(i * 6 + 3, 1) = -snd.at<float>(3, i) * fst.at<float>(1, i);
		A.at<float>(i * 6 + 3, 2) = -snd.at<float>(3, i) * fst.at<float>(2, i);
		A.at<float>(i * 6 + 3, 3) = -snd.at<float>(3, i) * fst.at<float>(3, i);

		A.at<float>(i * 6 + 3, 4) = 0;
		A.at<float>(i * 6 + 3, 5) = 0;
		A.at<float>(i * 6 + 3, 6) = 0;
		A.at<float>(i * 6 + 3, 7) = 0;

		A.at<float>(i * 6 + 3, 8) = 0;
		A.at<float>(i * 6 + 3, 9) = 0;
		A.at<float>(i * 6 + 3, 10) = 0;
		A.at<float>(i * 6 + 3, 11) = 0;

		A.at<float>(i * 6 + 3, 12) = snd.at<float>(0, i) * fst.at<float>(0, i);
		A.at<float>(i * 6 + 3, 13) = snd.at<float>(0, i) * fst.at<float>(1, i);
		A.at<float>(i * 6 + 3, 14) = snd.at<float>(0, i) * fst.at<float>(2, i);
		A.at<float>(i * 6 + 3, 15) = snd.at<float>(0, i) * fst.at<float>(3, i);

		A.at<float>(i * 6 + 4, 0) = 0;
		A.at<float>(i * 6 + 4, 1) = 0;
		A.at<float>(i * 6 + 4, 2) = 0;
		A.at<float>(i * 6 + 4, 3) = 0;

		A.at<float>(i * 6 + 4, 4) = -snd.at<float>(3, i) * fst.at<float>(0, i);
		A.at<float>(i * 6 + 4, 5) = -snd.at<float>(3, i) * fst.at<float>(1, i);
		A.at<float>(i * 6 + 4, 6) = -snd.at<float>(3, i) * fst.at<float>(2, i);
		A.at<float>(i * 6 + 4, 7) = -snd.at<float>(3, i) * fst.at<float>(3, i);

		A.at<float>(i * 6 + 4, 8) = 0;
		A.at<float>(i * 6 + 4, 9) = 0;
		A.at<float>(i * 6 + 4, 10) = 0;
		A.at<float>(i * 6 + 4, 11) = 0;

		A.at<float>(i * 6 + 4, 12) = snd.at<float>(1, i) * fst.at<float>(0, i);
		A.at<float>(i * 6 + 4, 13) = snd.at<float>(1, i) * fst.at<float>(1, i);
		A.at<float>(i * 6 + 4, 14) = snd.at<float>(1, i) * fst.at<float>(2, i);
		A.at<float>(i * 6 + 4, 15) = snd.at<float>(1, i) * fst.at<float>(3, i);

		A.at<float>(i * 6 + 5, 0) = 0;
		A.at<float>(i * 6 + 5, 1) = 0;
		A.at<float>(i * 6 + 5, 2) = 0;
		A.at<float>(i * 6 + 5, 3) = 0;

		A.at<float>(i * 6 + 5, 4) = 0;
		A.at<float>(i * 6 + 5, 5) = 0;
		A.at<float>(i * 6 + 5, 6) = 0;
		A.at<float>(i * 6 + 5, 7) = 0;

		A.at<float>(i * 6 + 5, 8) = -snd.at<float>(3, i) * fst.at<float>(0, i);
		A.at<float>(i * 6 + 5, 9) = -snd.at<float>(3, i) * fst.at<float>(1, i);
		A.at<float>(i * 6 + 5, 10) = -snd.at<float>(3, i) * fst.at<float>(2, i);
		A.at<float>(i * 6 + 5, 11) = -snd.at<float>(3, i) * fst.at<float>(3, i);

		A.at<float>(i * 6 + 5, 12) = snd.at<float>(2, i) * fst.at<float>(0, i);
		A.at<float>(i * 6 + 5, 13) = snd.at<float>(2, i) * fst.at<float>(1, i);
		A.at<float>(i * 6 + 5, 14) = snd.at<float>(2, i) * fst.at<float>(2, i);
		A.at<float>(i * 6 + 5, 15) = snd.at<float>(2, i) * fst.at<float>(3, i);


		/*
		A.at<float>(i * 3, 0) = -fst.at<float>(3, i) * snd.at<float>(0, i);
		A.at<float>(i * 3, 1) = -fst.at<float>(3, i) * snd.at<float>(1, i);
		A.at<float>(i * 3, 2) = -fst.at<float>(3, i) * snd.at<float>(2, i);
		A.at<float>(i * 3, 3) = -fst.at<float>(3, i) * snd.at<float>(3, i);

		A.at<float>(i * 3, 4) = 0;
		A.at<float>(i * 3, 5) = 0;
		A.at<float>(i * 3, 6) = 0;
		A.at<float>(i * 3, 7) = 0;

		A.at<float>(i * 3, 8) = 0;
		A.at<float>(i * 3, 9) = 0;
		A.at<float>(i * 3, 10) = 0;
		A.at<float>(i * 3, 11) = 0;

		A.at<float>(i * 3, 12) = fst.at<float>(0, i) * snd.at<float>(0, i);
		A.at<float>(i * 3, 13) = fst.at<float>(0, i) * snd.at<float>(1, i);
		A.at<float>(i * 3, 14) = fst.at<float>(0, i) * snd.at<float>(2, i);
		A.at<float>(i * 3, 15) = fst.at<float>(0, i) * snd.at<float>(3, i);

		A.at<float>(i * 3 + 1, 0) = 0;
		A.at<float>(i * 3 + 1, 1) = 0;
		A.at<float>(i * 3 + 1, 2) = 0;
		A.at<float>(i * 3 + 1, 3) = 0;

		A.at<float>(i * 3 + 1, 4) = -fst.at<float>(3, i) * snd.at<float>(0, i);
		A.at<float>(i * 3 + 1, 5) = -fst.at<float>(3, i) * snd.at<float>(1, i);
		A.at<float>(i * 3 + 1, 6) = -fst.at<float>(3, i) * snd.at<float>(2, i);
		A.at<float>(i * 3 + 1, 7) = -fst.at<float>(3, i) * snd.at<float>(3, i);

		A.at<float>(i * 3 + 1, 8) = 0;
		A.at<float>(i * 3 + 1, 9) = 0;
		A.at<float>(i * 3 + 1, 10) = 0;
		A.at<float>(i * 3 + 1, 11) = 0;

		A.at<float>(i * 3 + 1, 12) = fst.at<float>(1, i) * snd.at<float>(0, i);
		A.at<float>(i * 3 + 1, 13) = fst.at<float>(1, i) * snd.at<float>(1, i);
		A.at<float>(i * 3 + 1, 14) = fst.at<float>(1, i) * snd.at<float>(2, i);
		A.at<float>(i * 3 + 1, 15) = fst.at<float>(1, i) * snd.at<float>(3, i);

		A.at<float>(i * 3 + 2, 0) = 0;
		A.at<float>(i * 3 + 2, 1) = 0;
		A.at<float>(i * 3 + 2, 2) = 0;
		A.at<float>(i * 3 + 2, 3) = 0;

		A.at<float>(i * 3 + 2, 4) = 0;
		A.at<float>(i * 3 + 2, 5) = 0;
		A.at<float>(i * 3 + 2, 6) = 0;
		A.at<float>(i * 3 + 2, 7) = 0;

		A.at<float>(i * 3 + 2, 8) = -fst.at<float>(3, i) * snd.at<float>(0, i);
		A.at<float>(i * 3 + 2, 9) = -fst.at<float>(3, i) * snd.at<float>(1, i);
		A.at<float>(i * 3 + 2, 10) = -fst.at<float>(3, i) * snd.at<float>(2, i);
		A.at<float>(i * 3 + 2, 11) = -fst.at<float>(3, i) * snd.at<float>(3, i);

		A.at<float>(i * 3 + 2, 12) = fst.at<float>(2, i) * snd.at<float>(0, i);
		A.at<float>(i * 3 + 2, 13) = fst.at<float>(2, i) * snd.at<float>(1, i);
		A.at<float>(i * 3 + 2, 14) = fst.at<float>(2, i) * snd.at<float>(2, i);
		A.at<float>(i * 3 + 2, 15) = fst.at<float>(2, i) * snd.at<float>(3, i);
		*/
	}

	return A;
}

// apply transformation
/*
geomObj		matrix with input objects (one per column)
H			matrix representing the transformation
type		the type of the geometric object (for now: only point and line)
return		transformed objects (one per column)
*/
Mat Pcv5::applyH(Mat& geomObj, Mat& H, string type){
	// TODO !!!
	if (type.compare("point") == 0) {
		// TO DO !!!
		Mat newPoint = H * geomObj;// APPARENTLY THE RESULTS OF THIS LINE ARE NOT THA SMAE AS IN MATLAB
		return newPoint;
	}
	// if object is a line
	if (type.compare("line") == 0) {
		// TO DO !!!
		Mat newLine = H.inv().t() * geomObj;
		return newLine;
	}
	cout << "ERROR: Do not know how to move " << type << endl;
	return geomObj.clone();
}


// get the conditioning matrix of given points
/*
p		the points as matrix
return		the condition matrix (already allocated)
*/
Mat Pcv5::getCondition2D(Mat& p){
	// TODO !!!
	double tx, ty, sx, sy;//have to comment sx sy for std implementation

	tx = 0.0;
	for (int i = 0; i < p.cols; i++) {
		tx = tx + p.at<float>(0, i);
	}
	tx = tx / p.cols;

	ty = 0.0;
	for (int i = 0; i < p.cols; i++) {
		ty = ty + p.at<float>(1, i);
	}
	ty = ty / p.cols;

	//change after feedback "getCondition2D overwrites input points!"
	Mat newP = Mat::zeros(2, p.cols, CV_32FC1);//before this staement was: Mat newP = p;


	for (int i = 0; i < newP.cols; i++) {
		newP.at<float>(0, i) = p.at<float>(0, i) - tx;//before this staement was:newP.at<float>(0, i) = newP.at<float>(0, i) - tx;
		newP.at<float>(1, i) = p.at<float>(1, i) - ty;
	}
	
	sx = 0.0;
	for (int i = 0; i < newP.cols; i++) {
	sx = sx + abs(newP.at<float>(0, i));
	}
	sx = sx / p.cols;

	sy = 0.0;
	for (int i = 0; i < newP.cols; i++) {
	sy = sy + abs(newP.at<float>(1, i));
	}
	sy = sy / newP.cols;
	
	Mat T(3, 3, CV_32FC1);

	T.at<float>(0, 0) = 1/sx;
	T.at<float>(1, 0) = 0;
	T.at<float>(2, 0) = 0;

	T.at<float>(0, 1) = 0;
	T.at<float>(1, 1) = 1/sy;
	T.at<float>(2, 1) = 0;
	T.at<float>(0, 2) = -tx * sx;
	T.at<float>(1, 2) = -ty * sy;
	T.at<float>(2, 2) = 1;

	return T;
}

// get the conditioning matrix of given 3D points
/*
p	the points as matrix
T	the condition matrix (already allocated)
*/
Mat Pcv5::getCondition3D(Mat& p){
	double tx, ty, tz, sx, sy, sz;

	tx = 0.0;
	for (int i = 0; i < p.cols; i++) {
		tx = tx + p.at<float>(0, i);
	}
	tx = tx / p.cols;

	ty = 0.0;
	for (int i = 0; i < p.cols; i++) {
		ty = ty + p.at<float>(1, i);
	}
	ty = ty / p.cols;

	tz = 0.0;
	for (int i = 0; i < p.cols; i++) {
		tz = tz + p.at<float>(2, i);
	}
	tz = tz / p.cols;

	//change after feedback "getCondition2D overwrites input points!"
	Mat newP = Mat::zeros(3, p.cols, CV_32FC1);//before this staement was: Mat newP = p;

	for (int i = 0; i < newP.cols; i++) {
		newP.at<float>(0, i) = p.at<float>(0, i) - tx;//before this staement was:newP.at<float>(0, i) = newP.at<float>(0, i) - tx; 
		newP.at<float>(1, i) = p.at<float>(1, i) - ty;
		newP.at<float>(2, i) = p.at<float>(2, i) - tz;
	}

	sx = 0.0;
	for (int i = 0; i < newP.cols; i++) {
		sx = sx + abs(newP.at<float>(0, i));
	}
	sx = sx / newP.cols;

	sy = 0.0;
	for (int i = 0; i < newP.cols; i++) {
		sy = sy + abs(newP.at<float>(1, i));
	}
	sy = sy / newP.cols;

	sz = 0.0;
	for (int i = 0; i < newP.cols; i++) {
		sz = sz + abs(newP.at<float>(2, i));
	}
	sz = sz / newP.cols;

	Mat T(4, 4, CV_32FC1);
	T.at<float>(0, 0) = 1 / sx;
	T.at<float>(1, 0) = 0;
	T.at<float>(2, 0) = 0;
	T.at<float>(3, 0) = 0;

	T.at<float>(0, 1) = 0;
	T.at<float>(1, 1) = 1 / sy;
	T.at<float>(2, 1) = 0;
	T.at<float>(3, 1) = 0;

	T.at<float>(0, 2) = 0;
	T.at<float>(1, 2) = 0;
	T.at<float>(2, 2) = 1 / sz;
	T.at<float>(3, 2) = 0;

	T.at<float>(0, 3) = -tx / sx;
	T.at<float>(1, 3) = -ty / sy;
	T.at<float>(2, 3) = -tz / sz;
	T.at<float>(3, 3) = 1;

	return T;
	/*double tx, ty , tz;

	tx = 0.0;
	for (int i = 0; i < p.cols; i++) {
		tx = tx + p.at<float>(0, i);
	}
	tx = tx / p.cols;

	ty = 0.0;
	for (int i = 0; i < p.cols; i++) {
		ty = ty + p.at<float>(1, i);
	}
	ty = ty / p.cols;

	tz = 0.0;
	for (int i = 0; i < p.cols; i++) {
		tz = tz + p.at<float>(2, i);
	}
	tz = tz / p.cols;

	Mat newP = Mat::zeros(3, p.cols, CV_32FC1);//before this staement was: Mat newP = p;

	for (int i = 0; i < newP.cols; i++) {
		newP.at<float>(0, i) = p.at<float>(0, i) - tx;//before this staement was:newP.at<float>(0, i) = newP.at<float>(0, i) - tx; 
		newP.at<float>(1, i) = p.at<float>(1, i) - ty;
		newP.at<float>(2, i) = p.at<float>(2, i) - tz;
	}
	//The square root of mean distance implementation of the conditioning matrix
	Mat dist = Mat::zeros(newP.cols, 1, CV_32FC1);
	for (int i = 0; i < newP.cols; i++) {
		dist.at<float>(i, 0) = sqrt((newP.at<float>(0, i) *(newP.at<float>(0, i))) + (newP.at<float>(1, i) *newP.at<float>(1, i)) + (newP.at<float>(2, i) *newP.at<float>(2, i)));
	}
	double meanDist = 0.0;
	for (int i = 0; i < newP.cols; i++) {
		meanDist = meanDist + dist.at<float>(i, 0);
	}
	meanDist = meanDist / newP.cols;

	double scale = sqrt(3) / meanDist;

	Mat T(4, 4, CV_32FC1);
	T.at<float>(0, 0) = scale;
	T.at<float>(1, 0) = 0;
	T.at<float>(2, 0) = 0;
	T.at<float>(3, 0) = 0;

	T.at<float>(0, 1) = 0;
	T.at<float>(1, 1) = scale;
	T.at<float>(2, 1) = 0;
	T.at<float>(3, 1) = 0;

	T.at<float>(0, 2) = 0;
	T.at<float>(1, 2) = 0;
	T.at<float>(2, 2) = scale;
	T.at<float>(3, 2) = 0;

	T.at<float>(0, 3) = -tx * scale;
	T.at<float>(1, 3) = -ty * scale;
	T.at<float>(2, 3) = -tz * scale;
	T.at<float>(3, 3) = 1;

	return T;*/
}

// enforce rank of 2 on fundamental matrix
/*
F	the matrix to be changed
*/
void Pcv5::forceSingularity(Mat& F) {
	// TODO !!!
	cv::SVD svd(F, SVD::FULL_UV);

	Mat newF = Mat::eye(3, 3, CV_32FC1);

	Mat diagonals = Mat::zeros(3, 3, CV_32FC1);//here the 3rd element of diagonal is already zero
	diagonals.at<float>(0, 0) = svd.w.at<float>(0, 0);//singular value of F as the diagonal element
	diagonals.at<float>(1, 1) = svd.w.at<float>(1, 0);

	newF = svd.u * diagonals * svd.vt; // re calculate F from its SVD components

	F = newF; //set the cahnge to the original F
}

// calculates epipols from fundamental matrix
/*
F	the fundamental matrix
e1	first epipol
e2	second epipol
*/
void Pcv5::getEpipols(Mat& F, Mat& e1, Mat& e2) {
	// TODO !!!
	Mat firstEpiPole = Mat::zeros(3, 1, CV_32FC1);
	Mat secondEpiPole = Mat::zeros(3, 1, CV_32FC1);

	cv::SVD svd(F, SVD::FULL_UV);
	Mat temp = svd.vt.t();
	for (int i = 0; i < 3; i++) {
		firstEpiPole.at<float>(i, 0) = temp.at<float>(i, 2); // the right null space of F
		secondEpiPole.at<float>(i, 0) = svd.u.at<float>(i, 2);//the left null space of F
	}
	e1 = firstEpiPole;
	e2 = secondEpiPole;
}

// generates skew matrix from vector
/*
e		given vector
return		skew matrix
*/
Mat Pcv5::makeSkewMatrix(Mat& e) {
	// TODO !!!
	Mat skewMatrix = Mat::zeros(3, 3, CV_32FC1);

	if (e.rows == 3) {
		skewMatrix.at<float>(1, 0) = e.at<float>(2, 0);
		skewMatrix.at<float>(0, 1) = -1 * e.at<float>(2, 0);

		skewMatrix.at<float>(0, 2) = e.at<float>(1, 0);
		skewMatrix.at<float>(2, 0) = -1 * e.at<float>(1, 0);

		skewMatrix.at<float>(2, 1) = e.at<float>(0, 0);
		skewMatrix.at<float>(1, 2) = -1 * e.at<float>(0, 0);
	}
	else {
		if (e.cols == 3) {
			skewMatrix.at<float>(1, 0) = e.at<float>(0, 2);
			skewMatrix.at<float>(0, 1) = -1 * e.at<float>(0, 2);

			skewMatrix.at<float>(0, 2) = e.at<float>(0, 1);
			skewMatrix.at<float>(2, 0) = -1 * e.at<float>(0, 1);

			skewMatrix.at<float>(2, 1) = e.at<float>(0, 0);
			skewMatrix.at<float>(1, 2) = -1 * e.at<float>(0, 0);
		}
		else {
			cout << "ERROR: skew matrix could not be made because the initializer is not size 3." << endl << "Returning zeros!" << endl;
		}
	}
	return skewMatrix;
}

// generates 2 projection matrices by using fundamental matrix
/*
F	the fundamental matrix
P1	first projection matrix (standard P matrix)
P2	second projection matrix based on F
*/
void Pcv5::defineCameras(Mat& F, Mat& P1, Mat& P2) {
	// TODO !!!
	Mat standardP = Mat::zeros(3, 4, CV_32FC1);

	standardP.at<float>(0, 0) = 1;
	standardP.at<float>(1, 1) = 1;
	standardP.at<float>(2, 2) = 1;

	P1 = standardP;

	Mat secondP = Mat::zeros(3, 4, CV_32FC1);
	Mat e1 = Mat::zeros(3, 1, CV_32FC1);
	Mat e2 = Mat::zeros(3, 1, CV_32FC1);

	getEpipols(F, e1, e2);
	Mat skewE2 = Mat::zeros(3, 3, CV_32FC1);
	skewE2 = makeSkewMatrix(e2);

	Mat temp = skewE2 * F;
	for (int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++){
			secondP.at<float>(i, j) = temp.at<float>(i, j);
		}
	}
	secondP.at<float>(0, 3) = e2.at<float>(0,0);
	secondP.at<float>(1, 3) = e2.at<float>(1, 0);
	secondP.at<float>(2, 3) = e2.at<float>(2, 0);

	P2 = secondP;
}

// triangulates given set of image points based on projection matrices
/*
P1	projection matrix of first image
P2	projection matrix of second image
x1	image point set of first image
x2	image point set of second image
return	triangulated object points
*/
Mat Pcv5::linearTriangulation(Mat& P1, Mat& P2, Mat& x1, Mat& x2){
	// TODO !!!
	int numOfPoints = x1.cols;
	Mat X0 = Mat::zeros(4, numOfPoints, CV_32FC1);

	Mat designMatrix = Mat::zeros(4, 4,CV_32FC1);
	
	
	Mat P1_1 = P1.row(0);
	Mat P1_2 = P1.row(1);
	Mat P1_3 = P1.row(2);

	Mat P2_1 = P2.row(0);
	Mat P2_2 = P2.row(1);
	Mat P2_3 = P2.row(2);
	for (int i = 0; i < x1.cols; i++) {
		designMatrix.at<float>(0, 0) = (x1.at<float>(0, i) * P1_3.at<float>(0, 0)) - P1_1.at<float>(0, 0);
		designMatrix.at<float>(0, 1) = (x1.at<float>(0, i) * P1_3.at<float>(0, 1)) - P1_1.at<float>(0, 1);
		designMatrix.at<float>(0, 2) = (x1.at<float>(0, i) * P1_3.at<float>(0, 2)) - P1_1.at<float>(0, 2);
		designMatrix.at<float>(0, 3) = (x1.at<float>(0, i) * P1_3.at<float>(0, 3)) - P1_1.at<float>(0, 3);

		designMatrix.at<float>(1, 0) = (x1.at<float>(1, i) * P1_3.at<float>(0, 0)) - P1_2.at<float>(0, 0);
		designMatrix.at<float>(1, 1) = (x1.at<float>(1, i) * P1_3.at<float>(0, 1)) - P1_2.at<float>(0, 1);
		designMatrix.at<float>(1, 2) = (x1.at<float>(1, i) * P1_3.at<float>(0, 2)) - P1_2.at<float>(0, 2);
		designMatrix.at<float>(1, 3) = (x1.at<float>(1, i) * P1_3.at<float>(0, 3)) - P1_2.at<float>(0, 3);

		designMatrix.at<float>(2, 0) = (x2.at<float>(0, i) * P2_3.at<float>(0, 0)) - P2_1.at<float>(0, 0);
		designMatrix.at<float>(2, 1) = (x2.at<float>(0, i) * P2_3.at<float>(0, 1)) - P2_1.at<float>(0, 1);
		designMatrix.at<float>(2, 2) = (x2.at<float>(0, i) * P2_3.at<float>(0, 2)) - P2_1.at<float>(0, 2);
		designMatrix.at<float>(2, 3) = (x2.at<float>(0, i) * P2_3.at<float>(0, 3)) - P2_1.at<float>(0, 3);

		designMatrix.at<float>(3, 0) = (x2.at<float>(1, i) * P2_3.at<float>(0, 0)) - P2_2.at<float>(0, 0);
		designMatrix.at<float>(3, 1) = (x2.at<float>(1, i) * P2_3.at<float>(0, 1)) - P2_2.at<float>(0, 1);
		designMatrix.at<float>(3, 2) = (x2.at<float>(1, i) * P2_3.at<float>(0, 2)) - P2_2.at<float>(0, 2);
		designMatrix.at<float>(3, 3) = (x2.at<float>(1, i) * P2_3.at<float>(0, 3)) - P2_2.at<float>(0, 3);

		cv::SVD svd(designMatrix, SVD::FULL_UV);

		Mat tempW = Mat::zeros(4, 1, CV_32FC1);
		for (int i = 0; i < svd.w.rows; i++) {
			tempW.at<float>(i, 0) = svd.w.at<float>(i, 0);  // making sure that we have 9 singular values. 
		}
		float min = tempW.at<float>(0, 0);
		int minIndex = 0;   //find the minimum singular value
		for (int i = 0; i < tempW.rows; i++) {
			if (min > tempW.at<float>(i, 0)) {
				min = tempW.at<float>(i, 0);
				minIndex = i;
			}
		}

		//reshaping the vt row into a 4X1 homogenous 3D point vector
		//implementing change due to feedback "you use arbitrary flips of signs at many many places in your code"
		X0.at<float>(0,i) = (svd.vt.at<float>(minIndex, 0)); // befoe this statement was: X0.at<float>(0,i) = -1 * (svd.vt.at<float>(minIndex, 0));
		X0.at<float>(1,i) = (svd.vt.at<float>(minIndex, 1));
		X0.at<float>(2,i) = (svd.vt.at<float>(minIndex, 2));
		X0.at<float>(3,i) = (svd.vt.at<float>(minIndex, 3));
	}
	
	return X0;
}

// computes 3D homography
/*
X1	first set of points
X2	second set of points
H	computed homography
*/
Mat Pcv5::homography3D(Mat& X1, Mat& X2){
	// TODO !!!
	
	Mat T3D_1 = getCondition3D(X1);
	Mat T3D_2 = getCondition3D(X2);

	Mat new_X1 = applyH(X1, T3D_1, "point");
	Mat new_X2 = applyH(X2, T3D_2, "point");

	Mat A_homography = getDesignMatrix_homography3D(new_X1, new_X2);
	Mat conditionedH = solve_dlt(A_homography);
	
	//implementing feedback change "homography not deconditioned"
	decondition_homography3D(T3D_1, T3D_2, conditionedH);//before this statement was missing.

	return conditionedH;
}

// decondition a homography that was estimated from conditioned point clouds
/*
T_to	conditioning matrix of first set of points
T_from	conditioning matrix of second set of points
H	conditioned homography that has to be un-conditioned (in-place)
*/
void Pcv5::decondition_homography3D(Mat& T_to, Mat& T_from, Mat& H){
	// TODO !!!
	Mat deconditionedH = Mat::eye(3, 3, CV_32FC1);
	deconditionedH = T_to.inv() * H * T_from;
	H = deconditionedH;
}

/* *****************************
  GIVEN FUNCTIONS
***************************** */

// function loads input image, calls processing function, and saves result
/*
fname	path to input image
*/
void Pcv5::run(string imgPoints, string ctrPoints){

   // get corresponding points within the two images
    Mat x1, x2;
    int numberOfPointPairs = readMatchingPoints(imgPoints, x1, x2);
    
    // just some putput
    cout << "Number of defined point pairs: " << numberOfPointPairs << endl;
    
    // calculate fundamental matrix
    Mat F = getFundamentalMatrix(x1, x2);

	F.convertTo(F, CV_32FC1);//making sure the last element is 1
	F = F / F.at<float>(2, 2);

    cout << endl << "Fundamental matrix: " << F << endl;

    // define projection matrices
    Mat P1, P2;
    defineCameras(F, P1, P2);
   
    cout << endl << "Camera 1: " << P1 << endl;
    cout << endl << "Camera 2: " << P2 << endl;
    
    // linear triangulation of image points
    // resulting projective reconstruction
    Mat X = linearTriangulation(P1, P2, x1, x2);

    // save reconstructed points to file
    savePointList("projectiveReconstruction.asc", X);
    
    // read control points
    // clean re-use of x1- and x2-matrices
    Mat Xm;
    numberOfPointPairs = readControlPoints(ctrPoints, x1, x2, Xm);

    cout << endl << "Control points:" << endl;
    cout << "x1: " << x1 << endl;
    cout << "x2: " << x2 << endl;
    cout << "Xm: " << Xm << endl;

    // linear triangulation of image points
    Mat Xp = linearTriangulation(P1, P2, x1, x2);
    cout << "Reconstructed control points: " << Xp << endl;

    // Transform projective reconstructed points to euclidian reconstruction
    Mat H = homography3D(Xp, Xm);
	H.convertTo(H, CV_32FC1);//making the last enitiy '1'. this is done because it was done in testing part of exercise 2 as well.
	H = H / H.at<float>(3, 3);

    Mat X_final = applyH(X, H, "point");
    
    // save reconstructed points to file
    savePointList("euclidianReconstruction_new.asc", X_final);

}

// saves point list to file
/*
fname		file name
points		matrix of points
*/
void Pcv5::savePointList(string fname, Mat& points){

  // open file for write
  fstream file(fname.c_str(), ios::out);
  if(!file){
      cerr << "ERROR: cannot open file " << fname << endl;
	  cerr << "Press enter to continue..." << endl;
	  cin.get();
      exit(-1);
  }
  
  // if homogeneous points: norm and write points
  if (points.rows == 4){
      for(int i=0; i<points.cols; i++){
		file << points.at<float>(0, i)/points.at<float>(3, i) << "," << points.at<float>(1, i)/points.at<float>(3, i) << "," << points.at<float>(2, i)/points.at<float>(3, i) << endl;
      }
  }
  // if euclidian points: write points
  if (points.rows == 3){
      for(int i=0; i<points.cols; i++){
		file << points.at<float>(0, i) << "," << points.at<float>(1, i) << "," << points.at<float>(2, i) << endl;
      }
  }
  
  // close file
  file.close();

}

// read homologeous points from file
/*
fname	path of file containing point list
x1	pointer to matrix containing points of first image
x2	pointer to matrix containing points of second image
*/
int Pcv5::readMatchingPoints(string fname, Mat& x1, Mat& x2){

    // open file
    fstream file(fname.c_str(), ios::in);
    if (!file){
		cerr << "ERROR: Cannot open file " << fname << endl;
	    cerr << "Press enter to continue..." << endl;
	    cin.get();
        exit(-1);
    }

    // read points into list
    list<Scalar> points;
    int end, numberOfPointPairs = 0;
    string buffer;
    // read line by line
    while(getline(file, buffer)){

		// counting
		numberOfPointPairs++;
		
		// get semicolon
		end = buffer.find(';');
		// first point before semicolon
		string fst = buffer.substr(0, end);
		// second point after semicolon
		string snd = buffer.substr(end+1);
		
		// get x and y
		Scalar cur;
		end = fst.find(',');
		cur.val[0] = atof(fst.substr(0, end).c_str());
		cur.val[1] = atof(fst.substr(end+1).c_str());
		
		// get x and y
		end = snd.find(',');
		cur.val[2] = atof(snd.substr(0, end).c_str());
		cur.val[3] = atof(snd.substr(end+1).c_str());
		
		// push point pair to list
		points.push_back(cur);
    }
    
    // allocate memory for point matrices
    x1 = Mat(3, numberOfPointPairs, CV_32FC1);
    x2 = Mat(3, numberOfPointPairs, CV_32FC1);
    
    // fill point matrices
    int i=0;
    for(list<Scalar>::iterator p = points.begin(); p != points.end(); p++,i++){
		// x1
		x1.at<float>(0, i) = (*p).val[0];
		x1.at<float>(1, i) = (*p).val[1];
		x1.at<float>(2, i) = 1;
		// x2
		x2.at<float>(0, i) = (*p).val[2];
		x2.at<float>(1, i) = (*p).val[3];
		x2.at<float>(2, i) = 1;
    }
    return numberOfPointPairs;
}

// read control points from file
/*
fname	path of file containing point list
x1	pointer to matrix containing points of first image
x2	pointer to matrix containing points of second image
Xm	pointer to matrix containing object points
*/
int Pcv5::readControlPoints(string fname, Mat& x1, Mat& x2, Mat& Xm){

    // open file
    fstream file(fname.c_str(), ios::in);
    if (!file){
		cerr << "ERROR: Cannot open file " << fname << endl;
	    cerr << "Press enter to continue..." << endl;
	    cin.get();
        exit(-1);
    }

    // read points into list
    list< vector<double> > points;
    int pos1, pos2, numberOfPointPairs = 0;
    string buffer;
    // read line by line
    while(getline(file, buffer)){
      
		// counting
		numberOfPointPairs++;
		
		// get semicolons
		pos1 = buffer.find(';');
		pos2 = buffer.rfind(';');
		// first point before first semicolon
		string fst = buffer.substr(0, pos1);
		// second point between semicolons
		string snd = buffer.substr(pos1+1, pos2);
		// third point after last semicolon
		string thd = buffer.substr(pos2+1);
		
		// current point triplet
		vector<double> cur;

		// get x and y of first point
		pos1 = fst.find(',');
		cur.push_back( atof(fst.substr(0, pos1).c_str()) );
		cur.push_back( atof(fst.substr(pos1+1).c_str()) );
		// get x and y of second point
		pos1 = snd.find(',');
		cur.push_back( atof(snd.substr(0, pos1).c_str()) );
		cur.push_back( atof(snd.substr(pos1+1).c_str()) );
		// get x, y, and z of object point
		pos1 = thd.find(',');
		pos2 = thd.rfind(',');
		cur.push_back( atof(thd.substr(0, pos1).c_str()) );
		cur.push_back( atof(thd.substr(pos1+1, pos2).c_str()) );
		cur.push_back( atof(thd.substr(pos2+1).c_str()) );
		
		// push point triplet to list
		points.push_back(cur);
    }
    
    // allocate memory for point matrices
    x1 = Mat(3, numberOfPointPairs, CV_32FC1);
    x2 = Mat(3, numberOfPointPairs, CV_32FC1);
    Xm = Mat(4, numberOfPointPairs, CV_32FC1);
    
    // fill point matrices
    int i=0;
    for(list< vector<double> >::iterator p = points.begin(); p != points.end(); p++,i++){
		// x1
		x1.at<float>(0, i) = (*p).at(0);
		x1.at<float>(1, i) = (*p).at(1);
		x1.at<float>(2, i) = 1;
		// x2
		x2.at<float>(0, i) = (*p).at(2);
		x2.at<float>(1, i) = (*p).at(3);
		x2.at<float>(2, i) = 1;
		// Xm
		Xm.at<float>(0, i) = (*p).at(4);
		Xm.at<float>(1, i) = (*p).at(5);
		Xm.at<float>(2, i) = (*p).at(6);
		Xm.at<float>(3, i) = 1;
    }
    
    return numberOfPointPairs;
}
