#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"

// BGR
// TODO pokombinowac jeszcze z odcieniem 
const int RED_TRESH[3] = {80, 80, 90};
// TODO
const int YELLOW_TRESH[3] = {0,0,0};
// HORIZONTAL_FILTER
// -1	-1	-1
//  2	 2	 2
// -1	-1	-1
const std::vector<int> HORIZONTAL_FILTER = {-1,-1,-1,2,2,2,-1,-1,-1};

void treshold(cv::Mat& inImg, cv::Mat& outImg) {
	CV_Assert(inImg.depth() != sizeof(uchar));
	cv::Mat_<cv::Vec3b> _I = inImg;
	cv::Mat_<cv::Vec3b> Iout = outImg;

	for( int i = 0; i < (inImg.rows); ++i)
		for (int j = 0; j < (inImg.cols); ++j) {
			Iout(i, j)[0] = 0;
			Iout(i, j)[1] = 0;
			Iout(i, j)[2] = 0;
			if (_I(i, j)[2] > RED_TRESH[0] && _I(i, j)[0] < RED_TRESH[1] &&  _I(i, j)[1] < RED_TRESH[2]) {
				Iout(i, j)[2] = 255;
				Iout(i, j)[1] = 255;
				Iout(i, j)[0] = 255;
			}
		}
}

// nieoptymalne, zmienic na greyscale
// nazwa: apply_filter? dodatkowy argument?
void find_horizontal_lines(cv::Mat& inImg, cv::Mat& outImg, std::vector<int> filter) {
	CV_Assert(inImg.depth() != sizeof(uchar));
	cv::Mat_<cv::Vec3b> _I = inImg;
	cv::Mat_<cv::Vec3b> Iout = outImg;
	int mask_size = filter.size();

	int check = floor(mask_size / 2);
	for (int i = check; i < _I.rows - check; ++i)
		for (int j = check; j < _I.cols - check; ++j) {
			std::vector<int> vec;
			int k = i - check;
			int count = 0;
			for (k, count; count < mask_size; k++, count++) {
				int l = j - check;
				int count2 = 0;
				for (l, count2; count2 < mask_size; l++, count2++) {
					int pix_value;
					pix_value = int((_I(k, l)[0] + _I(k, l)[1] + _I(k, l)[2]) / 3);
					vec.push_back(pix_value);
				}
			}
			int new_pix_value = 0;
			for (int i = 0; i < 9; ++i) {
				new_pix_value += HORIZONTAL_FILTER[i] * vec[i];
			}

			Iout(i, j)[2] = new_pix_value;
			Iout(i, j)[1] = new_pix_value;
			Iout(i, j)[0] = new_pix_value;
		}
}

void rotate_by(cv::Mat& inImg, cv::Mat& outImg, int angle) {

}

int main(int argc, char* argv[]) {
	if (argc == 2 && (cv::imread(argv[1]).data != NULL)) {
		cv::Mat image = cv::imread(argv[1]);
		cv::Mat tmp = cv::imread(argv[1]);
		treshold(image, tmp);
		cv::Mat tmp2 = tmp.clone();
		find_horizontal_lines(tmp, tmp2, HORIZONTAL_FILTER);

		cv::imshow("Shell", image);
		cv::imshow("Shell_tresh", tmp);
		cv::imshow("Shell_horizontal", tmp2);
		cv::waitKey(-1);
		return 0;
	} else {
		std::cout << "Wrong file name!" << std::endl;
		std::cout << "Usage: shell_recognition <file_name>" << std::endl;

		return -1;
	}
}
