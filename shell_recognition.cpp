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

const std::vector<int> HORIZONTAL_FILTER = {-1,-1,-1,
											 2, 2, 2,
											-1,-1,-1};

const std::vector<int> HORIZONTAL_FILTER2 = {2,2,2,2,2,2,2};

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

std::vector<std::vector<int>> find_horizontal_lines(cv::Mat& inImg, cv::Mat& outImg) {
	CV_Assert(inImg.depth() != sizeof(uchar));
	cv::Mat_<cv::Vec3b> _I = inImg;
	cv::Mat_<cv::Vec3b> Iout = outImg;

	std::vector<std::vector<int>> lines;
	for (int i = 0; i < _I.rows; ++i) {
		int line_begin = -1;
		int last_pixel_ok = 0;
		int line_length = 0;
		for (int j = 0; j < _I.cols; ++j) {
			Iout(i, j)[0] = 0;
			Iout(i, j)[1] = 0;
			Iout(i, j)[2] = 0;
			if (_I(i, j)[2] == 255 && _I(i, j)[1] == 255 && _I(i, j)[0] == 255 && line_begin == -1) {
				line_begin = j;
				line_length++;
				last_pixel_ok = 1;
			} else if (_I(i, j)[2] == 255 && _I(i, j)[1] == 255 && _I(i, j)[0] == 255 && line_begin != -1 && last_pixel_ok == 1) {
				line_length++;
			} else if (line_begin != -1) {
				last_pixel_ok = 0;
				// TODO ogarniac ten treshold automatycznie
				if (line_length > 20) {
					std::vector<int> line;
					line.push_back(i);
					line.push_back(line_begin);
					line.push_back(line_length);
					lines.push_back(line);
					for (int k = 0; k < line_length; ++k) {
						Iout(i, line_begin+k)[0] = 255;
						Iout(i, line_begin+k)[1] = 255;
						Iout(i, line_begin+k)[2] = 255;
					}
				}
				line_begin = -1;
				line_length = 0;
			}
		}
	}
	//for (int i = 0; i < lines.size(); ++i) {
	//	std::cout << "I: " << lines[i][0] << " Length of line: " << lines[i][2] << std::endl;
	//}
	return lines;
}

std::vector<int> get_logos_possible_bounds(std::vector<std::vector<int>> lines) {
	for (int i = 0; i < lines.size(); ++i) {
		std::cout << "I: " << lines[i][0] << " Begin: " << lines[i][1] << " Length of line: " << lines[i][2] << std::endl;
	}

	int start_flag = -1;
	int tmp_length = 0;
	int tmp_begin = 0;
	int tmp_row = 0;
	int tmp_thickness = 0;
	for (int i = 0; i < lines.size(); ++i) {
		if (start_flag == -1) {
			tmp_length = lines[i][2];
			tmp_begin = lines[i][1];
			tmp_row = lines[i][0];
			tmp_thickness++;
			start_flag = 1;
			continue;
		}
		if (lines[i][2] > tmp_length && lines[i][1] < tmp_begin && lines[i][1] > (tmp_begin - tmp_begin*6/100)) {
			tmp_length = lines[i][2];
			tmp_begin = lines[i][1];
			tmp_row = lines[i][0];
			tmp_thickness++;
		}
		// ogarnac przypadek kiedy na jednej linii sa dwa prawdopodobne loga
		// (czyli tu else w trakcie, albo kolejnego szukac w liniach ponizej)

	}
	std::cout << " " << std::endl;
	std::cout << "I: " << tmp_row << " Begin: " << tmp_begin << " Length of line: "
		<< tmp_length << " Thickness: " << tmp_thickness << std::endl;

	int probable_bottom_line = tmp_thickness * 11;
	int lower_row_bound = tmp_row - tmp_thickness + probable_bottom_line - probable_bottom_line*10/100;
	int upper_row_bound = tmp_row - tmp_thickness + probable_bottom_line + probable_bottom_line*10/100;
	std::cout << "Lower row bound: " << lower_row_bound << std::endl;
	std::cout << "Upper row bound: " << upper_row_bound << std::endl;
	int bottom_begin = 0;
	int bottom_row = 0;
	int bottom_length = 0;
	for (int i = 0; i < lines.size(); ++i) {
		if (lines[i][0] > (lower_row_bound) && lines[i][0] < (upper_row_bound) &&
			   	lines[i][1] > (tmp_begin - tmp_begin*10/100) &&
				lines[i][1] < (tmp_begin + tmp_begin*10/100) &&
			   	lines[i][2] > tmp_length) {
			bottom_begin = lines[i][1];
			bottom_row = lines[i][0];
			bottom_length = lines[i][2];

			std::cout << "I: " << lines[i][0] << " Begin: " << lines[i][1] << " Length of line: " << lines[i][2] << std::endl;
			break;
		}
	}

	// Calculate borders
	int left_x = bottom_begin - bottom_begin*30/100;
	int right_x = bottom_begin + bottom_length + bottom_begin*30/100;
	int upper_y = tmp_row - tmp_row*15/100;
	int lower_y = bottom_row + bottom_row*6/100;
	std::vector<int> bounds = {left_x, right_x, upper_y, lower_y};
	return bounds;
}
// nieoptymalne, zmienic na greyscale
// nazwa: apply_filter? dodatkowy argument?
void apply_filter(cv::Mat& inImg, cv::Mat& outImg, std::vector<int> filter) {
	CV_Assert(inImg.depth() != sizeof(uchar));
	cv::Mat_<cv::Vec3b> _I = inImg;
	cv::Mat_<cv::Vec3b> Iout = outImg;

	int mask_size = sqrt(filter.size());
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
				new_pix_value += filter[i] * vec[i];
			}

			Iout(i, j)[2] = new_pix_value;
			Iout(i, j)[1] = new_pix_value;
			Iout(i, j)[0] = new_pix_value;
		}
}

void draw_bounds(cv::Mat& inImg, std::vector<int> bounds) {
	CV_Assert(inImg.depth() != sizeof(uchar));
	cv::Mat_<cv::Vec3b> _I = inImg;
	for (int i = 0; i < _I.rows; ++i)
		for (int j = 0; j < _I.cols; ++j) {
			if(i == bounds[2] && j > bounds[0] && j < bounds[1] ||
				i == bounds[3] && j > bounds[0] && j < bounds[1] ||
				j == bounds[0] && i > bounds[2] && i < bounds[3] ||
				j == bounds[1] && i > bounds[2] && i < bounds[3]) {
				_I(i, j)[2] = 0;
				_I(i, j)[1] = 0;
				_I(i, j)[0] = 255;
			}
		}

}

void rotate_by(cv::Mat& inImg, cv::Mat& outImg, int angle) {

}

int main(int argc, char* argv[]) {
	if (argc == 2 && (cv::imread(argv[1]).data != NULL)) {
		cv::Mat image = cv::imread(argv[1]);
		cv::Mat tmp = cv::imread(argv[1]);
		std::vector<std::vector<int>> horizontal_lines;
		std::vector<int> logo_bounds;
		treshold(image, tmp);
		cv::Mat tmp2 = tmp.clone();
		//find_horizontal_lines(tmp, tmp2);
		horizontal_lines = find_horizontal_lines(tmp, tmp2);
		logo_bounds = get_logos_possible_bounds(horizontal_lines);
		draw_bounds(image, logo_bounds);

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
