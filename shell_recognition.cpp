#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/opencv.hpp"  //getRotationMatrix2D, warpAffine

// BGR
// TODO pokombinowac jeszcze z odcieniem 
const int RED_TRESH[3] = {80, 80, 90};
// TODO
const int YELLOW_TRESH[3] = {0,0,0};

const std::vector<int> HORIZONTAL_FILTER = {-1,-1,-1,
											 2, 2, 2,
											-1,-1,-1};

const std::vector<int> HORIZONTAL_FILTER2 = {2,2,2,2,2,2,2};
const int MIN_THICKNESS = 3;

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
				if (line_length > 5) {
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

std::vector<int> get_logos_possible_borders(std::vector<std::vector<int>> lines) {
	for (int i = 0; i < lines.size(); ++i) {
		std::cout << "I: " << lines[i][0] << " Begin: " << lines[i][1] << " Length of line: " << lines[i][2] << std::endl;
	}

	int start_flag = -1;
	int tmp_length = 0;
	int tmp_begin = 0;
	int tmp_row = 0;
	int tmp_thickness = 0;

	int iter = 0;
	int search_end = 0;
	while (search_end != 1) {
		if (iter >= lines.size()) {
			search_end = 1;
		}
		std::cout << "Loop" << std::endl;
		int i = 0;
		for (i = 0; i < lines.size(); ++i) {
			if (start_flag == -1) {
				tmp_length = lines[i][2];
				tmp_begin = lines[i][1];
				tmp_row = lines[i][0];
				tmp_thickness++;
				start_flag = 1;
				continue;
			}
			if (lines[i][2] > tmp_length && lines[i][0] < tmp_row + 5
					&& lines[i][1] < tmp_begin && lines[i][1] > (tmp_begin - tmp_begin*6/100)) {
				tmp_length = lines[i][2];
				tmp_begin = lines[i][1];
				tmp_row = lines[i][0];
				tmp_thickness++;
			}
		}
		if (tmp_thickness > MIN_THICKNESS) {
		   	search_end = 1;
		} else {
			iter = i+1;
			start_flag = 0;
		}
		// ogarnac przypadek kiedy na jednej linii sa dwa prawdopodobne loga
		// (czyli tu else w trakcie, albo kolejnego szukac w liniach ponizej)

	}
	std::vector<int> b_rbl_tr;
	if (tmp_thickness < MIN_THICKNESS) {
		return b_rbl_tr;
	}
	std::cout << " " << std::endl;
	std::cout << "I: " << tmp_row << " Begin: " << tmp_begin << " Length of line: "
		<< tmp_length << " Thickness: " << tmp_thickness << std::endl;

	int probable_bottom_line = tmp_thickness * 11;
	int lower_row_bound = tmp_row - tmp_thickness + probable_bottom_line - probable_bottom_line*10/100;
	int upper_row_bound = tmp_row - tmp_thickness + probable_bottom_line + probable_bottom_line*10/100;
	std::cout << "Lower row bound: " << lower_row_bound << std::endl;
	std::cout << "Upper row bound: " << upper_row_bound << std::endl;
	// bottom_row_begin_length_tmp_row
	for (int i = 0; i < lines.size()-1; ++i) {
		int double_line = 0;
		if (lines[i][0] == lines[i+1][0] && lines[i+1][1] < (lines[i][1] + lines[i][2] * 1.5)) {
			double_line = lines[i][2] + lines[i+1][2];
		} else {
			double_line = lines[i][2];
		}
		if (lines[i][0] > (lower_row_bound) /*  && lines[i][0] < (upper_row_bound)*/ &&
			   	lines[i][1] > (tmp_begin - tmp_begin*20/100) &&
				lines[i][1] < (tmp_begin + tmp_begin*10/100) &&
			   	double_line > tmp_length) {
			b_rbl_tr.push_back(lines[i][0]);
			b_rbl_tr.push_back(lines[i][1]);
			b_rbl_tr.push_back(lines[i][2]);
			b_rbl_tr.push_back(tmp_row);

			std::cout << "I: " << lines[i][0] << " Begin: " << lines[i][1] << " Length of line: " << lines[i][2] << std::endl;
			return b_rbl_tr;
		}
	}
}

std::vector<int> calculate_bounds(std::vector<int> factors) {
	// Calculate borders
	int left_x = factors[1] - factors[1]*30/100;
	int right_x = factors[1] + factors[2] + factors[1]*30/100;
	int upper_y = factors[3] - factors[3]*15/100;
	int lower_y = factors[0] + factors[0]*6/100;
	std::vector<int> bounds = {left_x, right_x, upper_y, lower_y};
	return bounds;
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
	CV_Assert(inImg.depth() != sizeof(uchar));
	cv::Mat src = inImg;

	// get rotation matrix for rotating the image around its center in pixel coordinates
	cv::Point2f center((src.cols-1)/2.0, (src.rows-1)/2.0);
	cv::Mat rot = cv::getRotationMatrix2D(center, angle, 1.0);

	// determine bounding rectangle, center not relevant
	cv::Rect2f bbox = cv::RotatedRect(cv::Point2f(), src.size(), angle).boundingRect2f();

	// adjust transformation matrix
	rot.at<double>(0,2) += bbox.width/2.0 - src.cols/2.0;
	rot.at<double>(1,2) += bbox.height/2.0 - src.rows/2.0;

	//cv::Mat dst;
	cv::warpAffine(src, outImg, rot, bbox.size());
	//cv::imwrite("rotated_im.png", dst);
}

int main(int argc, char* argv[]) {
	if (argc == 2 && (cv::imread(argv[1]).data != NULL)) {
		cv::Mat image = cv::imread(argv[1]);
		///cv::Mat processed = cv::imread(argv[1]);
		//cv::Mat tmp = cv::imread(argv[1]);
		cv::Mat rotated;
		//rotate_by(image, rotated, 5*2);
		//cv::imwrite("proba2.png", rotated);


		// Processing
		for (int iter = 0; iter < 1; ++iter) {
			std::vector<std::vector<int>> horizontal_lines;
			std::vector<int> logo_borders;
			std::vector<int> logo_bounds;

			// Tu uporzadkowac
			cv::Mat rotated;
			rotate_by(image, rotated, 5*iter);
			cv::Mat processed = rotated.clone();
			cv::Mat tmp = rotated.clone();
			treshold(rotated, tmp);
			cv::imshow("Shell_tresh", tmp);
			cv::imwrite("tresh.png", tmp);
			cv::Mat tmp2 = rotated.clone();

			horizontal_lines = find_horizontal_lines(tmp, tmp2);
			logo_borders = get_logos_possible_borders(horizontal_lines);
			std::cout << "logo_borders size: " << logo_borders.size() << std::endl;
			if (!logo_borders.empty()) {
				logo_bounds = calculate_bounds(logo_borders);
				draw_bounds(processed, logo_bounds);
				rotate_by(processed, processed, -5*iter);
				cv::imshow("Obrysowane znalezione logo", processed);
				cv::imwrite("result.png", processed);
				break;
			}
		}

		// Show results
		//cv::imshow("Oryginal", image);
		//cv::imshow("Shell_tresh", tmp);
		//cv::imshow("Shell_horizontal", tmp2);
		cv::waitKey(-1);
		return 0;
	} else {
		std::cout << "Wrong file name!" << std::endl;
		std::cout << "Usage: shell_recognition <file_name>" << std::endl;

		return -1;
	}
}


// nieoptymalne, zmienic na greyscale
// nazwa: apply_filter? dodatkowy argument?
/*  
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
*/
