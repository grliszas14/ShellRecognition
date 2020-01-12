#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/opencv.hpp"  //getRotationMatrix2D, warpAffine

// RGB
const int RED_TRESH[3] = {100, 80, 80};
const int YELLOW_TRESH[3] = {150,90,80};
const int MIN_THICKNESS = 2;

void treshold(cv::Mat& inImg, cv::Mat& outImg) {
	CV_Assert(inImg.depth() != sizeof(uchar));
	cv::Mat_<cv::Vec3b> _I = inImg;
	cv::Mat_<cv::Vec3b> Iout = outImg;

	for( int i = 0; i < (inImg.rows); ++i)
		for (int j = 0; j < (inImg.cols); ++j) {
			Iout(i, j)[0] = 0;
			Iout(i, j)[1] = 0;
			Iout(i, j)[2] = 0;
			if (_I(i, j)[2] > RED_TRESH[0] && _I(i, j)[1] < RED_TRESH[1] &&  _I(i, j)[0] < RED_TRESH[2]) {
				Iout(i, j)[2] = 255;
				Iout(i, j)[1] = 255;
				Iout(i, j)[0] = 255;
			}
		}
}

double check_yellow_percentage(cv::Mat& inImg) {
	CV_Assert(inImg.depth() != sizeof(uchar));
	cv::Mat_<cv::Vec3b> _I = inImg;
	double count_yellow_pixels = 0;
	double pixels = inImg.rows * inImg.cols;

	for( int i = 0; i < (inImg.rows); ++i)
		for (int j = 0; j < (inImg.cols); ++j) {
			//std::cout << "R:" << int(_I(i, j)[2]) << " G:" << int(_I(i, j)[1]) << " B:" << int(_I(i, j)[0]) << std::endl;
			if (_I(i, j)[2] > YELLOW_TRESH[0] && _I(i, j)[1] > YELLOW_TRESH[1] &&  _I(i, j)[0] < YELLOW_TRESH[2]) {
				count_yellow_pixels++;
			}
		}
	/*  TO SAMO W HSV
	cv::Mat HSV;
	cv::cvtColor(inImg, HSV, CV_BGR2HSV);
	for( int i = 0; i < (HSV.rows); ++i)
		for (int j = 0; j < (HSV.cols); ++j) {
			//std::cout << "R:" << int(_I(i, j)[2]) << " G:" << int(_I(i, j)[1]) << " B:" << int(_I(i, j)[0]) << std::endl;
			if (_I(i, j)[0] > 28 && _I(i, j)[0] > 33) {
				count_yellow_pixels++;
			}
		}
*/
	//std::cout << "Yellow pixels: " << count_yellow_pixels << std::endl;
	//std::cout << "All pixels: " << pixels << std::endl;
	return count_yellow_pixels / pixels;
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

std::vector<int> get_logos_possible_borders(std::vector<std::vector<int>> lines,
		std::vector<std::vector<int>> &visited, int &last_size, bool &processed_angle) {
	//for (int i = 0; i < lines.size(); ++i) {
	//	std::cout << "I: " << lines[i][0] << " Begin: " << lines[i][1] << " Length of line: " << lines[i][2] << std::endl;
	//}

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
		int i = 0;
		for (i = 0; i < lines.size(); ++i) {
			if (start_flag == -1) {
				std::vector<int> tmp = {lines[i][0], lines[i][1], lines[i][2]};
				if (std::find(visited.begin(), visited.end(), tmp) != visited.end()) {
					continue;
				}

				tmp_length = lines[i][2];
				tmp_begin = lines[i][1];
				tmp_row = lines[i][0];

				visited.push_back(tmp);
				tmp_thickness++;
				start_flag = 1;
				continue;
			}
			if (lines[i][2] > tmp_length && lines[i][0] < tmp_row + 15
					&& lines[i][1] < tmp_begin && lines[i][1] > (tmp_begin - tmp_begin*6/100)) {
				std::vector<int> tmp = {lines[i][0], lines[i][1], lines[i][2]};
				if (std::find(visited.begin(), visited.end(), tmp) != visited.end()) {
					continue;
				}
				tmp_length = lines[i][2];
				tmp_begin = lines[i][1];
				tmp_row = lines[i][0];
				visited.push_back(tmp);
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
	//if (visited.size() == 6) {
	//	processed_angle = true;
	//}
	if (last_size == visited.size()) {
		//std::cout << "END this angle" << std::endl;
		processed_angle = true;
	}

	last_size = visited.size();
	std::vector<int> b_rbl_tr;
	if (tmp_thickness < MIN_THICKNESS) {
		return b_rbl_tr;
	}
	//std::cout << " " << std::endl;
	//std::cout << "I: " << tmp_row << " Begin: " << tmp_begin << " Length of line: "
		//<< tmp_length << " Thickness: " << tmp_thickness << std::endl;

	int probable_bottom_line = tmp_thickness * 11;
	int lower_row_bound = tmp_row - tmp_thickness + probable_bottom_line - probable_bottom_line*10/100;
	int upper_row_bound = tmp_row - tmp_thickness + probable_bottom_line + probable_bottom_line*10/100;
	//std::cout << "Lower row bound: " << lower_row_bound << std::endl;
	//std::cout << "Upper row bound: " << upper_row_bound << std::endl;
	// bottom_row_begin_length_tmp_row
	for (int i = 0; i < lines.size()-1; ++i) {
		int double_line = 0;
		if (lines[i][0] == lines[i+1][0] && lines[i+1][1] < (lines[i][1] + int(lines[i][2] * 1.5))) {
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
			b_rbl_tr.push_back(double_line);
			b_rbl_tr.push_back(tmp_row);

			//std::cout << "I: " << lines[i][0] << " Begin: " << lines[i][1] << " Length of line: " << lines[i][2] << std::endl;
			return b_rbl_tr;
		}
	}
}

std::vector<int> calculate_bounds(std::vector<int> factors) {
	// Calculate borders
	int width = factors[2];
	int height = factors[0] - factors[3];
	int left_x = factors[1] - width*35/100;
	int right_x = factors[1] + factors[2] + width*35/100;
	int upper_y = factors[3] - height*12/100;
	int lower_y = factors[0] + height*6/100;
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

void fill_black(cv::Mat& inImg, std::vector<int> bounds) {
	CV_Assert(inImg.depth() != sizeof(uchar));
	cv::Mat_<cv::Vec3b> _I = inImg;
	for (int i = 0; i < _I.rows; ++i)
		for (int j = 0; j < _I.cols; ++j) {
			if ( i > bounds[2] && i < bounds[3] && j > bounds[0] && j < bounds[1]) {
				_I(i, j)[2] = 0;
				_I(i, j)[1] = 0;
				_I(i, j)[0] = 0;
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

double calculate_moment(cv::Mat& I, int p, int q) {
	double mpq = 0;
	cv::Mat_<cv::Vec3b> _I = I;
	for (int i = 0; i < I.rows; ++i)
		for (int j = 0; j < I.cols; ++j) {
			if ((_I(i, j)[0] + _I(i, j)[1] + _I(i, j)[2]) / 3 > 127) {
				mpq = mpq + (pow(i, p) * pow(j, q) * 1);
			}
		}

	return mpq;
}

double calculate_M11(double m11, double m10, double m01, double m00) {
	double M11;
	M11 = m11 - (m10 * m01 / m00);
	return M11;
}

double calculate_M12(double m12, double m11, double m01, double m00, double m02, double m10) {
	double M12;
	M12 = m12 - 2 * m11 * m01 / m00 - m02 * m10 / m00 + 2 * m10 * pow(m01/m00, 2);
	return M12;
}

double calculate_M21(double m21, double m11, double m10, double m00, double m20, double m01) {
	double M21;
	M21 = m21 - 2 * m11 * m10 / m00 - m20 * m01 / m00 + 2 * m01 * pow(m10/m00, 2);
	return M21;
}

double calculate_M02(double m02, double m01, double m00) {
	double M02;
	M02 = m02 - (m01 * m01 / m00);
	return M02;
}

double calculate_M20(double m20, double m10, double m00) {
	double M20;
	M20 = m20 -( m10 * m10 / m00);
	return M20;
}

double calculate_M30(double m30, double m20, double m10, double m00) {
	double M30;
	M30 = m30 - 3 * m20 * m10 / m00 + 2 * m10 * pow(m10/m00, 2);
	return M30;
}

double calculate_M03(double m03, double m02, double m01, double m00) {
	double M03;
	M03 = m03 - 3 * m02 * m01 / m00 + 2 * m01 * pow(m01/m00, 2);
	return M03;
}

double calculate_M1(cv::Mat& I) {
	double M1;
	double m02, m01, m00, m20, m10, M20, M02;
	m02 = calculate_moment(I, 0, 2);
	m01 = calculate_moment(I, 0, 1);
	m00 = calculate_moment(I, 0, 0);
	m20 = calculate_moment(I, 2, 0);
	m10 = calculate_moment(I, 1, 0);
	M20 = calculate_M20(m20, m10, m00);
	M02 = calculate_M02(m02, m01, m00);
	M1 = (M20 + M02) / (m00*m00);
	return M1;
}

double calculate_M2(cv::Mat& I) {
	double M2;
	double M20, M02, M11, m20, m10, m02, m01, m11, m00;

	m20 = calculate_moment(I, 2, 0);
	m10 = calculate_moment(I, 1, 0);
	m02 = calculate_moment(I, 0, 2);
	m01 = calculate_moment(I, 0, 1);
	m11 = calculate_moment(I, 1, 1);
	m00 = calculate_moment(I, 0, 0);
	M11 = calculate_M11(m11, m10, m01, m00);
	M02 = calculate_M02(m02, m01, m00);
	M20 = calculate_M20(m20, m10, m00);
	M2 = (pow(M20 - M02, 2) + 4 * pow(M11, 2)) / pow(m00,4);
	return M2;
}

double calculate_M3(cv::Mat& I) {
	double M3;
	double M30, M12, M21, M03, m21, m03, m02, m12, m11, m01, m30, m20, m10, m00;

	m03 = calculate_moment(I, 0, 3);
	m02 = calculate_moment(I, 0, 2);
	m12 = calculate_moment(I, 1, 2);
	m11 = calculate_moment(I, 1, 1);
	m01 = calculate_moment(I, 0, 1);
	m30 = calculate_moment(I, 3, 0);
	m20 = calculate_moment(I, 2, 0);
	m10 = calculate_moment(I, 1, 0);
	m00 = calculate_moment(I, 0, 0);
	m21 = calculate_moment(I, 2, 1);
	M03 = calculate_M03(m03, m02, m01, m00);
	M21 = calculate_M21(m21, m11, m10, m00, m20, m01);
	M12 = calculate_M12(m12, m11, m01, m00, m02, m10);
	M30 = calculate_M30(m30, m20, m10, m00);
	M3 = (pow(M30 - 3 * M12, 2) + pow(3 * M21 - M03, 2)) / pow(m00,5);
	return M3;
}

double calculate_M4(cv::Mat& I) {
	double M4;
	double M30, M12, M21, M03, m21, m03, m02, m12, m11, m01, m30, m20, m10, m00;

	m03 = calculate_moment(I, 0, 3);
	m02 = calculate_moment(I, 0, 2);
	m12 = calculate_moment(I, 1, 2);
	m11 = calculate_moment(I, 1, 1);
	m01 = calculate_moment(I, 0, 1);
	m30 = calculate_moment(I, 3, 0);
	m20 = calculate_moment(I, 2, 0);
	m10 = calculate_moment(I, 1, 0);
	m00 = calculate_moment(I, 0, 0);
	m21 = calculate_moment(I, 2, 1);
	M03 = calculate_M03(m03, m02, m01, m00);
	M21 = calculate_M21(m21, m11, m10, m00, m20, m01);
	M12 = calculate_M12(m12, m11, m01, m00, m02, m10);
	M30 = calculate_M30(m30, m20, m10, m00);
	M4 = (pow(M30 + M12, 2) + pow(M21 + M03, 2)) / pow(m00,5);
	return M4;
}

double calculate_M5(cv::Mat& I) {
	double M5;
	double M30, M12, M21, M03, m21, m03, m02, m12, m11, m01, m30, m20, m10, m00;

	m03 = calculate_moment(I, 0, 3);
	m02 = calculate_moment(I, 0, 2);
	m12 = calculate_moment(I, 1, 2);
	m11 = calculate_moment(I, 1, 1);
	m01 = calculate_moment(I, 0, 1);
	m30 = calculate_moment(I, 3, 0);
	m20 = calculate_moment(I, 2, 0);
	m10 = calculate_moment(I, 1, 0);
	m00 = calculate_moment(I, 0, 0);
	m21 = calculate_moment(I, 2, 1);
	M03 = calculate_M03(m03, m02, m01, m00);
	M21 = calculate_M21(m21, m11, m10, m00, m20, m01);
	M12 = calculate_M12(m12, m11, m01, m00, m02, m10);
	M30 = calculate_M30(m30, m20, m10, m00);


	M5 = ((M30 - 3 * M12) * (M30 + M12) * (pow(M30 + M12, 2) - 3 * pow(M21 + M03, 2)) + (3 * M21 - M03) * (M21 + M03) * (3 * pow(M30 + M12, 2) - pow(M21 + M03, 2))) / pow(m00,10);
	return M5;
}

double calculate_M6(cv::Mat& I) {
	double M6;
	double M30, M12, M21, M03, M11, M20, M02, m21, m03, m02, m12, m11, m01, m30, m20, m10, m00;

	m03 = calculate_moment(I, 0, 3);
	m02 = calculate_moment(I, 0, 2);
	m12 = calculate_moment(I, 1, 2);
	m11 = calculate_moment(I, 1, 1);
	m01 = calculate_moment(I, 0, 1);
	m30 = calculate_moment(I, 3, 0);
	m20 = calculate_moment(I, 2, 0);
	m10 = calculate_moment(I, 1, 0);
	m00 = calculate_moment(I, 0, 0);
	m21 = calculate_moment(I, 2, 1);
	M03 = calculate_M03(m03, m02, m01, m00);
	M21 = calculate_M21(m21, m11, m10, m00, m20, m01);
	M12 = calculate_M12(m12, m11, m01, m00, m02, m10);
	M30 = calculate_M30(m30, m20, m10, m00);
	M11 = calculate_M11(m11, m10, m01, m00);
	M02 = calculate_M02(m02, m01, m00);
	M20 = calculate_M20(m20, m10, m00);

	M6 = ((M20 - M02) * (pow(M30 + M12, 2) - pow(M21 + M03, 2) ) + 4 * M11 * (M30 + M12) * (M21 + M03)) / pow(m00, 7);
	return M6;
}

double calculate_M7(cv::Mat& I) {
	double M7, m02, m01, m00, m20, m10, m11, M20, M02, M11;
	m02 = calculate_moment(I, 0, 2);
	m01 = calculate_moment(I, 0, 1);
	m00 = calculate_moment(I, 0, 0);
	m20 = calculate_moment(I, 2, 0);
	m10 = calculate_moment(I, 1, 0);
	m11 = calculate_moment(I, 1, 1);
	M20 = calculate_M20(m20, m10, m00);
	M02 = calculate_M02(m02, m01, m00);
	M11 = calculate_M11(m11, m10, m01, m00);
	M7 = (M02 * M20 - M11 * M11) / powl(m00, 4);

	return M7;
}

std::vector<int> calculate_left_down_square(std::vector<int> bounds) {
	std::vector<int> ld_square;
	int minX = bounds[0];
	int maxX = bounds[0] + ((bounds[1] - bounds[0]) / 2);
	int minY = bounds[3] - ((bounds[3] - bounds[2]) / 2);
	int maxY = bounds[3];

	ld_square.push_back(minX);
	ld_square.push_back(maxX);
	ld_square.push_back(minY);
	ld_square.push_back(maxY);

	return ld_square;
}

std::vector<int> calculate_right_down_square(std::vector<int> bounds) {
	std::vector<int> rd_square;
	int minX = bounds[1] - ((bounds[1] - bounds[0]) / 2);
	int maxX = bounds[1];
	int minY = bounds[3] - ((bounds[3] - bounds[2]) / 2);
	int maxY = bounds[3];

	rd_square.push_back(minX);
	rd_square.push_back(maxX);
	rd_square.push_back(minY);
	rd_square.push_back(maxY);

	return rd_square;
}

std::vector<int> calculate_upper_half_rect(std::vector<int> bounds) {
	std::vector<int> rect;
	int minX = bounds[0];
	int maxX = bounds[1];
	int minY = bounds[2];
	int maxY = bounds[2] + ((bounds[3] - bounds[2]) / 2);

	rect.push_back(minX);
	rect.push_back(maxX);
	rect.push_back(minY);
	rect.push_back(maxY);

	return rect;
}

int main(int argc, char* argv[]) {
	if (argc == 3 && (cv::imread(argv[1]).data != NULL)) {
		int iterations = std::stoi(argv[2]);
		cv::Mat image = cv::imread(argv[1]);
		cv::Mat result = cv::imread(argv[1]);
		cv::Mat rotated;
		//rotate_by(image, rotated, -5*4);
		//cv::imwrite("proba9.png", rotated);
		bool found = false;
		std::vector<std::vector<int>> visited;

///*
		// Processing
		for (int iter = 0; iter < iterations; ++iter) {
			std::cout << iter << std::endl;
			std::vector<std::vector<int>> horizontal_lines;
			std::vector<int> logo_borders;
			std::vector<int> logo_bounds;

			// Tu uporzadkowac
			cv::Mat rotated;
			rotate_by(image, rotated, 5*iter);
			rotate_by(result, result, 5*iter);
			cv::Mat processed = rotated.clone();
			cv::Mat tresh = rotated.clone();
			treshold(rotated, tresh);
			//cv::imshow("Shell_tresh", tresh);
			cv::imwrite("tresh.png", tresh);
			cv::Mat tmp2 = rotated.clone();

			horizontal_lines = find_horizontal_lines(tresh, tmp2);

			bool processed_angle = false;
			int last_size = visited.size();
			while (processed_angle != true) {
				logo_borders = get_logos_possible_borders(horizontal_lines, visited, last_size, processed_angle);
				if (!logo_borders.empty()) {
					logo_bounds = calculate_bounds(logo_borders);

					std::vector<int> ld_square = calculate_left_down_square(logo_bounds);
					cv::Mat left_down_corner(tresh, cv::Range(ld_square[2], ld_square[3]), cv::Range(ld_square[0], ld_square[1]));
					double ldc_M1 = calculate_M1(left_down_corner);
					double ldc_M7 = calculate_M7(left_down_corner);
					//std::cout << "Left down M1: " << ldc_M1 << std::endl;
					//std::cout << "Left down M7: " << ldc_M7 << std::endl;
					if (ldc_M1 < 0.40 || ldc_M1 > 0.76) continue;
					if (ldc_M7 < 0.025 || ldc_M7 > 0.065) continue;


					std::vector<int> rd_square = calculate_right_down_square(logo_bounds);
					cv::Mat right_down_corner(tresh, cv::Range(rd_square[2], rd_square[3]), cv::Range(rd_square[0], rd_square[1]));
					double rdc_M1 = calculate_M1(right_down_corner);
					double rdc_M7 = calculate_M7(right_down_corner);
					//std::cout << "Right down M1: " << rdc_M1 << std::endl;
					//std::cout << "Right down M7: " << rdc_M7 << std::endl;
					if (rdc_M1 < 0.40 || rdc_M1 > 0.76) continue;
					if (rdc_M7 < 0.025 || rdc_M7 > 0.065) continue;

					std::vector<int> uh_rect = calculate_upper_half_rect(logo_bounds);
					cv::Mat upper_half(tresh, cv::Range(uh_rect[2], uh_rect[3]), cv::Range(uh_rect[0], uh_rect[1]));
					double uh_M1 = calculate_M1(upper_half);
					//std::cout << "Upper half M1: " << uh_M1 << std::endl;
					if (uh_M1 < 0.6) continue;

					cv::Mat logo(processed, cv::Range(uh_rect[2], rd_square[3]), cv::Range(uh_rect[0], uh_rect[1]));
					double yellow_percentage = check_yellow_percentage(logo);
					//std::cout << "Yellow percentage: " << yellow_percentage << std::endl;
					if (yellow_percentage < 0.25) continue;

					found = true;

					draw_bounds(result, logo_bounds);
					fill_black(processed, logo_bounds);
				}
			}
			rotate_by(processed, processed, -5*iter);
			rotate_by(result, result, -5*iter);
			image = processed.clone();

		}
//		*/
		if (found) {
			std::cout << "Shell logo found!" << std::endl;
		} else {
			std::cout << "Shell logo not found!" << std::endl;
		}
		// Show results
		cv::imshow("Result", result);
		cv::imwrite("result.png", result);
		cv::waitKey(-1);
		return 0;
	} else {
		std::cout << "Usage: shell_recognition <file_name> <iterations>" << std::endl;

		return -1;
	}
}
