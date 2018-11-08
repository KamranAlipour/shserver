#pragma once
#ifndef OPENCV_LIBS

#define OPENCV_LIBS
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/videoio.hpp> // to load video capture
#include <opencv2/imgproc.hpp> // to resize input video stream
//#include <opencv2/aruco.hpp> // Marker Tracking
#include <opencv2/calib3d.hpp> // for Rodrigues function in Aruco ViewMatrix Calculation

#endif // !OPENCV_LIBS