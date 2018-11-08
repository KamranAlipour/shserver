#include <Poco/Net/HTTPClientSession.h>
#include <Poco/Net/HTTPRequest.h>
#include <Poco/Net/HTTPResponse.h>
#include <Poco/Net/Context.h>
#include <Poco/Net/SSLManager.h>
#include <Poco/Net/SSLException.h>
#include <Poco/URI.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp> // to resize input video stream
#include <opencv2/calib3d.hpp>

using namespace cv;
using namespace std;
using namespace Poco;
using namespace Net;

#define THETA_FRAME_SIZE 550000

class RicohTheta
{
	uint8_t Image[THETA_FRAME_SIZE];
	Mat rawData;
	Mat img;
	streambuf* sb;
	string thetaUrl;
	string sessionId;
	string url;
	URI uri;
	HTTPClientSession session;
	HTTPRequest req; // request
	string reqBody; //request body
	Poco::Net::HTTPResponse res; // response
public:
	RicohTheta(string URL)
	{
		string executeCmd = "/osc/commands/execute";
		rawData = Mat(1, 500000, CV_8UC1, Image);
		thetaUrl = URL;
		url = thetaUrl + executeCmd;
		uri = Poco::URI(url);
		session.setHost(uri.getHost());
		session.setPort(uri.getPort());
		req.setMethod(HTTPRequest::HTTP_POST);
		req.setURI(uri.getPathAndQuery());
		req.setVersion(HTTPMessage::HTTP_1_1);
		req.setContentType("application/json");
		req.setKeepAlive(true);
	}

	void startSession()
	{
		printf("Starting the Theta Cam Connection ...\n");
		while (true)
		{
			try
			{
				reqBody = "{ \"name\" : \"camera.startSession\", \"parameters\": {} }";
				req.setContentLength(reqBody.length());

				session.sendRequest(req) << reqBody;

				std::istream& rs = session.receiveResponse(res);

				string resp(istreambuf_iterator<char>(rs), {});

				size_t sesIDpo = resp.find("\"results\":{ \"sessionId\":\"");
				sessionId = resp.substr(sesIDpo + 25, 8);
				return;
			}
			catch (const Poco::Net::SSLException& e)
			{
				std::cerr << e.what() << ": " << e.message() << std::endl;
			}
			catch (const std::exception& e)
			{
				std::cerr << e.what() << std::endl;
			}
		}
	}

	void startLivePreview()
	{
		printf("Starting Live Preview from Theta Cam ...\n");
		while (true)
		{
			try
			{

				reqBody = "{ \"name\": \"camera._getLivePreview\", \"parameters\": { \"sessionId\": \"" + sessionId + "\" } }";
				req.setContentLength(reqBody.length());

				session.sendRequest(req) << reqBody;

				std::istream& rs = session.receiveResponse(res);
				sb = rs.rdbuf();
				return;
			}
			catch (const Poco::Net::SSLException& e)
			{
				std::cerr << e.what() << ": " << e.message() << std::endl;
			}
			catch (Poco::TimeoutException& e)
			{
				startSession();
				continue;
			}
		}
		//catch (const std::exception& e)
		//{
		//	std::cerr << e.what() << std::endl;
		//}
		//updateLiveFrame();
	}

	void updateLiveFrame()
	{
		uint8_t ch1 = 0;
		uint8_t ch2 = 0;

		bool isLoadStart = false;
		int counter = 0;
		while (true)
		{
			try
			{
				ch1 = sb->sbumpc();
				ch2 = sb->sbumpc();
			}
			catch (Poco::Net::ConnectionAbortedException e)
			{
				startLivePreview();
				return;
			}
			catch (Poco::Net::ConnectionResetException e)
			{
				startSession();
				startLivePreview();
				return;
			}
			
			if (!isLoadStart)
			{
				if (ch1 == 0xFF && ch2 == 0xD8)
				{
					// mjpeg start! ( [0xFF 0xD8 ... )
					Image[counter++] = ch1;
					Image[counter++] = ch2;

					isLoadStart = true;
				}
			}
			else
			{
				//frame.data[counter++] = (counter % 3 == 2) ? 0xFF : 0x0;
				Image[counter++] = ch1;
				Image[counter++] = ch2;

				if (ch1 == 0xFF && ch2 == 0xD9)
				{
					img = imdecode(rawData, IMREAD_ANYCOLOR);
					counter = 0;
					isLoadStart = false;
					return;
				}
			}
		}
	}

	Mat getLiveFrame()
	{
		return img;
	}

	/*
	int main()
	{
	string thetaUrl = "http://192.168.1.1:80";
	string sessionID = startSession(thetaUrl);
	getLivePreview(thetaUrl, sessionID);
	return 0;
	}
	*/
};
