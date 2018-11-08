#include "client_http.hpp"
#include "server_http.hpp"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/filesystem.hpp>

using namespace std;
// Added for the json-example:
using namespace boost::property_tree;

using HttpServer = SimpleWeb::Server<SimpleWeb::HTTP>;
using HttpClient = SimpleWeb::Client<SimpleWeb::HTTP>;

#include "SphericalHarmonics.h"
using namespace cv;

int main(int argc, char* argv[])
{
	char test_flag[] = "image"; // a flag which defines if the algorithm should use sample test image as input instead of the wifi theta cam stream
	bool is_test = false;
	
	HttpServer server_interface;
	server_interface.config.port = 6060;

	HttpServer server;
	server.config.port = 8080;

	float LSR00 = 1.0; // Red Light Scale 00
	float LSR = 1.0; // Red Light Scale all except 00
	float LSG00 = 1.0; // Green Light Scale 00
	float LSG = 1.0; // Green Light Scale all except 00
	float LSB00 = 1.0; // Blue Light Scale 00
	float LSB = 1.0; // Blue Light Scale all except 00

	if (argv[1] && strcmp(test_flag, argv[1]) == 0)
	{
		printf("running the test sample input image...\n");
		is_test = true;
	}
	server_interface.default_resource["GET"] = [](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request)
	{
		try {
			auto web_root_path = boost::filesystem::canonical("web");
			//cout << "WEB ROOT PATH: " << web_root_path << endl;
			auto path = boost::filesystem::canonical(web_root_path / request->path);
			// Check if path is within web_root_path
			if (distance(web_root_path.begin(), web_root_path.end()) > distance(path.begin(), path.end()) ||
				!equal(web_root_path.begin(), web_root_path.end(), path.begin()))
				throw invalid_argument("path must be within root path");
			if (boost::filesystem::is_directory(path))
				path /= "index.html";
			//cout << "INDEX HTML PATH: " << path << endl;
			SimpleWeb::CaseInsensitiveMultimap header;

			// Uncomment the following line to enable Cache-Control
			// header.emplace("Cache-Control", "max-age=86400");

#ifdef HAVE_OPENSSL
			//    Uncomment the following lines to enable ETag
			//    {
			//      ifstream ifs(path.string(), ifstream::in | ios::binary);
			//      if(ifs) {
			//        auto hash = SimpleWeb::Crypto::to_hex_string(SimpleWeb::Crypto::md5(ifs));
			//        header.emplace("ETag", "\"" + hash + "\"");
			//        auto it = request->header.find("If-None-Match");
			//        if(it != request->header.end()) {
			//          if(!it->second.empty() && it->second.compare(1, hash.size(), hash) == 0) {
			//            response->write(SimpleWeb::StatusCode::redirection_not_modified, header);
			//            return;
			//          }
			//        }
			//      }
			//      else
			//        throw invalid_argument("could not read file");
			//    }
#endif
			auto ifs = make_shared<ifstream>();
			ifs->open(path.string(), ifstream::in | ios::binary | ios::ate);

			if (*ifs) {
				auto length = ifs->tellg();
				ifs->seekg(0, ios::beg);

				header.emplace("Content-Length", to_string(length));
				response->write(header);

				// Trick to define a recursive function within this scope (for example purposes)
				class FileServer {
				public:
					static void read_and_send(const shared_ptr<HttpServer::Response> &response, const shared_ptr<ifstream> &ifs) {
						// Read and send 128 KB at a time
						static vector<char> buffer(131072); // Safe when server is running on one thread
						streamsize read_length;
						if ((read_length = ifs->read(&buffer[0], static_cast<streamsize>(buffer.size())).gcount()) > 0) {
							response->write(&buffer[0], read_length);
							if (read_length == static_cast<streamsize>(buffer.size())) {
								response->send([response, ifs](const SimpleWeb::error_code &ec) {
									if (!ec)
										read_and_send(response, ifs);
									else
										cerr << "Connection interrupted" << endl;
								});
							}
						}
					}
				};
				FileServer::read_and_send(response, ifs);
			}
			else
				throw invalid_argument("could not read file");
		}
		catch (const exception &e) {
			response->write(SimpleWeb::StatusCode::client_error_bad_request, "Could not open path " + request->path + ": " + e.what());
		}
	};

	server_interface.resource["^/GetLightFactors$"]["GET"] = [&LSR00, &LSG00, &LSB00, &LSR, &LSG, &LSB](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request)
	{
		stringstream stream;

		stream << "{\"LSR00\":"<< LSR00 << ", \"LSG00\":" << LSG00 << ", \"LSB00\":" << LSB00 << ", \"LSR\":" << LSR << ", \"LSG\":" << LSG << ", \"LSB\":" << LSB << "}";

		response->write(stream);
	};
	server.default_resource["GET"] = [&LSR00, &LSG00, &LSB00, &LSR, &LSG, &LSB](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> /*request*/)
	{
		char sb[500];
		sprintf(sb, "%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f",
			LSR00*lightcoeffs[0][0][2], LSG00*lightcoeffs[1][0][2], LSB00*lightcoeffs[2][0][2], //L00
			LSR*lightcoeffs[0][1][1], LSG*lightcoeffs[1][1][1], LSB*lightcoeffs[2][1][1], //L1m1
			LSR*lightcoeffs[0][1][2], LSG*lightcoeffs[1][1][2], LSB*lightcoeffs[2][1][2], //L10
			LSR*lightcoeffs[0][1][3], LSG*lightcoeffs[1][1][3], LSB*lightcoeffs[2][1][3], //L11
			LSR*lightcoeffs[0][2][0], LSG*lightcoeffs[1][2][0], LSB*lightcoeffs[2][2][0], //L2m2
			LSR*lightcoeffs[0][2][1], LSG*lightcoeffs[1][2][1], LSB*lightcoeffs[2][2][1], //L2m1
			LSR*lightcoeffs[0][2][2], LSG*lightcoeffs[1][2][2], LSB*lightcoeffs[2][2][2], //L20
			LSR*lightcoeffs[0][2][3], LSG*lightcoeffs[1][2][3], LSB*lightcoeffs[2][2][3], //L21
			LSR*lightcoeffs[0][2][4], LSG*lightcoeffs[1][2][4], LSB*lightcoeffs[2][2][4] //L22
			);
		response->write(sb);
	};

	server_interface.resource["^/SetLight$"]["POST"] = [&LSR00, &LSG00, &LSB00, &LSR, &LSG, &LSB](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request)
	{
		// Retrieve string:
		// request->content.string() is a convenience function for:
		// stringstream ss;
		// ss << request->content.rdbuf();
		// auto content=ss.str();


		std::string::size_type sz;
		//cout << "SET LIGHT REQUEST RECEIVED:" << content << endl;
		ptree pt;
		read_json(request->content, pt);
		//cout << "SET LIGHT REQUEST PARSED:" << endl;
		auto lsr00_srt=pt.get<string>("LSR00");
		LSR00 = std::stof(lsr00_srt, &sz);
		//cout << "SET LSR00 :" << LSR00 << endl;

		auto lsg00_srt = pt.get<string>("LSG00");
		LSG00 = std::stof(lsg00_srt, &sz);
		//cout << "SET LSG00 :" << LSG00 << endl;;

		auto lsb00_srt = pt.get<string>("LSB00");
		LSB00 = std::stof(lsb00_srt, &sz);
		//cout << "SET LSB00 :" << LSB00 << endl;;

		auto lsr_srt = pt.get<string>("LSR");
		LSR = std::stof(lsr_srt, &sz);
		//cout << "SET LSR :" << LSR << endl;;

		auto lsg_srt = pt.get<string>("LSG");
		LSG = std::stof(lsg_srt, &sz);
		//cout << "SET LSG :" << LSG << endl;;

		auto lsb_srt = pt.get<string>("LSB");
		LSB = std::stof(lsb_srt, &sz);
		//cout << "SET LSB :" << LSB << endl;;

		// Alternatively, use one of the convenience functions, for instance:
		// response->write(content);

		*response << "HTTP/1.1 200 OK\r\nContent-Length: " << 100 << "\r\n\r\n"
			<< (LSR00+ LSG00+ LSB00 + LSR + LSG + LSB);
	};

	server_interface.resource["^/SetImage$"]["POST"] = [](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request)
	{
		// Retrieve string:
		auto content = request->content.string();
		// request->content.string() is a convenience function for:
		// stringstream ss;
		// ss << request->content.rdbuf();
		// auto content=ss.str();

		*response << "HTTP/1.1 200 OK\r\nContent-Length: " << content.length() << "\r\n\r\n"
			<< content;
		input_image_filename = content;
		reload_image = true;

		// Alternatively, use one of the convenience functions, for instance:
		// response->write(content);
	};

	server_interface.resource["^/SaveFrame"]["POST"] = [](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request)
	{
		// Retrieve string:
		auto content = request->content.string();
		// request->content.string() is a convenience function for:

		*response << "HTTP/1.1 200 OK\r\nContent-Length: " << content.length() << "\r\n\r\n"
			<< content;
		save_frame_path = content;
		save_frame = true;
	};

	server_interface.resource["^/DrawImage$"]["POST"] = [](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request)
	{
		// Retrieve string:
		ptree pt;
		std::ostringstream ss;
		ss << request->content.rdbuf();
		std::string s = ss.str();
		std::ofstream out;
		out.open("output.jpg");
		out << request->content.rdbuf() << endl;
		out.close();
		// request->content.string() is a convenience function for:
		// stringstream ss;
		// ss << request->content.rdbuf();
		// auto content=ss.str();

		*response << "HTTP/1.1 200 OK\r\nContent-Length: " << 0 << "\r\n\r\n";
		// Alternatively, use one of the convenience functions, for instance:
		// response->write(content);
	};

	thread server_interface_thread([&server_interface]()
	{
		// Start server
		server_interface.start();
	});

	thread server_thread([&server]()
	{
		// Start server
		server.start();
	});

	/*
	thread client_thread([&server]()
	{
		// Wait for server to start so that the client can connect
		this_thread::sleep_for(chrono::seconds(1));

		// Client examples
		HttpClient client("localhost:8080");
		auto r = client.request("GET");
		for (;;)
		{
			this_thread::sleep_for(chrono::milliseconds(500));
			client.io_service->run();
			try
			{
				r.get();
				if (r)
				cout << r->content.rdbuf() << endl;
			}
			catch (const SimpleWeb::system_error &e)
			{
				cerr << "Client request error: " << e.what() << endl;
			}
			
		}
	});
	*/
	if (is_test)
	{
		input_image_filename = argv[2];
		file_spherical_harmonics();
	}
	else
		webcam_spherical_harmonics();
	
	//client_thread.join();
	server_thread.join();
	
	//future<void> SphericalHarmonics(async([]()
	//{
	//	webcam_spherical_harmonics();
	//}
	//));
	
	//for (;;)
	//{
		//printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L00 ", lightcoeffs[0][0][2], lightcoeffs[1][0][2], lightcoeffs[2][0][2]);
		//printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L1m1", lightcoeffs[0][1][1], lightcoeffs[1][1][1], lightcoeffs[2][1][1]);
		//printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L10 ", lightcoeffs[0][1][2], lightcoeffs[1][1][2], lightcoeffs[2][1][2]);
		//printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L11 ", lightcoeffs[0][1][3], lightcoeffs[1][1][3], lightcoeffs[2][1][3]);
		//printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L2m2", lightcoeffs[0][2][0], lightcoeffs[1][2][0], lightcoeffs[2][2][0]);
		//printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L2m1", lightcoeffs[0][2][1], lightcoeffs[1][2][1], lightcoeffs[2][2][1]);
		//printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L20 ", lightcoeffs[0][2][2], lightcoeffs[1][2][2], lightcoeffs[2][2][2]);
		//printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L21 ", lightcoeffs[0][2][3], lightcoeffs[1][2][3], lightcoeffs[2][2][3]);
		//printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L22 ", lightcoeffs[0][2][4], lightcoeffs[1][2][4], lightcoeffs[2][2][4]);
	//}
	
	
	//RicohTheta theta(thetaUrl);
	//theta.startSession();
	//theta.startLivePreview();
	/*
	VideoCapture cap;
	if (!cap.open(0))
		return 0;
	for (;;)
	{
		Mat frame;
		cap >> frame;
		if (frame.empty()) break;
		imshow("this is you, smile!", frame);
		if (waitKey(10) == 27) break;
	}
	
	return 0;
	*/
	
}