//hts2beta_v2.exe
//_ないやつ
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <windows.h>
#include <filesystem>
#include <iomanip>
#include <chrono>
#include <omp.h>

int do_hts_ali_one(std::string filepath);

int do_hts2beta_one(std::string file_in_data, std::string file_out_data, int arg, int flg);
int do_hts_deadpixel_one(std::string filepath, int arg);
int do_hts_fvxx_one(std::string filepath, std::string &file_out_data, int pl, int arg);
int do_f_filter_one(std::string file_in_data, std::string file_out_data, int pl, int arg);

PROCESS_INFORMATION do_prg(std::string path, std::string command, std::string output_log, HANDLE &h_out);
void process_combination_buffer(std::uintmax_t file_size[4], std::vector<std::vector<int>> &combination, int buffer);
void process_combination_memory(std::uintmax_t file_size[4], std::vector<std::vector<int>> &combination, int memory);
void process_combination(std::uintmax_t file_size[4], std::vector<int> &combination, int buffer, int memory, int &num_thread);

void console_error_out(std::string filename);

bool sort_filesize(const std::pair<int, std::uintmax_t>& left, const std::pair<int, std::uintmax_t>& right) {
	return left.second > right.second;
}
int main(int argc, char**argv) {
	if (argc != 8) {
		fprintf(stderr, "usage:prg input-folder intermediate-folder buffer[GB] memory[GB] output-folder pl flg\n");
		exit(1);
	}
	std::string file_in_data = argv[1];
	std::string file_intermadeitae_data = argv[2];
	int buffer = std::stoi(argv[3]);
	int memory = std::stoi(argv[4]);
	std::string file_out_data = argv[5];
	int pl = std::stoi(argv[6]);
	int flg = std::stoi(argv[7]);
	if (flg != 0 && flg != 1 && flg != 2 && flg != 3 && flg != 4) {
		fprintf(stderr, "flg=%d permitted value = 0 or 1 or 2 or 3 or 4\n", flg);
		fprintf(stderr, "flg=0 : normal\n");
		fprintf(stderr, "flg=1 : thin 1 vph=3 cut\n");
		fprintf(stderr, "flg=2 : thick vph=3 cut\n");
		fprintf(stderr, "flg=3 : (data size) small flg=0 large flg=2 (ECC mode)\n");
		fprintf(stderr, "flg=4 : (data size) small flg=0 large flg=1 (Shifter mode)\n");
		return 1;
	}
	bool ret;
	ret = std::filesystem::create_directories(file_intermadeitae_data);
	ret = std::filesystem::create_directories(file_out_data);

	std::filesystem::copy(file_in_data + "\\Beta_EachImagerParam.json", file_intermadeitae_data + "\\Beta_EachImagerParam.json", std::filesystem::copy_options::overwrite_existing);
	std::filesystem::copy(file_in_data + "\\Beta_EachShotParam.json", file_intermadeitae_data + "\\Beta_EachShotParam.json", std::filesystem::copy_options::overwrite_existing);
	std::filesystem::copy(file_in_data + "\\Beta_EachViewParam.json", file_intermadeitae_data + "\\Beta_EachViewParam.json", std::filesystem::copy_options::overwrite_existing);

	auto start = std::chrono::system_clock::now(); // 計測開始時間

	std::string file_name = file_in_data + "\\DATA";
	size_t DATA_size = 0;
	for (const auto & entry : std::filesystem::directory_iterator(file_name)) {
		for (const auto & entry2 : std::filesystem::directory_iterator(entry.path())) {
			DATA_size += std::filesystem::file_size(entry2.path());
		}
	}
	printf("DATA size %5.1lf[GB]\n", DATA_size*1.0 / 1000 / 1000 / 1000);
	//datasizeが大きい場合はflg=2に変更
	if (flg == 3) {
		if (DATA_size*1.0 / 1000 / 1000 / 1000 > 120) {
			flg = 2;
		}
		else {
			flg = 0;
		}
	}
	else if (flg == 4) {
		if (DATA_size*1.0 / 1000 / 1000 / 1000 > 120) {
			flg = 1;
		}
		else {
			flg = 0;
		}
	}

	int thread_num = 0;
	double beta_size = (DATA_size * 1.0 / 1000 / 1000 / 1000) * 2 / 3;
	if (memory > beta_size) {
		thread_num = 4;
		printf("thread num=%d\n", thread_num);
#pragma omp parallel sections 
		{
#pragma omp section 
			{
				if (do_hts2beta_one(file_in_data, file_intermadeitae_data, 0, flg))exit(1);
			}
#pragma omp section 
			{
				if (do_hts2beta_one(file_in_data, file_intermadeitae_data, 1, flg))exit(1);
			}
#pragma omp section 
			{
				if (do_hts2beta_one(file_in_data, file_intermadeitae_data, 2, flg))exit(1);
			}
#pragma omp section 
			{
				if (do_hts2beta_one(file_in_data, file_intermadeitae_data, 3, flg))exit(1);
			}
		}
	}
	else if (memory * 2 > beta_size) {
		thread_num = 2;
		printf("thread num=%d\n", thread_num);
#pragma omp parallel sections 
		{
#pragma omp section 
			{
				if (do_hts2beta_one(file_in_data, file_intermadeitae_data, 0, flg))exit(1);
			}
#pragma omp section 
			{
				if (do_hts2beta_one(file_in_data, file_intermadeitae_data, 2, flg))exit(1);
			}
		}
#pragma omp parallel sections 
		{
#pragma omp section 
			{
				if (do_hts2beta_one(file_in_data, file_intermadeitae_data, 1, flg))exit(1);
			}
#pragma omp section 
			{
				if (do_hts2beta_one(file_in_data, file_intermadeitae_data, 3, flg))exit(1);
			}
		}
	}
	else {
		thread_num = 1;
		printf("thread num=%d\n", thread_num);
		if (do_hts2beta_one(file_in_data, file_intermadeitae_data, 0, flg))exit(1);
		if (do_hts2beta_one(file_in_data, file_intermadeitae_data, 1, flg))exit(1);
		if (do_hts2beta_one(file_in_data, file_intermadeitae_data, 2, flg))exit(1);
		if (do_hts2beta_one(file_in_data, file_intermadeitae_data, 3, flg))exit(1);
	}

	if (do_hts_ali_one(file_intermadeitae_data))exit(1);

	//Byte単位で出る
	std::uintmax_t file_size[4];
	file_size[0] = std::filesystem::file_size(file_intermadeitae_data + "\\beta_thick_0.dat");
	file_size[1] = std::filesystem::file_size(file_intermadeitae_data + "\\beta_thick_1.dat");
	file_size[2] = std::filesystem::file_size(file_intermadeitae_data + "\\beta_thin_0.dat");
	file_size[3] = std::filesystem::file_size(file_intermadeitae_data + "\\beta_thin_1.dat");

	printf("beta thick 0 : %.1lf[GB]\n", file_size[0] * 1. / (1000 * 1000 * 1000));
	printf("beta thick 1 : %.1lf[GB]\n", file_size[1] * 1. / (1000 * 1000 * 1000));
	printf("beta thin  0 : %.1lf[GB]\n", file_size[2] * 1. / (1000 * 1000 * 1000));
	printf("beta thin  1 : %.1lf[GB]\n", file_size[3] * 1. / (1000 * 1000 * 1000));
	std::vector<int> combination;

	process_combination(file_size, combination, buffer, memory, thread_num);
	printf("number of thread = %d\n", thread_num);
	//process_combination_buffer(file_size, combination, buffer);
	//process_combination_memory(file_size, combination, memory);
#pragma omp parallel for num_threads(thread_num) schedule(dynamic,1)
	for (int i = 0; i < combination.size(); i++) {
		int  arg = combination[i];
#pragma omp critical
		printf("process %d, order %d start\n", arg, i);
		if (do_hts_deadpixel_one(file_intermadeitae_data, arg))exit(1);
		std::string file_out_tmp;
		if (do_hts_fvxx_one(file_intermadeitae_data, file_out_tmp, pl, arg))exit(1);
		if (do_f_filter_one(file_out_tmp, file_out_data, pl, arg))exit(1);
		//if (do_f_filter_one2(file_out_tmp, file_out_data, pl, arg))exit(1);
#pragma omp critical
		printf("process %d, order %d end\n", arg, i);
	}

	auto end = std::chrono::system_clock::now(); // 計測開始時間
	double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start).count(); //処理に要した時間をミリ秒に変換

	printf("hts2fvxx time:%.0lf[s]\n", elapsed);
	std::filesystem::copy(file_intermadeitae_data + "\\Beta_EachImagerParam.json", file_out_data + "\\Beta_EachImagerParam.json", std::filesystem::copy_options::overwrite_existing);
	std::filesystem::copy(file_intermadeitae_data + "\\Beta_EachShotParam.json", file_out_data + "\\Beta_EachShotParam.json", std::filesystem::copy_options::overwrite_existing);
	std::filesystem::copy(file_intermadeitae_data + "\\Beta_EachViewParam.json", file_out_data + "\\Beta_EachViewParam.json", std::filesystem::copy_options::overwrite_existing);
	std::filesystem::copy(file_intermadeitae_data + "\\ali.json", file_out_data + "\\ali.json", std::filesystem::copy_options::overwrite_existing);
	std::filesystem::copy(file_intermadeitae_data + "\\ali.pdf", file_out_data + "\\ali.pdf", std::filesystem::copy_options::overwrite_existing);
	std::filesystem::copy(file_intermadeitae_data + "\\deadpixel.json", file_out_data + "\\deadpixel.json", std::filesystem::copy_options::overwrite_existing);

	//printf("delete start [%s]\n", file_intermadeitae_data.c_str());
	//std::filesystem::remove_all(file_intermadeitae_data);
	//printf("delete start [%s]\n", file_in_data.c_str());
	//std::filesystem::remove_all(file_in_data);
	//printf("delete end   [%s]\n", file_in_data.c_str());
	return 0;
}

int do_hts_ali_one(std::string filepath) {
	std::stringstream command[2];

	//command[0] << "hts_beta_ali_v2.exe " << filepath << "\\beta_thick_0.dat ali.json --th_angle2 0.02 --th_pos2 0.015";
	//command[0] << "hts_beta_ali_v2.exe " << filepath << "\\beta_thick_0.dat ali.json --th_angle2 0.02 --th_pos2 0.020";//2022/04/04
	command[0] << "hts_beta_ali_v2.exe " << filepath << "\\beta_thick_0.dat ali.json --th_angle2 0.04 --th_pos2 0.030";//2023/03/14
	command[1] << "root_macro ali.json ali.pdf plot_hts_beta_ali";

	PROCESS_INFORMATION p;
	std::stringstream file_out_log;
	HANDLE h_out;
	//hts2beta start
	auto start = std::chrono::system_clock::now(); // 計測開始時間
	file_out_log << filepath << "\\log0.txt";
	p = do_prg(filepath, command[0].str(), file_out_log.str(), h_out);

	DWORD rc;
	WaitForSingleObject(p.hProcess, INFINITE);
	GetExitCodeProcess(p.hProcess, &rc);
	if (rc != 0) {
		TerminateProcess(p.hProcess, 0);
		CloseHandle(p.hProcess);
		CloseHandle(h_out);
		console_error_out(file_out_log.str());
		return 1;
	}

	p = do_prg(filepath, command[1].str(), file_out_log.str(), h_out);
	WaitForSingleObject(p.hProcess, INFINITE);
	GetExitCodeProcess(p.hProcess, &rc);
	if (rc != 0) {
		TerminateProcess(p.hProcess, 0);
		CloseHandle(p.hProcess);
		CloseHandle(h_out);
		console_error_out(file_out_log.str());
		return 1;
	}

	auto end = std::chrono::system_clock::now(); // 計測開始時間
	double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start).count(); //処理に要した時間をミリ秒に変換

	printf("hts_ali time:%.0lf[s]\n", elapsed);
	return 0;
}

int do_hts2beta_one(std::string file_in_data, std::string file_out_data, int arg, int flg) {
	//すでにbetaがあればreturn 0
	std::stringstream beta_file_name;
	if (arg == 0)	beta_file_name << file_out_data << "\\beta_thick_0.dat";
	else if (arg == 1)	beta_file_name << file_out_data << "\\beta_thick_1.dat";
	else if (arg == 2)	beta_file_name << file_out_data << "\\beta_thin_0.dat";
	else if (arg == 3)	beta_file_name << file_out_data << "\\beta_thin_1.dat";
	if (std::filesystem::exists(beta_file_name.str()))return 0;

	std::stringstream command;
	//if (arg == 0)	command << "_hts2beta_v2.exe " << file_out_data << "\\beta_thick_0.dat -alternate0 --outputformat 1 -removeoutofrange -ignorelackfile";
	//else if (arg == 1)command << "_hts2beta_v2.exe " << file_out_data << "\\beta_thick_1.dat -alternate1 --outputformat 1 -removeoutofrange -ignorelackfile";
	//else if (arg == 2)	command << "_hts2beta_v2.exe " << file_out_data << "\\beta_thin_0.dat -thin16_base --outputformat 1 -removeoutofrange -ignorelackfile";
	//else if (arg == 3)	command << "_hts2beta_v2.exe " << file_out_data << "\\beta_thin_1.dat -thin16_outer --outputformat 1 -removeoutofrange -ignorelackfile";
	if (arg == 0 && flg == 2)	command << "hts2beta_v2.exe " << file_out_data << "\\beta_thick_0.dat -alternate0 --outputformat 1 -removeoutofrange  -ignorelackfile --volcut {\\\"Space\\\":[1000],\\\"Th\\\":[3]}";
	else if (arg == 0)	command << "hts2beta_v2.exe " << file_out_data << "\\beta_thick_0.dat -alternate0 --outputformat 1 -removeoutofrange  -ignorelackfile";
	else if (arg == 1 && flg == 2)command << "hts2beta_v2.exe " << file_out_data << "\\beta_thick_1.dat -alternate1 --outputformat 1 -removeoutofrange  -ignorelackfile --volcut {\\\"Space\\\":[1000],\\\"Th\\\":[3]}";
	else if (arg == 1)command << "hts2beta_v2.exe " << file_out_data << "\\beta_thick_1.dat -alternate1 --outputformat 1 -removeoutofrange  -ignorelackfile";
	else if (arg == 2)	command << "hts2beta_v2.exe " << file_out_data << "\\beta_thin_0.dat -thin16_base --outputformat 1 -removeoutofrange  -ignorelackfile";
	//else if (arg == 2 && flg == 0)	command << "_hts2beta_v2.exe " << file_out_data << "\\beta_thin_0.dat -thin16_base --outputformat 1 -removeoutofrange";
	//else if (arg == 2 && flg == 1)	command << "_hts2beta_v2.exe " << file_out_data << "\\beta_thin_0.dat -thin16_base --outputformat 1 -removeoutofrange --volcut {\\\"Space\\\":[1000],\\\"Th\\\":[3]}";
	else if (arg == 3 && flg == 1)	command << "hts2beta_v2.exe " << file_out_data << "\\beta_thin_1.dat -thin16_outer --outputformat 1 -removeoutofrange  -ignorelackfile --volcut {\\\"Space\\\":[1000],\\\"Th\\\":[3]}";
	else if (arg == 3)	command << "hts2beta_v2.exe " << file_out_data << "\\beta_thin_1.dat -thin16_outer --outputformat 1 -removeoutofrange  -ignorelackfile";
	else {
		fprintf(stderr, "arg error\n");
		exit(1);
	}

	PROCESS_INFORMATION p;
	std::stringstream file_out_log;
	HANDLE h_out;
	//hts2beta start
	auto hts2beta_start = std::chrono::system_clock::now(); // 計測開始時間
	file_out_log << file_out_data << "\\log" << arg << ".txt";
	p = do_prg(file_in_data, command.str(), file_out_log.str(), h_out);

	DWORD rc;
	WaitForSingleObject(p.hProcess, INFINITE);
	GetExitCodeProcess(p.hProcess, &rc);
	if (rc != 0) {
		TerminateProcess(p.hProcess, 0);
		CloseHandle(p.hProcess);
		CloseHandle(h_out);
		console_error_out(file_out_log.str());
		return 1;

	}
	CloseHandle(p.hProcess);
	CloseHandle(h_out);

	auto hts2beta_end = std::chrono::system_clock::now(); // 計測開始時間
	double elapsed = std::chrono::duration_cast<std::chrono::seconds>(hts2beta_end - hts2beta_start).count(); //処理に要した時間をミリ秒に変換

	printf("hts2beta[%d] time:%.0lf[s]\n", arg, elapsed);
	return 0;
}
int do_hts_deadpixel_one(std::string filepath, int arg) {
	std::stringstream command;
	if (arg == 0)	command << "hts_beta_deadpixel_v2.exe " << filepath << "\\beta_thick_0.dat deadpixel.json";
	else if (arg == 1)command << "hts_beta_deadpixel_v2.exe " << filepath << "\\beta_thick_1.dat deadpixel.json";
	else if (arg == 2)	command << "hts_beta_deadpixel_v2.exe " << filepath << "\\beta_thin_0.dat deadpixel.json";
	else if (arg == 3)	command << "hts_beta_deadpixel_v2.exe " << filepath << "\\beta_thin_1.dat deadpixel.json";
	else {
		fprintf(stderr, "arg error\n");
		exit(1);
	}


	PROCESS_INFORMATION p;
	std::stringstream file_out_log;
	HANDLE h_out;
	//hts2beta start
	auto hts2beta_start = std::chrono::system_clock::now(); // 計測開始時間
	file_out_log << filepath << "\\log" << arg << ".txt";
	p = do_prg(filepath, command.str(), file_out_log.str(), h_out);

	DWORD rc;
	WaitForSingleObject(p.hProcess, INFINITE);
	GetExitCodeProcess(p.hProcess, &rc);
	int flg = 0;
	if (rc != 0) {
		TerminateProcess(p.hProcess, 0);
		CloseHandle(p.hProcess);
		CloseHandle(h_out);
		console_error_out(file_out_log.str());
		return 1;

	}
	CloseHandle(p.hProcess);
	CloseHandle(h_out);

	auto hts2beta_end = std::chrono::system_clock::now(); // 計測開始時間
	double elapsed = std::chrono::duration_cast<std::chrono::seconds>(hts2beta_end - hts2beta_start).count(); //処理に要した時間をミリ秒に変換

	printf("hts_deadpixel[%d] time:%.0lf[s]\n", arg, elapsed);
	return 0;
}
int do_hts_fvxx_one(std::string filepath, std::string &file_out_data, int pl, int arg) {
	std::stringstream out_file;
	out_file << filepath << "\\out" << arg;
	bool ret = std::filesystem::create_directories(out_file.str());
	file_out_data = out_file.str();

	std::stringstream command;
	std::string key;
	if (arg == 0)key = "thick_0";
	else if (arg == 1)key = "thick_1";
	else if (arg == 2)key = "thin_0";
	else if (arg == 3)key = "thin_1";
	else {
		fprintf(stderr, "arg error\n");
		exit(1);
	}
	command << "hts_beta_fvxx_v2.exe " << filepath << "\\beta_" << key << ".dat " << pl << " --directory " << out_file.str() << "\\ -numvola2pxpy ";

	PROCESS_INFORMATION p;
	std::stringstream file_out_log;
	HANDLE h_out;
	//hts2beta start
	auto hts2beta_start = std::chrono::system_clock::now(); // 計測開始時間
	file_out_log << filepath << "\\log" << arg << ".txt";
	p = do_prg(filepath, command.str(), file_out_log.str(), h_out);

	DWORD rc;
	WaitForSingleObject(p.hProcess, INFINITE);
	GetExitCodeProcess(p.hProcess, &rc);
	int flg = 0;
	if (rc != 0) {
		TerminateProcess(p.hProcess, 0);
		CloseHandle(p.hProcess);
		CloseHandle(h_out);
		console_error_out(file_out_log.str());
		return 1;

	}
	CloseHandle(p.hProcess);
	CloseHandle(h_out);

	auto hts2beta_end = std::chrono::system_clock::now(); // 計測開始時間
	double elapsed = std::chrono::duration_cast<std::chrono::seconds>(hts2beta_end - hts2beta_start).count(); //処理に要した時間をミリ秒に変換

	printf("hts_fvxx[%d] time:%.0lf[s]\n", arg, elapsed);
	//成功していたらbeta 削除
	//std::filesystem::remove(filepath + "\\beta_" + key + ".dat ");

	return 0;
}
int do_f_filter_one(std::string file_in_data, std::string file_out_data, int pl, int arg) {
	std::stringstream command[2];
	std::string key;
	if (arg == 0)key = "_thick_0";
	else if (arg == 1)key = "_thick_1";
	else if (arg == 2)key = "_thin_0";
	else if (arg == 3)key = "_thin_1";
	else {
		fprintf(stderr, "arg error\n");
		exit(1);
	}

	command[0] << "f_filter " << pl * 10 + 1 << " "
		<< file_in_data << "\\f" << std::setw(4) << std::setfill('0') << pl * 10 + 1 << ".vxx --o "
		<< file_out_data << "\\f" << std::setw(4) << std::setfill('0') << pl * 10 + 1 << key << ".vxx"
		<< " --ghost-rl 5 5 0.01 0.01 0.05 5 --view 5000 500";

	command[1] << "f_filter " << pl * 10 + 2 << " "
		<< file_in_data << "\\f" << std::setw(4) << std::setfill('0') << pl * 10 + 2 << ".vxx --o "
		<< file_out_data << "\\f" << std::setw(4) << std::setfill('0') << pl * 10 + 2 << key << ".vxx"
		<< " --ghost-rl 5 5 0.01 0.01 0.05 5 --view 5000 500";

	PROCESS_INFORMATION p[2];
	std::stringstream file_out_log[2];
	HANDLE h_out[2];
	//hts2beta start
	auto start = std::chrono::system_clock::now(); // 計測開始時間
	file_out_log[0] << file_in_data << "\\log" << arg << "_0.txt";
	file_out_log[1] << file_in_data << "\\log" << arg << "_1.txt";
	p[0] = do_prg(file_in_data, command[0].str(), file_out_log[0].str(), h_out[0]);
	Sleep(10 * 1000);
	p[1] = do_prg(file_in_data, command[1].str(), file_out_log[1].str(), h_out[1]);

	DWORD rc;
	for (int i = 0; i < 2; i++) {
		WaitForSingleObject(p[i].hProcess, INFINITE);
		GetExitCodeProcess(p[i].hProcess, &rc);

		if (rc != 0) {
			for (int j = i; j < 2; j++) {
				TerminateProcess(p[j].hProcess, 0);
				CloseHandle(p[j].hProcess);
				CloseHandle(h_out[j]);
			}
			console_error_out(file_out_log[i].str());
			printf("error arg%d\n", arg);
			return 1;
		}
		CloseHandle(p[i].hProcess);
		CloseHandle(h_out[i]);
	}
	auto end = std::chrono::system_clock::now(); // 計測開始時間
	double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start).count(); //処理に要した時間をミリ秒に変換

	printf("hts_f_filter[%d] time:%.0lf[s]\n", arg, elapsed);
	//成功していたらfvxx 削除
	//printf("start %s\n", file_in_data.c_str());
	//std::filesystem::remove_all(file_in_data);
	//printf("end %s\n", file_in_data.c_str());

	//std::stringstream del_file[2];
	//bool res;
	//del_file[0] << file_in_data << "\\f" << std::setw(4) << std::setfill('0') << pl * 10 + 1 << ".vxx";
	//del_file[1] << file_in_data << "\\f" << std::setw(4) << std::setfill('0') << pl * 10 + 2 << ".vxx";

	//printf("start %s\n", del_file[0].str().c_str());
	//res=std::filesystem::remove(del_file[0].str());
	//printf("end %s result %d\n", del_file[0].str().c_str(),res);

	//printf("start %s\n", del_file[1].str().c_str());
	//res = std::filesystem::remove(del_file[1].str());
	//printf("end %s result %d\n", del_file[1].str().c_str(), res);

	return 0;

}
PROCESS_INFORMATION do_prg(std::string path, std::string command, std::string output_log, HANDLE &h_out) {
	SECURITY_ATTRIBUTES sec_attr;
	ZeroMemory(&sec_attr, sizeof(sec_attr));
	sec_attr.bInheritHandle = TRUE;

	//HANDLE h_out = CreateFile(TEXT("CONOUT$"), GENERIC_WRITE, FILE_SHARE_WRITE, NULL, OPEN_EXISTING, 0, 0);
	//HANDLE h_err = CreateFile(TEXT("CONOUT$"), GENERIC_WRITE, FILE_SHARE_WRITE, NULL, OPEN_EXISTING, 0, 0);
	h_out = CreateFile(output_log.c_str(), FILE_APPEND_DATA, FILE_SHARE_WRITE, &sec_attr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);


	PROCESS_INFORMATION p;
	STARTUPINFO s;

	ZeroMemory(&s, sizeof(s));
	s.cb = sizeof(s);
	s.dwFlags = STARTF_USESTDHANDLES;
	s.hStdOutput = h_out;
	s.hStdError = h_out;

	LPSTR str = const_cast<char *>(command.c_str());
	//printf("%s\n", command.c_str());
	int ret = CreateProcess(
		NULL, // 実行可能モジュールの名
		str, // コマンドラインの文字列
		&sec_attr, // セキュリティ記述子
		&sec_attr,// セキュリティ記述子
		TRUE, // ハンドルの継承オプション
		NORMAL_PRIORITY_CLASS, // 作成のフラグ
		NULL,// 新しい環境ブロック
		path.c_str(), // カレントディレクトリの名前
		&s, // スタートアップ情報
		&p // プロセス情報
	);
	if (!ret)
	{
		printf("miss %s\n", command.c_str());
		exit(1);
	}
	return p;

}

void process_combination(std::uintmax_t file_size[4], std::vector<int> &combination, int buffer, int memory, int &num_thread) {
	//ディスクを[buffer]GBと仮定
	//fvxx=3.24*beta
	//必要buffer=4.24*beta[GB]
	std::vector<std::pair<int, std::uintmax_t>>file_id;

	for (int i = 0; i < 4; i++) {
		file_id.push_back(std::make_pair(i, file_size[i]));
	}
	sort(file_id.begin(), file_id.end(), sort_filesize);
	for (auto itr = file_id.begin(); itr != file_id.end(); itr++) {
		printf("id:%d filesize:%.1lf[GB]\n", itr->first, itr->second*1. / (1000 * 1000 * 1000));
		combination.push_back(itr->first);
	}

	if (buffer > (file_id[0].second + file_id[1].second + file_id[2].second + file_id[3].second)*4.24 / (1000 * 1000 * 1000)) {
		if (memory > (file_id[0].second + file_id[1].second + file_id[2].second + file_id[3].second) / (1000 * 1000 * 1000)) {
			num_thread = 4;
		}
		else if (memory > (file_id[0].second + file_id[1].second + file_id[2].second) / (1000 * 1000 * 1000)) {
			num_thread = 3;
		}
		else if (memory > (file_id[0].second + file_id[1].second) / (1000 * 1000 * 1000)) {
			num_thread = 2;
		}
		else {
			num_thread = 1;
		}
	}
	else if (buffer > (file_id[0].second + file_id[1].second + file_id[2].second)*4.24 / (1000 * 1000 * 1000)) {
		if (memory > (file_id[0].second + file_id[1].second + file_id[2].second) / (1000 * 1000 * 1000)) {
			num_thread = 3;
		}
		else if (memory > (file_id[0].second + file_id[1].second) / (1000 * 1000 * 1000)) {
			num_thread = 2;
		}
		else {
			num_thread = 1;
		}
	}
	else if (buffer > (file_id[0].second + file_id[1].second)*4.24 / (1000 * 1000 * 1000)) {
		if (memory > (file_id[0].second + file_id[1].second) / (1000 * 1000 * 1000)) {
			num_thread = 2;
		}
		else {
			num_thread = 1;
		}
	}
	else {
		num_thread = 1;
	}
}
void process_combination_buffer(std::uintmax_t file_size[4], std::vector<std::vector<int>> &combination, int buffer) {
	//ディスクを[buffer]GBと仮定
	//fvxx=3.24*beta
	//必要buffer=4.24*beta[GB]
	int remain = buffer;
	int file_size_GB[4];
	std::multimap<int, int> sort;

	for (int i = 0; i < 4; i++) {
		file_size_GB[i] = 1 + file_size[i] / (1000 * 1000 * 1000);
		sort.insert(std::make_pair(file_size_GB[i], i));
		remain -= file_size_GB[i];
	}
	int num = 0;
	while (num < 4) {
		//printf("num %d\n", num);
		std::vector<int> procces_pair;
		for (auto itr = sort.begin(); itr != sort.end(); itr++) {
			if (itr == sort.begin()) {
				procces_pair.push_back(itr->second);
				remain -= itr->first * 3;
				num++;
			}
			else {
				if (remain > itr->first * 3) {
					procces_pair.push_back(itr->second);
					remain -= itr->first * 3;
					num++;
				}
			}
		}
		combination.push_back(procces_pair);
		std::multimap<int, int> sort_tmp = sort;
		sort.clear();
		remain = buffer;
		bool flg = true;
		for (auto itr = sort_tmp.begin(); itr != sort_tmp.end(); itr++) {
			flg = true;
			for (auto itr2 = procces_pair.begin(); itr2 != procces_pair.end(); itr2++) {
				if (*itr2 == itr->second)flg = false;
			}
			if (flg) {
				sort.insert(*itr);
				remain -= itr->first;
			}
		}
	}
}
void process_combination_memory(std::uintmax_t file_size[4], std::vector<std::vector<int>> &combination, int memory) {

	//ディスクを[buffer]GBと仮定
	//fvxx=2*beta
	std::vector<std::vector<int>> ret;
	int remain = memory;
	double file_size_GB[4];
	std::map<int, double> beta_size;

	for (int i = 0; i < 4; i++) {
		file_size_GB[i] = file_size[i] * 1.0 / (1000 * 1000 * 1000);
		beta_size.insert(std::make_pair(i, file_size_GB[i]));
	}
	for (int i = 0; i < combination.size(); i++) {
		int now = 0;
		int cnt = 0;
		while (now < combination[i].size()) {
			std::vector<int>procces_pair;
			cnt = 0;
			remain = memory;
			for (int j = 0; j + now < combination[i].size(); j++) {
				if (cnt == 0) {
					remain -= int(beta_size[combination[i][now + j]]) + 1;
					procces_pair.push_back(combination[i][now + j]);
					cnt++;
				}
				else {
					if (remain < beta_size[combination[i][now + j]])break;
					else {
						remain -= int(beta_size[combination[i][now + j]]) + 1;
						procces_pair.push_back(combination[i][now + j]);
						cnt++;
					}
				}
			}
			ret.push_back(procces_pair);
			now += cnt;
		}
	}

	combination = ret;
	//process listの出力
	for (int i = 0; i < combination.size(); i++) {
		printf("process %d\n", i);
		for (int j = 0; j < combination[i].size(); j++) {
			printf("id:%d beta_size:%.1lf[GB]\n", combination[i][j], beta_size[combination[i][j]]);
		}
	}
}
void console_error_out(std::string filename) {
	printf("%s\n", filename.c_str());
	std::ifstream ifs(filename);
	std::string str;

	if (ifs.fail()) {
		std::cerr << "Failed to open file." << std::endl;
	}

	std::cout << ":::::::::::::::::::Error message::::::::::::::" << std::endl;

	HANDLE hStdout;
	WORD wAttributes;
	CONSOLE_SCREEN_BUFFER_INFO csbi;//構造体です

	hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
	GetConsoleScreenBufferInfo(hStdout, &csbi);
	wAttributes = FOREGROUND_RED;
	SetConsoleTextAttribute(hStdout, wAttributes);

	while (getline(ifs, str)) {
		std::cout << str << std::endl;
	}
	wAttributes = FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE;
	SetConsoleTextAttribute(hStdout, wAttributes);
	std::cout << "::::::::::::::::::::::::::::::::::::::::::" << std::endl;

}
void do_ADAPT(std::string adapt_path, std::string in_mfile) {
	PROCESS_INFORMATION p;
	STARTUPINFO s;
	ZeroMemory(&s, sizeof(s));
	s.cb = sizeof(s);

	std::string command = adapt_path + " " + in_mfile;
	LPSTR str = const_cast<char *>(command.c_str());
	int ret = CreateProcess(
		NULL, // 実行可能モジュールの名
		str, // コマンドラインの文字列
		NULL, // セキュリティ記述子
		NULL,// セキュリティ記述子
		FALSE, // ハンドルの継承オプション
		NULL, // 作成のフラグ
		NULL,// 新しい環境ブロック
		NULL, // カレントディレクトリの名前
		&s, // スタートアップ情報
		&p // プロセス情報
	);
	if (!ret)
	{
		printf("miss %s\n", command.c_str());
	}
	else
	{
		CloseHandle(p.hThread);

		//メモ帳が終了するまで待つ
		WaitForSingleObject(p.hProcess, INFINITE);
		CloseHandle(p.hProcess);
	}

}
