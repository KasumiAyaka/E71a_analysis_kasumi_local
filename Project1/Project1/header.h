#pragma once
#include <fstream>
#include <iostream>
#include <ios>
#include <iomanip> 
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <map>

#include <list>
#include <set>
#include <unordered_map>
#include <utility>
#include <random>
#include <cmath>
#include <limits>
#include <chrono>
#include <cassert>


namespace {
	std::vector<std::string> StringSplit(std::string str) {
		std::stringstream ss{ str };
		std::vector<std::string> v;
		std::string buf;
		while (std::getline(ss, buf, ' ')) {
			if (buf != "") {
				v.push_back(buf);
			}
		}
		return v;
	}
	std::vector<std::string> StringSplit_with_tab(std::string str) {
		std::vector<std::string> v;

		std::vector<std::string> str_v = StringSplit(str);
		for (int i = 0; i < str_v.size(); i++) {
			std::stringstream ss{ str_v[i] };
			std::string buf;
			while (std::getline(ss, buf, '\t')) {
				if (buf != "") {
					v.push_back(buf);
				}
			}
		}
		return v;
	}
}

//Mfile txt
namespace mfile0 {
	struct M_Header {
		std::string head[3];
		int num_all_plate;
		std::vector<int> all_pos;
	};
	struct M_Base {
		int64_t  rawid, group_id;
		int pos, ph, flg_i[4];
		double x, y, z, flg_d[2];
		float ax, ay;
	};
	struct M_Chain {
		int64_t chain_id;
		int  nseg, pos0, pos1;
		std::vector<M_Base> basetracks;
	};
	struct Mfile {
		M_Header header;
		std::vector<M_Chain> chains;
	};

	void read_mfile(std::string file, Mfile& mfile);
	void read_mfile(std::string file_path, Mfile& mfile, int nseg_thr);
	void write_mfile(std::string filename, const Mfile& mfile, int output = 1);
	void write_mfile_header(std::ofstream& ofs, const M_Header& header, int output = 1);
	void write_mfile_chain(std::ofstream& ofs, const M_Chain& chains);
	double chain_ax(M_Chain chain);
	double chain_ay(M_Chain chain);
	double angle_diff_dev_rad(M_Chain chain, double ax, double ay);
	double angle_diff_dev_lat(M_Chain chain, double ax, double ay);
	void apply_vph_correction(mfile0::Mfile& m, std::string filename_corr);
	void write_mfile_chain_IVE(std::ofstream& ofs, const M_Chain& chains);
	double chain_fit_dist(M_Chain chain, double slope[3], double intercept[3], double& dist);

	struct M_Chain_inf {
		int64_t chainID;
		int  nseg, pos0, pos1, groupID;
		double ax, ay, radial_deviation, lateral_deviation, d_radial_deviation, d_lateral_deviation, vph_ratio, vph_slope, vph_slope_acc;
		std::pair<double, double> x_up, y_up, x_down, y_down;
		bool outflg_up, outflg_down;
	};
	M_Chain_inf chain2inf(mfile0::M_Chain c);
	void chain_inf_flg(M_Chain_inf& inf, mfile0::M_Chain c, int pl0, int pl1, double range[4], std::map<int, double> z);
	void write_chaininf(std::string filename, std::vector<M_Chain_inf> chain_inf);

	void set_header(int pl0, int pl1, mfile0::Mfile& m);
}
//Mfile bin
namespace mfile1 {
	struct MFileHeader
	{
		uint64_t filetype;
		uint64_t filesize;
		uint64_t reserved;
		uint32_t offset1;
		uint32_t offset2;
		uint8_t name[64];
		uint8_t object[256];
	};

	struct MFileInfoHeader
	{
		uint16_t classsize1;
		uint16_t classsize2;
		uint32_t reserved1;
		uint64_t Nchain;
		uint64_t Nbasetrack;
		uint64_t reserved2;
	};

	class MFileBase {
	public:
		int pos, group_id;
		uint64_t rawid;
		int ph;
		float ax, ay;
		double x, y, z;
	};
	class MFileBase1 : public MFileBase
	{
	public:
		int16_t tmp[4]; // 計算で使うが、ファイルに書く必要のない一時パラメータ
	};


	class MFileChain {
	public:
		uint64_t chain_id;
		int nseg, pos0, pos1;
	};

	class MFileChain1 : public MFileChain
	{
	public:
		int16_t tmp[4]; // 計算で使うが、ファイルに書く必要のない一時パラメータ
	};

	class MFile {

	public:
		MFileHeader header;
		MFileInfoHeader info_header;
		std::vector<MFileChain1> chains;
		std::vector<std::vector<MFileBase1>> all_basetracks;
	public:
		//i番目のchain
		double chain_ax(int i);
		double chain_ay(int i);
		double angle_diff_dev_rad(int i, double ax, double ay);
		double angle_diff_dev_lat(int i, double ax, double ay);

	};
	class MFile_minimum {

	public:
		MFileHeader header;
		MFileInfoHeader info_header;
		std::vector<MFileChain> chains;
		std::vector<std::vector<MFileBase>> all_basetracks;
	public:
		//i番目のchain
		double chain_ax(int i);
		double chain_ay(int i);
		double angle_diff_dev_rad(int i, double ax, double ay);
		double angle_diff_dev_lat(int i, double ax, double ay);

	};
	void converter(const MFile& old, mfile0::Mfile& mfile);
	void converter(const MFile_minimum& old, mfile0::Mfile& mfile);
	void converter(const mfile0::Mfile& old, MFile& mfile);
	void read_mfile(std::string filepath, MFile& mfile);
	void read_mfile(std::string filepath, MFile_minimum& mfile);
	void read_mfile_txt(std::string filepath, MFile& mfile);

	void write_mfile(std::string filepath, MFile& mfile);
	void write_mfile(std::string filepath, MFile_minimum& mfile);

	void read_mfile_extension(std::string filename, mfile0::Mfile& m);
	void write_mfile_extension(std::string filename, mfile0::Mfile& m);

	double chain_ax(std::vector<MFileBase>& b);
	double chain_ay(std::vector<MFileBase>& b);
	double angle_diff_dev_rad(std::vector<MFileBase>& b);
	double angle_diff_dev_lat(std::vector<MFileBase>& b);
	double chain_vph(std::vector<MFileBase>& b);
	double chain_ph(std::vector<MFileBase>& b);
	double angle_diff_mom_iron(std::vector<MFileBase>& b, int num_thr = 5);

}