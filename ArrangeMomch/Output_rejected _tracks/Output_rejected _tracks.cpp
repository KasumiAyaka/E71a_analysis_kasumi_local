// 2026/11/06
// Output_rejected_tracks
// kasumi
//
// Output information of partner tracks (muon excluded) with reject flags
// collected from multiple input lists.
//
// output format :
//   groupid chainid vertex-material pid ax ay x y z vph nseg npl md oa MCS-momentum Range-momentum reject-flg
//
// reject-flg :
//   7 : reject by poor convergence from the vertex point (track not listed anywhere
//       in the param file; only assigned when "accept" is specified in param.txt)
//   6 : reject by mip nseg-MCS cut (nseg < 26 + 4.5*ln(MCS-momentum/500); applied
//       only to mip tracks, i.e. particle_flg != 2212 and != 13)
//   5 : reject by 10 seg cut (nseg < threshold given in param.txt; applied only to
//       mip tracks, i.e. particle_flg != 2212 and != 13)
//   4 : reject by range-mom cut
//   3 : reject by manual check
//   2 : reject by VE
//   1 : reject by viewer check
//   0 : remain / accepted (see "accept" in param.txt)
// (if more than one of the above applies to a track, the smallest reject-flg wins)
//
// md(minimum distance), oa(opening angle) are calculated between
// muon(particle_flg==13 or chainid==0) and the partner track.

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}

struct ListPaths {
	std::string viewer_check, VE, manual_check, range_mom_cut, accept;
	int nseg_cut_thr = -1;// nseg < nseg_cut_thr -> reject_flg = 5 (disabled if < 0)
};

double minimum_distance_fixed(matrix_3D::vector_3D pos0, matrix_3D::vector_3D pos1, matrix_3D::vector_3D dir0, matrix_3D::vector_3D dir1, double z_range[2], double extra[2], double refz) {
	double extra0_distance, extra1_distance, delta;
	matrix_3D::vector_3D pos;
	pos.x = pos1.x - pos0.x;
	pos.y = pos1.y - pos0.y;
	pos.z = pos1.z - pos0.z;
	// almost parallel case
	if (opening_angle(dir0, dir1) < 0.0001) {
		extra0_distance = (pos1.z + pos0.z) / 2 - pos0.z;
		extra1_distance = (pos1.z + pos0.z) / 2 - pos1.z;
	}
	else {
		delta = dot(dir0, dir0) * dot(dir1, dir1) - pow(dot(dir0, dir1), 2.);
		extra0_distance = (+1 * dot(pos, dir0) * dot(dir1, dir1) - dot(dir0, dir1) * dot(pos, dir1)) / delta;
		extra1_distance = (-1 * dot(pos, dir1) * dot(dir0, dir0) + dot(dir0, dir1) * dot(pos, dir0)) / delta;
	}
	if (z_range[0] > z_range[1]) {
		double tmp_d = z_range[0];
		z_range[0] = z_range[1];
		z_range[1] = tmp_d;
	}

	matrix_3D::vector_3D extra0 = addition(pos0, const_multiple(dir0, extra0_distance));
	matrix_3D::vector_3D extra1 = addition(pos1, const_multiple(dir1, extra1_distance));

	if (extra0.z < refz + z_range[0] || extra1.z < refz + z_range[0]) {
		extra0_distance = refz - pos0.z + z_range[0];
		extra1_distance = refz - pos1.z + z_range[0];
	}
	else if (extra0.z > refz + z_range[1] || extra1.z > refz + z_range[1]) {
		extra0_distance = refz - pos0.z + z_range[1];
		extra1_distance = refz - pos1.z + z_range[1];
	}

	extra[0] = extra0_distance;
	extra[1] = extra1_distance;
	extra0 = addition(pos0, const_multiple(dir0, extra0_distance));
	extra1 = addition(pos1, const_multiple(dir1, extra1_distance));

	return distance(extra0, extra1);
}

ListPaths read_param_file(std::string paramfile);
void add_flg(std::map<Lst, int>& flgmap, int gid, int cid, int flg);
void read_list_viewercheck(std::string file, std::map<Lst, int>& flgmap);
void read_list_simple(std::string file, std::map<Lst, int>& flgmap, int flg);
bool is_muon(Momentum_recon::Mom_chain& c);
void get_vertex_side_point(Momentum_recon::Mom_chain& c, int vertex_pl, double& ax, double& ay, double& x, double& y, double& z, double& vph);
void Output_partner_information(std::vector<Momentum_recon::Event_information>& momch, std::map<Lst, int>& flgmap, int nseg_cut_thr, bool accept_specified, std::string output);

int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "===============================================================================\n");
		fprintf(stderr, " usage : prg in.momch param.txt out.txt\n\n");
		fprintf(stderr, " param.txt : one \"key path\" per line (space/tab separated, rest of the\n");
		fprintf(stderr, "             line after the key is treated as the path). unused keys may\n");
		fprintf(stderr, "             be omitted. lines starting with # are comments.\n\n");
		fprintf(stderr, "   viewer_check <path>    -> reject_flg = 1 (list format : (\"%%d%%04d1\",groupid,chainid))\n");
		fprintf(stderr, "   VE           <path>    -> reject_flg = 2 (list format : groupid chainid)\n");
		fprintf(stderr, "   manual_check <path>    -> reject_flg = 3 (list format : groupid chainid)\n");
		fprintf(stderr, "   range_mom_cut<path>    -> reject_flg = 4 (list format : groupid chainid)\n");
		fprintf(stderr, "   10seg_cut    <threshold(int)>  -> reject_flg = 5 if nseg < threshold\n");
		fprintf(stderr, "                                     (applied only to mip tracks, i.e. particle_flg != 2212 and != 13)\n");
		fprintf(stderr, "   (reject_flg = 6 if nseg < 26 + 4.5*ln(MCS-momentum/500), applied only to mip tracks)\n");
		fprintf(stderr, "   accept       <path>    -> reject_flg = 0 (list format : groupid chainid)\n");
		fprintf(stderr, "                              when \"accept\" is specified, any partner not listed\n");
		fprintf(stderr, "                              anywhere in the param file gets reject_flg = 7\n");
		fprintf(stderr, "                              (poor convergence from the vertex point).\n");
		fprintf(stderr, "                              when \"accept\" is omitted, unlisted partners get reject_flg = 0\n");
		fprintf(stderr, "                              (same behavior as before).\n");
		fprintf(stderr, "   (if more than one of the above applies to a track, the smallest reject_flg wins.\n");
		fprintf(stderr, "    tracks not matched by anything get reject_flg = 0 unless \"accept\" is specified)\n");
		fprintf(stderr, "===============================================================================\n");
		exit(1);
	}
	std::string in_momch = argv[1];
	std::string param_file = argv[2];
	std::string out_txt = argv[3];

	ListPaths lp = read_param_file(param_file);

	std::map<Lst, int> flgmap;
	if (!lp.viewer_check.empty())  read_list_viewercheck(lp.viewer_check, flgmap);
	if (!lp.VE.empty())            read_list_simple(lp.VE, flgmap, 2);
	if (!lp.manual_check.empty())  read_list_simple(lp.manual_check, flgmap, 3);
	if (!lp.range_mom_cut.empty()) read_list_simple(lp.range_mom_cut, flgmap, 4);
	if (lp.nseg_cut_thr >= 0)      std::cout << "* 10seg_cut : nseg < " << lp.nseg_cut_thr << " -> reject_flg = 5" << std::endl;
	bool accept_specified = !lp.accept.empty();
	if (accept_specified) {
		read_list_simple(lp.accept, flgmap, 0);
		std::cout << "* accept is specified : partners not listed in the param file -> reject_flg = 7" << std::endl;
	}

	std::cout << "* reject flg list size = " << flgmap.size() << std::endl;

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);
	std::cout << "* #of event = " << momch.size() << std::endl;

	Output_partner_information(momch, flgmap, lp.nseg_cut_thr, accept_specified, out_txt);
}

ListPaths read_param_file(std::string paramfile) {
	ListPaths lp;
	std::ifstream ifs(paramfile);
	if (!ifs) {
		std::cerr << "Failed to open param file : " << paramfile << std::endl;
		exit(1);
	}
	std::string line;
	while (std::getline(ifs, line)) {
		while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
		if (line.empty() || line[0] == '#') continue;

		size_t p0 = line.find_first_of(" \t");
		if (p0 == std::string::npos) continue;
		std::string key = line.substr(0, p0);
		size_t p1 = line.find_first_not_of(" \t", p0);
		if (p1 == std::string::npos) continue;
		std::string val = line.substr(p1);

		if (key == "viewer_check")       lp.viewer_check = val;
		else if (key == "VE")            lp.VE = val;
		else if (key == "manual_check")  lp.manual_check = val;
		else if (key == "range_mom_cut") lp.range_mom_cut = val;
		else if (key == "accept")        lp.accept = val;
		else if (key == "10seg_cut") {
			try {
				lp.nseg_cut_thr = std::stoi(val);
			}
			catch (...) {
				std::cerr << "Failed to parse 10seg_cut threshold : " << val << std::endl;
			}
		}
	}
	return lp;
}

void add_flg(std::map<Lst, int>& flgmap, int gid, int cid, int flg) {
	Lst k{ gid, cid };
	auto itr = flgmap.find(k);
	if (itr == flgmap.end()) {
		flgmap.insert(std::make_pair(k, flg));
	}
	else if (flg < itr->second) {
		itr->second = flg;
	}
}

// list format : printf("%d%04d1", groupid, chainid) -> reject_flg = 1 (viewer check)
void read_list_viewercheck(std::string file, std::map<Lst, int>& flgmap) {
	std::ifstream ifs(file);
	if (!ifs) {
		std::cerr << "Failed to open list(viewer check) : " << file << std::endl;
		return;
	}
	std::string str;
	int cnt = 0;
	while (std::getline(ifs, str)) {
		std::vector<std::string> v = StringSplit_with_tab(str);
		if (v.empty()) continue;
		std::string num = v[0];
		if (num.size() < 6) continue;// group(>=1digit) + chain(4digit) + tail(1digit)

		try {
			std::string chain_s = num.substr(num.size() - 5, 4);
			std::string group_s = num.substr(0, num.size() - 5);
			int gid = std::stoi(group_s);
			int cid = std::stoi(chain_s);
			add_flg(flgmap, gid, cid, 1);
			cnt++;
		}
		catch (...) {
			continue;
		}
	}
	std::cout << "* viewer check list : " << cnt << " entries" << std::endl;
}

// list format : groupid chainid -> reject_flg = flg
void read_list_simple(std::string file, std::map<Lst, int>& flgmap, int flg) {
	std::ifstream ifs(file);
	if (!ifs) {
		std::cerr << "Failed to open list(flg=" << flg << ") : " << file << std::endl;
		return;
	}
	std::string str;
	int cnt = 0;
	while (std::getline(ifs, str)) {
		std::vector<std::string> v = StringSplit_with_tab(str);
		if (v.size() < 2) continue;
		try {
			int gid = std::stoi(v[0]);
			int cid = std::stoi(v[1]);
			add_flg(flgmap, gid, cid, flg);
			cnt++;
		}
		catch (...) {
			continue;
		}
	}
	std::cout << "* list(flg=" << flg << ") : " << cnt << " entries" << std::endl;
}

bool is_muon(Momentum_recon::Mom_chain& c) {
	return (c.particle_flg == 13 || c.chainid == 0);
}

// pick the basetrack of the chain nearest to the vertex plate
void get_vertex_side_point(Momentum_recon::Mom_chain& c, int vertex_pl, double& ax, double& ay, double& x, double& y, double& z, double& vph) {
	int dpl = c.base.begin()->pl;
	int upl = c.base.rbegin()->pl;
	if (dpl - vertex_pl > 0) {
		ax = c.base.begin()->ax;
		ay = c.base.begin()->ay;
		x = c.base.begin()->x;
		y = c.base.begin()->y;
		z = c.base.begin()->z;
		vph = c.base.begin()->m[0].ph + c.base.begin()->m[1].ph;
	}
	else {
		ax = c.base.rbegin()->ax;
		ay = c.base.rbegin()->ay;
		x = c.base.rbegin()->x;
		y = c.base.rbegin()->y;
		z = c.base.rbegin()->z;
		vph = c.base.rbegin()->m[0].ph + c.base.rbegin()->m[1].ph;
	}
	(void)upl;
}

void Output_partner_information(std::vector<Momentum_recon::Event_information>& momch, std::map<Lst, int>& flgmap, int nseg_cut_thr, bool accept_specified, std::string output) {
	std::ofstream ofs(output);

	int n_out = 0, n_no_muon = 0;
	for (auto& ev : momch) {

		Momentum_recon::Mom_chain* muon = nullptr;
		for (auto& c : ev.chains) {
			if (is_muon(c)) {
				muon = &c;
				break;
			}
		}

		double mu_ax = 0, mu_ay = 0, mu_x = 0, mu_y = 0, mu_z = 0, mu_vph = 0;
		bool has_muon = (muon != nullptr) && (!muon->base.empty());
		if (has_muon) {
			get_vertex_side_point(*muon, ev.vertex_pl, mu_ax, mu_ay, mu_x, mu_y, mu_z, mu_vph);
		}
		else {
			n_no_muon++;
		}

		// same z-range convention as the existing vertex-clustering tools
		double zrange[2] = { 0,0 };
		if (ev.vertex_pl <= 15 || (ev.vertex_pl >= 16 && ev.vertex_pl % 2 == 0)) {
			zrange[0] = -1000;
		}
		else if (ev.vertex_pl % 2 == 1) {
			zrange[0] = -3200;
		}

		for (auto& c : ev.chains) {
			if (muon != nullptr && &c == muon) continue;
			if (c.base.empty()) continue;

			double ax, ay, x, y, z, vph;
			get_vertex_side_point(c, ev.vertex_pl, ax, ay, x, y, z, vph);

			int nseg = (int)c.base.size();
			int npl = c.base.rbegin()->pl - c.base.begin()->pl + 1;

			double md = -1, oa = -1;
			if (has_muon) {
				matrix_3D::vector_3D pos0, pos1, dir0, dir1;
				pos0.x = mu_x; pos0.y = mu_y; pos0.z = mu_z;
				dir0.x = mu_ax; dir0.y = mu_ay; dir0.z = 1;
				pos1.x = x; pos1.y = y; pos1.z = z;
				dir1.x = ax; dir1.y = ay; dir1.z = 1;

				double extra[2];
				double zr[2] = { zrange[0], zrange[1] };
				md = minimum_distance_fixed(pos0, pos1, dir0, dir1, zr, extra, mu_z);
				oa = matrix_3D::opening_angle(dir0, dir1);
			}

			// [0]:muon-hypothesis momentum, [1]:proton-hypothesis momentum
			int mom_idx = (c.particle_flg == 2212) ? 1 : 0;
			double mcs_mom = c.ecc_mcs_mom[mom_idx];
			double range_mom = c.ecc_range_mom[mom_idx];

			int reject_flg = 0;
			bool has_flg = false;
			Lst k{ ev.groupid, c.chainid };
			auto fitr = flgmap.find(k);
			if (fitr != flgmap.end()) {
				reject_flg = fitr->second;
				has_flg = true;
			}
			else if (accept_specified) {
				// not listed anywhere in the param file -> poor convergence from the vertex point
				reject_flg = 7;
				has_flg = true;
			}
			// reject_flg = 5 and 6 take top priority (5 first, then 6), applied last
			// and overriding whatever was determined above by the other conditions.
			bool is_mip = (c.particle_flg != 2212 && c.particle_flg != 13);
			if (nseg_cut_thr >= 0 && is_mip && nseg < nseg_cut_thr) {
				reject_flg = 5;
				has_flg = true;
			}
			// mip cut : nseg < 26 + 4.5*ln(MCS-momentum/500) -> reject_flg = 6
			else if (is_mip && mcs_mom > 0) {
				double thr6 = 26.0 + 4.5 * std::log(mcs_mom / 500.0);
				if (nseg < thr6) {
					reject_flg = 6;
					has_flg = true;
				}
			}

			ofs << std::right << std::fixed
				<< std::setw(6)  << ev.groupid << " "
				<< std::setw(5)  << c.chainid << " "
				<< std::setw(3)  << ev.vertex_material << " "
				<< std::setw(5)  << c.particle_flg << " "
				<< std::setw(9)  << std::setprecision(4) << ax << " "
				<< std::setw(9)  << std::setprecision(4) << ay << " "
				<< std::setw(11) << std::setprecision(1) << x << " "
				<< std::setw(11) << std::setprecision(1) << y << " "
				<< std::setw(11) << std::setprecision(1) << z << " "
				<< std::setw(7)  << std::setprecision(0) << vph << " "
				<< std::setw(4)  << nseg << " "
				<< std::setw(4)  << npl << " "
				<< std::setw(9)  << std::setprecision(2) << md << " "
				<< std::setw(7)  << std::setprecision(4) << oa << " "
				<< std::setw(9)  << std::setprecision(1) << mcs_mom << " "
				<< std::setw(9)  << std::setprecision(1) << range_mom << " "
				<< std::setw(2)  << reject_flg << std::endl;
			n_out++;
		}
	}
	std::cout << "* #of output(partner) track = " << n_out << std::endl;
	std::cout << "* #of event without muon = " << n_no_muon << std::endl;
}
