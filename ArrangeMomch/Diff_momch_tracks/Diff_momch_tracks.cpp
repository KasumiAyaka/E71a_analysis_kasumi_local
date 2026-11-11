// 2026/11/06
// Diff_momch_tracks
// kasumi
//
// Compare two momch files (before / after a chain of processing steps) and
// list the partner tracks (muon excluded) that existed in "before" but no
// longer exist in "after", with the same per-track feature columns as
// Output_rejected_tracks.exe.
//
// output format :
//   groupid chainid vertex-material pid ax ay x y z vph nseg npl md oa MCS-momentum Range-momentum water_pene_tag
//
// water_pene_tag (only meaningful if a water_pene.txt-style file is given) :
//   0 : not listed in water_pene.txt
//   1 : listed in water_pene.txt without "p"  (fully removed by viewer check)
//   2 : listed in water_pene.txt with "p"     (reclassified as proton, not a plain reject)
//  -1 : no water_pene.txt given

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}

double minimum_distance_fixed(matrix_3D::vector_3D pos0, matrix_3D::vector_3D pos1, matrix_3D::vector_3D dir0, matrix_3D::vector_3D dir1, double z_range[2], double extra[2], double refz) {
	double extra0_distance, extra1_distance, delta;
	matrix_3D::vector_3D pos;
	pos.x = pos1.x - pos0.x;
	pos.y = pos1.y - pos0.y;
	pos.z = pos1.z - pos0.z;
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

bool is_muon(Momentum_recon::Mom_chain& c) {
	return (c.particle_flg == 13 || c.chainid == 0);
}

// pick the basetrack of the chain nearest to the vertex plate
void get_vertex_side_point(Momentum_recon::Mom_chain& c, int vertex_pl, double& ax, double& ay, double& x, double& y, double& z, double& vph) {
	int dpl = c.base.begin()->pl;
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
}

std::map<Lst, int> read_water_pene(std::string filename) {
	std::map<Lst, int> ret;
	std::ifstream ifs(filename);
	std::string str;
	while (std::getline(ifs, str)) {
		std::vector<std::string> v = StringSplit(str);
		if (v.size() != 1 && v.size() != 2) continue;
		try {
			int id_all = std::stoi(v[0]);
			Lst k{ id_all / 100000, (id_all % 100000) / 10 };
			ret[k] = (v.size() == 2) ? 2 : 1;
		}
		catch (...) { continue; }
	}
	return ret;
}

int main(int argc, char** argv) {
	if (argc != 4 && argc != 5) {
		fprintf(stderr, "usage : prg before.momch after.momch out.txt [water_pene.txt]\n");
		exit(1);
	}
	std::string file_before = argv[1];
	std::string file_after = argv[2];
	std::string file_out = argv[3];
	std::string file_water_pene = (argc == 5) ? argv[4] : "";

	std::vector<Momentum_recon::Event_information> before = Momentum_recon::Read_Event_information_extension(file_before);
	std::vector<Momentum_recon::Event_information> after = Momentum_recon::Read_Event_information_extension(file_after);

	std::set<Lst> surviving;
	for (auto& ev : after) {
		for (auto& c : ev.chains) {
			surviving.insert(Lst{ ev.groupid, c.chainid });
		}
	}

	std::map<Lst, int> water_pene;
	bool has_water_pene = !file_water_pene.empty();
	if (has_water_pene) water_pene = read_water_pene(file_water_pene);

	std::ofstream ofs(file_out);

	int n_out = 0, n_before = 0, n_after = 0;
	for (auto& ev : before) {
		Momentum_recon::Mom_chain* muon = nullptr;
		for (auto& c : ev.chains) {
			if (is_muon(c)) { muon = &c; break; }
		}
		double mu_ax = 0, mu_ay = 0, mu_x = 0, mu_y = 0, mu_z = 0, mu_vph = 0;
		bool has_muon = (muon != nullptr) && (!muon->base.empty());
		if (has_muon) {
			get_vertex_side_point(*muon, ev.vertex_pl, mu_ax, mu_ay, mu_x, mu_y, mu_z, mu_vph);
		}

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
			n_before++;

			Lst k{ ev.groupid, c.chainid };
			if (surviving.count(k) != 0) continue; // still present in "after" -> not removed

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

			int mom_idx = (c.particle_flg == 2212) ? 1 : 0;
			double mcs_mom = c.ecc_mcs_mom[mom_idx];
			double range_mom = c.ecc_range_mom[mom_idx];

			int wp_tag = -1;
			if (has_water_pene) {
				auto itr = water_pene.find(k);
				wp_tag = (itr != water_pene.end()) ? itr->second : 0;
			}

			ofs << std::right << std::fixed
				<< std::setw(6) << ev.groupid << " "
				<< std::setw(5) << c.chainid << " "
				<< std::setw(3) << ev.vertex_material << " "
				<< std::setw(5) << c.particle_flg << " "
				<< std::setw(9) << std::setprecision(4) << ax << " "
				<< std::setw(9) << std::setprecision(4) << ay << " "
				<< std::setw(11) << std::setprecision(1) << x << " "
				<< std::setw(11) << std::setprecision(1) << y << " "
				<< std::setw(11) << std::setprecision(1) << z << " "
				<< std::setw(7) << std::setprecision(0) << vph << " "
				<< std::setw(4) << nseg << " "
				<< std::setw(4) << npl << " "
				<< std::setw(9) << std::setprecision(2) << md << " "
				<< std::setw(7) << std::setprecision(4) << oa << " "
				<< std::setw(9) << std::setprecision(1) << mcs_mom << " "
				<< std::setw(9) << std::setprecision(1) << range_mom << " "
				<< std::setw(2) << wp_tag << std::endl;
			n_out++;
		}
	}
	for (auto& ev : after) {
		for (auto& c : ev.chains) {
			if (is_muon(c)) continue;
			n_after++;
		}
	}

	std::cout << "* #of partner tracks in before-momch = " << n_before << std::endl;
	std::cout << "* #of partner tracks in after-momch  = " << n_after << std::endl;
	std::cout << "* #of removed partner tracks         = " << n_out << std::endl;
}
