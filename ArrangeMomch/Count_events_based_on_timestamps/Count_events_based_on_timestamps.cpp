// 2026/04/28
// Count_events_based_on_timestamps
// kasumi
// Timestamp-based event counting

#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
#include <iostream>
#include <random>

#include <list>
#include <cassert>
#include <filesystem>

struct Key {
	int64_t utime;
	int bunch;
};
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.utime, lhs.bunch) < std::tie(rhs.utime, rhs.bunch);
}

struct Event {
	int tnum, material, vpl, flg;
};
struct Materials {
	int wt, bs, fe, em, env, pene, side, samet, stop, oth;
	// water,base,iton,emulsion,envelop,penetrate(sandmuon),sideout,sametimestamp,stop
};


void MeasureProcessingTime(std::chrono::system_clock::time_point& start, std::chrono::system_clock::time_point& end) {
	auto dur = end - start;        // 要した時間を計算
	auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
	// 要した時間をミリ秒（1/1000秒）に変換して表示
	std::cout << msec << " milli sec \n";
	if (msec / 1000 < 60) {
		std::cout << msec / 1000 << "sec\n";
	}
	else if (msec / 1000 / 60 < 60) {
		std::cout << msec / 1000 / 60 << "min\n";
	}
	else if (msec / 1000 / 3600 < 24) {
		std::cout << msec / 1000 / 3600 << "h\n";
	}
	else {
		std::cout << (msec / 1000 / 3600) / 24 << "day" << (msec / 1000 / 3600) % 24 << "h\n";
	}
};

double minimum_distance_fixed(matrix_3D::vector_3D pos0, matrix_3D::vector_3D pos1, matrix_3D::vector_3D dir0, matrix_3D::vector_3D dir1, double z_range[2], double extra[2], double refz) {
	double extra0_distance, extra1_distance, delta;
	matrix_3D::vector_3D pos;
	pos.x = pos1.x - pos0.x;
	pos.y = pos1.y - pos0.y;
	pos.z = pos1.z - pos0.z;
	//ほぼ平行な場合
	if (opening_angle(dir0, dir1) < 0.0001) {
		extra0_distance = (pos1.z + pos0.z) / 2 - pos0.z;
		extra1_distance = (pos1.z + pos0.z) / 2 - pos1.z;
	}
	else {
		delta = dot(dir0, dir0) * dot(dir1, dir1) - pow(dot(dir0, dir1), 2.);
		extra0_distance = (+1 * dot(pos, dir0) * dot(dir1, dir1) - dot(dir0, dir1) * dot(pos, dir1)) / delta;
		extra1_distance = (-1 * dot(pos, dir1) * dot(dir0, dir0) + dot(dir0, dir1) * dot(pos, dir0)) / delta;
	}
	//range[0]:小,range[1]:大
	if (z_range[0] > z_range[1]) {
		double tmp_d = z_range[0];
		z_range[0] = z_range[1];
		z_range[1] = tmp_d;
	}

	matrix_3D::vector_3D extra0 = addition(pos0, const_multiple(dir0, extra0_distance));
	matrix_3D::vector_3D extra1 = addition(pos1, const_multiple(dir1, extra1_distance));

	if (extra0.z < refz + z_range[0] || extra1.z < refz + z_range[0]) {//2025/8/20 fixed
		extra0_distance = refz - pos0.z + z_range[0];
		extra1_distance = refz - pos1.z + z_range[0];
	}
	else if (extra0.z > refz + z_range[1] || extra1.z > refz + z_range[1]) {//2025/8/20 fixed
		extra0_distance = refz - pos0.z + z_range[1];
		extra1_distance = refz - pos1.z + z_range[1];
	}

	extra[0] = extra0_distance;
	extra[1] = extra1_distance;
	extra0 = addition(pos0, const_multiple(dir0, extra0_distance));
	extra1 = addition(pos1, const_multiple(dir1, extra1_distance));

	return distance(extra0, extra1);

}

void Set_events_by_timestamp(std::vector<Momentum_recon::Event_information>& momch, std::multimap<Key, Event>& map, std::set<Key>& set);

void Count_events(std::multimap<Key, Event>& map, std::set<Key>& set);
void Count_events_listed_timestamp(std::multimap<Key, Event>& map, std::set<Key>& set, std::set<int64_t>& list);


int main(int argc, char** argv) {
	if (argc < 2 || argc>6) {
		fprintf(stderr, "===============================================================================\n");
		fprintf(stderr, " usage:prg in.momch [in.txt]\n");
		fprintf(stderr, "===============================================================================\n");
		exit(1);
	}
	std::string in_momch = argv[1];// input momch
	std::string in_list = "";
	std::set<int64_t> timestamplist;
	if (argc > 2) {
		in_list = argv[2];
		int64_t utime;
		std::ifstream ifs(in_list);
		while (ifs >> utime) {
			timestamplist.insert(utime);
		}
	}

	//bool result = std::filesystem::create_directories(out_momch);
	//assert(result);
	//assert(std::filesystem::exists(out_momch));//directryの存在確認
	//assert(std::filesystem::is_directory(out_momch));//指定されたパスがディレクトリを指しているかを確認する

	auto start = std::chrono::system_clock::now();//for measure working time
	// read momch
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);
	std::multimap<Key, Event>map;
	std::set<Key> set;
	Set_events_by_timestamp(momch, map, set);
	if (argc == 2) {
		Count_events(map, set);
	}
	else {
		Count_events_listed_timestamp(map, set, timestamplist);
	}
	auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
	MeasureProcessingTime(start, end);
}

void Set_events_by_timestamp(std::vector<Momentum_recon::Event_information>& momch, std::multimap<Key,Event> &map, std::set<Key>& set) {

	Key k;
	Event e;
	for (auto& ev : momch) {
		k.utime = ev.unix_time;
		k.bunch = int(ev.weight);
		e.material = ev.vertex_material;
		e.tnum = ev.chains.size();
		e.vpl = ev.vertex_pl;
		map.insert(std::make_pair(k, e));
		set.insert(k);
	}

	std::cout << " # of events(timestamp, bunch based) = " << set.size() << "/" << momch.size() << std::endl;



}
void Count_events(std::multimap<Key, Event>& map, std::set<Key>& set) {

	Materials m,sum;
	sum = { 0 };
	for (auto itr = set.begin(); itr!=set.end(); itr++) {
		m = { 0 };
		auto p=map.equal_range(*itr);
		for (auto q = p.first; q != p.second; q++) {
			if (q->second.material == 0) {
				m.wt++;
			}
			else if (q->second.material == 1) {
				m.bs++;
			}
			else if (q->second.material == 2) {
				m.fe++;
			}
			else if (q->second.material == 5) {
				m.em++;
			}
			else if (q->second.material == 6) {
				m.env++;
			}
			else if (q->second.material == 7) {
				m.stop++;
			}
			else if (q->second.material == -2) {
				m.pene++;
			}
			else if (q->second.material == -3) {
				m.side++;
			}
			else if (q->second.material == -5) {
				m.samet++;
			}
			else {
				if (q->second.material % 10 == 2) {
					m.fe++;
				}

				m.oth++;
			}
			//std::wcout << q->second.material << std::endl;
	}
		std::cout <<itr->utime<<", "<< itr->bunch
			<< "  ( water, base, iron, emulsion, envelop, pene, side, same timestamp, stop, others ) = " 
			<< std::setw(3) << m.wt << ", " 
			<< std::setw(3) << m.bs << ", "
			<< std::setw(3) << m.fe << ", "
			<< std::setw(3) << m.em << ", "
			<< std::setw(3) << m.env << ", "
			<< std::setw(3) << m.pene << ", "
			<< std::setw(3) << m.side << ", "
			<< std::setw(3) << m.samet << ", "
			<< std::setw(3) << m.stop << ", "
			<< std::setw(3) << m.oth << " " << std::endl;

		if (m.wt > 0) {
			sum.wt++;
		}
		else {// except water events
			if (m.fe > 0) {
				sum.fe++;
			}
			else {
				if (m.bs > 0 || m.em > 0 || m.env > 0 || m.stop > 0) {
					sum.stop++;
				}
				else {
					if (m.samet > 0) {
						sum.samet++;
					}
					else {
						if (m.pene > 0) {
							sum.pene++;
						}
						else {
							if (m.side > 0) {
								sum.side++;
							}
							else {
								std::cout << " unexpected caterory : " << itr->utime << std::endl;
								sum.oth++;
							}
						}
					}
				}
			}
		}
	}
	std::cout << " # of events (uniq / all)= " << set.size() << " / " << map.size() << std::endl;
	std::cout << "  ( water, iron, stop, pene, side, same timestamp ) = "
		<< std::setw(3) << sum.wt << ", "
		<< std::setw(3) << sum.fe << ", "
		<< std::setw(3) << sum.stop << ", "
		<< std::setw(3) << sum.pene << ", "
		<< std::setw(3) << sum.side << ", "
		<< std::setw(3) << sum.samet << ", "
		<< std::setw(3) << sum.oth << " " << std::endl;
	std::wcout << std::setw(3) << sum.wt + sum.fe + sum.stop + sum.pene + sum.side + sum.samet + sum.oth << std::endl;

}
void Count_events_listed_timestamp(std::multimap<Key, Event>& map, std::set<Key>& set,std::set<int64_t>&list) {

	Materials m, sum;
	sum = { 0 };
	for (auto itr = set.begin(); itr != set.end(); itr++) {
		m = { 0 };
		auto p = map.equal_range(*itr);
		if (list.count(itr->utime) == 0)continue;
		for (auto q = p.first; q != p.second; q++) {
			if (q->second.material == 0) {
				m.wt++;
			}
			else if (q->second.material == 1) {
				m.bs++;
			}
			else if (q->second.material == 2) {
				m.fe++;
			}
			else if (q->second.material == 5) {
				m.em++;
			}
			else if (q->second.material == 6) {
				m.env++;
			}
			else if (q->second.material == 7) {
				m.stop++;
			}
			else if (q->second.material == -2) {
				m.pene++;
			}
			else if (q->second.material == -3) {
				m.side++;
			}
			else if (q->second.material == -5) {
				m.samet++;
			}
			else {
				m.oth++;
			}
		}
		std::cout << itr->utime << ", " << itr->bunch
			<< "  ( water, base, iron, emulsion, envelop, pene, side, same timestamp, stop, others ) = "
			<< std::setw(3) << m.wt << ", "
			<< std::setw(3) << m.bs << ", "
			<< std::setw(3) << m.fe << ", "
			<< std::setw(3) << m.em << ", "
			<< std::setw(3) << m.env << ", "
			<< std::setw(3) << m.pene << ", "
			<< std::setw(3) << m.side << ", "
			<< std::setw(3) << m.samet << ", "
			<< std::setw(3) << m.stop << ", " << 
			std::setw(3) << m.oth << " " << std::endl;

		if (m.wt > 0) {
			sum.wt++;
		}
		else {// except water events
			if (m.fe > 0) {
				sum.fe++;
			}
			else {
				if (m.bs > 0 || m.em > 0 || m.env > 0 || m.stop > 0) {
					sum.stop++;
				}
				else {
					if (m.samet > 0) {
						sum.samet++;
					}
					else {
						if (m.pene > 0) {
							sum.pene++;
						}
						else {
							if (m.side > 0) {
								sum.side++;
							}
							else {
								std::cout << " unexpected caterory : " << itr->utime << std::endl;
								sum.oth++;
							}
						}
					}
				}
			}
		}
	}
	std::cout << " # of events (uniq / all)= " << set.size() << " / " << map.size() << std::endl;
	std::cout << "  ( water, iron, stop, pene, side, same timestamp ) = "
		<< std::setw(3) << sum.wt << ", "
		<< std::setw(3) << sum.fe << ", "
		<< std::setw(3) << sum.stop << ", "
		<< std::setw(3) << sum.pene << ", "
		<< std::setw(3) << sum.side << ", "
		<< std::setw(3) << sum.samet << ", " 
	<< std::setw(3) << sum.oth << " " << std::endl;

}
