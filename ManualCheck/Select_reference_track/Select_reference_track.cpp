std::vector<Difference> Calc_difference(std::vector<Momentum_recon::Event_information>& momch_ev, std::map<int, std::map<int, Momentum_recon::Mom_chain>>& momch_mom_map, std::map<int, std::vector<Fiducial_Area::Fiducial_Area>>& area) {
	std::vector<Difference> ret;

	for (auto& ev : momch_ev) {
		Difference diff;
		diff.groupid = ev.groupid;

		std::map<int, Momentum_recon::Mom_chain> mom_map, mom_map_inv;
		if (momch_mom_map.count(diff.groupid) == 0) {
			fprintf(stderr, "event %d partner not found\n", diff.groupid);
		}
		if (momch_mom_map.count(diff.groupid * -1) == 0) {
			fprintf(stderr, "event %d partner inv not found\n", diff.groupid);
		}
		mom_map = momch_mom_map.at(diff.groupid);
		mom_map_inv = momch_mom_map.at(-1 * diff.groupid);
		for (auto& c : ev.chains) {
			if (c.chainid == 0)continue;

			diff.vertex_pl = ev.vertex_pl;
			diff.direction = c.direction;
			diff.chainid = c.chainid;
			diff.edge_pl[0] = c.base.begin()->pl;
			diff.edge_pl[1] = c.base.rbegin()->pl;
			diff.edge_out_flg = judge_momch_edgeout(diff.direction, c, area, 0, 4);
			diff.add_nseg = 0;
			diff.target_pl[0] = -1;
			diff.target_pl[1] = -1;
			diff.s_pb_angle = 0;
			diff.s_pb_position = 0;
			diff.mcs_pb[0] = -1;
			diff.mcs_pb[1] = -1;
			diff.mcs_pb_error[0][0] = -1;
			diff.mcs_pb_error[0][1] = -1;
			diff.mcs_pb_error[1][0] = -1;
			diff.mcs_pb_error[1][1] = -1;
			std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack> divide_pair;
			//forward
			if (c.direction == 1) {
				diff.target_pl[0] = diff.vertex_pl;

				for (auto& b : c.base) {
					if (b.pl > diff.vertex_pl) {
						if (diff.target_pl[1] < 0) {
							diff.target_pl[1] = b.pl;
						}
						diff.add_nseg += 1;
					}
				}
				if (diff.target_pl[1] < 0)continue;
				for (auto& pair : c.base_pair) {
					if (pair.first.pl == diff.target_pl[0] && pair.second.pl == diff.target_pl[1]) {
						Calc_difference(diff, pair);
						divide_pair = pair;
						break;
					}
				}

			}
			//backward
			else if (c.direction == -1) {
				for (auto b = c.base.rbegin(); b != c.base.rend(); b++) {
					if (b->pl <= diff.vertex_pl) {
						if (diff.target_pl[1] < 0) {
							diff.target_pl[0] = b->pl;
							diff.target_pl[1] = std::next(b, -1)->pl;
						}
						diff.add_nseg += 1;
					}
				}
				if (diff.target_pl[1] < 0)continue;
				for (auto& pair : c.base_pair) {
					if (pair.first.pl == diff.target_pl[0] && pair.second.pl == diff.target_pl[1]) {
						Calc_difference(diff, pair);
						divide_pair = pair;
						break;
					}
				}

			}

			double pb_error[2];
			Momentum_recon::Mom_chain mom, mom_inv;
			if (mom_map.count(c.chainid) == 0) {
				fprintf(stderr, "event %d chain %d partner not found\n", diff.groupid, c.chainid);
			}
			else {
				mom = mom_map.at(c.chainid);
				diff.mcs_pb[0] = mom.Get_muon_mcs_pb();
				mom.Get_muon_pb_mcs_error(pb_error);
				diff.mcs_pb_error[0][0] = pb_error[0];
				diff.mcs_pb_error[0][1] = pb_error[1];
			}

			if (mom_map_inv.count(c.chainid) == 0) {
				fprintf(stderr, "event %d chain %d partner inv not found\n", diff.groupid, c.chainid);
				continue;
			}
			else {
				mom_inv = mom_map_inv.at(c.chainid);
				diff.mcs_pb[1] = mom_inv.Get_muon_mcs_pb();
				mom_inv.Get_muon_pb_mcs_error(pb_error);
				diff.mcs_pb_error[1][0] = pb_error[0];
				diff.mcs_pb_error[1][1] = pb_error[1];
			}
			if (count_pb_angle_diff(mom) < 2 && count_pb_angle_diff(mom_inv) < 2) {
				Calc_sigma_pb(diff, 2, divide_pair);
			}
			else if (count_pb_angle_diff(mom) < count_pb_angle_diff(mom_inv)) {
				Calc_sigma_pb(diff, 1, divide_pair);
			}
			else {
				Calc_sigma_pb(diff, 0, divide_pair);
			}
			ret.push_back(diff);
		}
	}
	return ret;

}
