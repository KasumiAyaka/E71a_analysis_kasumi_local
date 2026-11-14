// 2026/11/09
// Partner_penetrate_eval
// kasumi
//
// For every partner chain that the viewer-check group mfile (m_ev_group.all,
// produced by event_group_viewer_check.exe) found at least one vertex-side
// "attach candidate" chain for (the chains the viewer overlays on the "...1"
// id variant to judge whether the partner is actually one half of a single
// physical track that penetrates the vertex plate), evaluate the candidate
// against the partner using the ORIGINAL momch (which still has momentum):
//   - position offset  between partner's vertex-side point and the
//     candidate's nearest-to-vertex point
//   - angle offset     between the two direction vectors
//   - candidate's own md (minimum distance to muon), same formula as
//     Diff_momch_tracks
//   - momentum of the candidate, IF the candidate's near-vertex basetrack
//     (pl,rawid) can be matched back to an actual chain in the momch
//     (event_group_viewer_check.exe renumbers chain_id/group_id when it
//     writes m_ev_group.all, so the candidate's true identity is recovered
//     only via its basetrack rawid, not its chain_id field)
//
// output format (one row per partner x candidate pair) :
//   groupid chainid partner_vph partner_mcsmom partner_rangemom partner_md partner_oa
//   cand_idx cand_vph cand_pl pos_diff angle_diff cand_md
//   cand_matched cand_chainid cand_mcsmom cand_rangemom mom_diff
//   water_pene_tag
//   pc_md partner_vph_slope cand_vph_slope partner_linklet_std cand_linklet_std
//   dal dar rms_l rms_r chi2_al chi2_ar sigma_vph chi2_vph
//
// cand_matched : 1 if the candidate's near-vertex basetrack was found inside
//                a real momch chain in this event (so cand_mcsmom/rangemom/
//                mom_diff are meaningful), 0 otherwise (-1 sentinels used)
// water_pene_tag (partner, same scheme as Diff_momch_tracks) :
//   0 : not listed in water_pene.txt   1 : listed without "p" (reject)
//   2 : listed with "p" (proton reclass)   -1 : no water_pene.txt given
//
// pc_md               : minimum distance between the partner's line and the
//                        candidate's line (same extrapolation formula as
//                        partner_md/cand_md, but partner-vs-candidate instead
//                        of muon-vs-X) -- near 0 if they are collinear (one
//                        physical track), large if they are unrelated tracks
// {partner,cand}_vph_slope : least-squares slope of VPH vs distance-from-vertex
//                        (plates) along that chain's own basetracks. Positive
//                        slope = VPH grows away from the vertex (Bragg-like,
//                        reaction-origin suggestive); near 0/negative = flat or
//                        falling (more consistent with a through-going track)
// {partner,cand}_linklet_std : stddev of the plate-to-plate direction-angle
//                        difference ("linklet" scattering angle) along that
//                        chain's own basetracks, sorted by plate. Reflects MCS
//                        strength (~1/momentum); if partner+candidate are one
//                        physical track they should have similar momentum and
//                        thus similar linklet_std. -1 sentinel if <3 basetracks
//                        (not enough segments for a stddev)
//
// The following 8 columns reproduce, for the partner/candidate junction, the
// exact same lateral/radial-angle-jump and VPH-sigma methodology that
// 01_check_upstream/src/Check_upstream_base_analysis.cpp (Calc_divide_angle_lateral/
// radial, Calc_divide_vph) uses to judge whether a split point WITHIN a single
// muon chain looks like one continuous track. Here the "split point" is the
// partner-candidate boundary instead of a point inside one chain:
//   dal, dar       : lateral/radial decomposition of the direction jump from
//                     the partner's vertex-side basetrack to the candidate's
//                     near-vertex basetrack (same decomposition formula)
//   rms_l, rms_r   : rms of that same lateral/radial decomposition computed
//                     over the PARTNER's own internal consecutive-basetrack
//                     pairs -- i.e. the partner's own intrinsic scattering
//                     scale, used as the reference noise level. NAN if the
//                     partner has <3 basetracks (not enough internal pairs)
//   chi2_al, chi2_ar : (dal/rms_l)^2, (dar/rms_r)^2 -- large means the jump at
//                     the junction is much bigger than the partner's own
//                     typical segment-to-segment scattering (suggests a real
//                     break, i.e. two different tracks); small is consistent
//                     with one continuous track. Same thr_ang~100 convention
//                     as the reference code (compare chi2_al+chi2_ar to ~100)
//   sigma_vph, chi2_vph : (cand_vph - partner_vph_mean)/partner_vph_sd and its
//                     square, where partner_vph_mean/sd come from the
//                     partner's own basetracks (both films, same formula as
//                     Calc_divide_vph). Same thr_vph~25 convention
//
// Two more files are written alongside out.txt (same stem + suffix), holding
// SELF-ONLY (never cross partner/candidate) chi2 distributions, for directly
// comparing "does the candidate look internally as straight/as VPH-consistent
// as the partner" rather than testing the partner-candidate junction itself:
//   <stem>_self_radlat.txt : groupid chainid_partner role cand_idx cand_chainid
//       wp_tag split_idx chi2_al chi2_ar
//     -- role 0 = computed from the partner's OWN internal consecutive-pairs
//        only (one row per internal split point, cand_idx=-1); role 1 = same
//        computation but using ONLY that one candidate's own basetracks
//        (cand_idx = its index, matching out.txt's cand_idx)
//   <stem>_self_vph.txt : groupid chainid_partner role cand_idx cand_chainid
//       wp_tag base_idx chi2_vph
//     -- same role split, one row per basetrack instead of per split point
//
// NOTE: chi2_al/chi2_ar/chi2_vph above are each individually a self-chi2
// built from a chain's OWN rms/sd, normalized by that SAME quantity -- so
// their MEAN over that chain's own points is tautologically ~1 per chi2
// term regardless of track quality (verified: mean chi2_al+chi2_ar across
// 1958 chains = 2.000 +/- 0.004; mean chi2_vph = 1.0-2.0, bounded by the
// finite-sample bias correction only). They are only informative as a
// per-split-point DISTRIBUTION (looking for one outlier point), never as a
// single scalar compared across chains. For a scalar partner-vs-candidate
// comparison, use the file below instead.
//
//   <stem>_chain_scale.txt : groupid chainid_partner role cand_idx
//       cand_chainid wp_tag n_basetrack rms_l rms_r vph_mean vph_sd
//     -- one row per chain (role 0 = partner, role 1 = one candidate), the
//        RAW (not self-normalized) scattering scale and VPH statistics over
//        ALL of that chain's own basetracks. These actually vary with the
//        track's real momentum/particle species, so partner-vs-candidate
//        here (unlike the chi2 mean above) is a meaningful comparison: if
//        one physical track, same momentum => similar rms_l/rms_r, same
//        particle => similar vph_mean/sd.
//
// <stem>_boundary.txt -- a SEPARATE, topology-aware definition of the
// partner-candidate junction point, requested after the original pos_diff/
// angle_diff/pc_md/cand_pl/cand_vph above turned out to independently pick
// "nearest basetrack to vertex_pl" on EACH side, which (when both chains
// happen to have a basetrack right at vertex_pl) could compare the SAME
// plate against itself instead of the two plates straddling the vertex.
// Columns: groupid chainid_partner cand_idx cand_chainid wp_tag valid
//   p_pl cand_pl pl_diff  dal lateral_pos_diff pc_md_boundary
//   p_vph_mean_chain cand_vph_mean_chain vph_mean_diff
// Partner's point = partner_boundary_point(): largest pl<=vertex_pl if the
//   partner has one, else smallest pl overall (always pl>vertex_pl then).
// Candidate's point = candidate_boundary_point(): searched on the OPPOSITE
//   side of vertex_pl from the partner's point (smallest pl>=vertex_pl+1,
//   or largest pl<=vertex_pl) -- i.e. the two points are meant to straddle
//   the vertex plate from either side, never both landing on it.
// valid=0 (all following columns -1) if the candidate chain has no
//   basetrack on the required side at all; these rows are excluded from
//   any aggregate/plot, per instruction.
// dal, lateral_pos_diff : lateral_radial_diff() applied to the direction
//   jump and to the position-difference vector respectively, both relative
//   to the partner's boundary direction (so lateral_pos_diff is the
//   component of the position offset perpendicular to the partner's own
//   track direction, not the raw 3D distance).
// p_vph_mean_chain, cand_vph_mean_chain : VPH mean over ALL of that chain's
//   own basetracks (vec_mean_sd over p_vph_list / c_vph_list), NOT the
//   single boundary basetrack's VPH -- i.e. the whole-chain mean, as
//   requested, not a single-point value.

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#include <map>
#include <set>
#include <algorithm>
#include <cstdio>

struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}
struct PlRaw {
	int pl; int64_t rawid;
};
bool operator<(const PlRaw& lhs, const PlRaw& rhs) {
	return std::tie(lhs.pl, lhs.rawid) < std::tie(rhs.pl, rhs.rawid);
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

// VPH of one basetrack = average of its two films' DECODED ph (ph % 10000
// each), same convention as Calc_divide_vph in Check_upstream_base_analysis.cpp
// ("(m[0].ph % 10000 + m[1].ph % 10000) / 2"), so it is directly comparable
// to a candidate's single merged mfile0::M_Base.ph (also % 10000).
void get_vertex_side_point(Momentum_recon::Mom_chain& c, int vertex_pl, double& ax, double& ay, double& x, double& y, double& z, double& vph, int& pl) {
	int dpl = c.base.begin()->pl;
	if (dpl - vertex_pl > 0) {
		pl = c.base.begin()->pl;
		ax = c.base.begin()->ax; ay = c.base.begin()->ay;
		x = c.base.begin()->x; y = c.base.begin()->y; z = c.base.begin()->z;
		vph = ((double)(c.base.begin()->m[0].ph % 10000) + (double)(c.base.begin()->m[1].ph % 10000)) / 2.0;
	}
	else {
		pl = c.base.rbegin()->pl;
		ax = c.base.rbegin()->ax; ay = c.base.rbegin()->ay;
		x = c.base.rbegin()->x; y = c.base.rbegin()->y; z = c.base.rbegin()->z;
		vph = ((double)(c.base.rbegin()->m[0].ph % 10000) + (double)(c.base.rbegin()->m[1].ph % 10000)) / 2.0;
	}
}

void zrange_for_vertex(int vertex_pl, double zrange[2]) {
	zrange[0] = 0; zrange[1] = 0;
	if (vertex_pl <= 15 || (vertex_pl >= 16 && vertex_pl % 2 == 0)) {
		zrange[0] = -1000;
	}
	else if (vertex_pl % 2 == 1) {
		zrange[0] = -3200;
	}
}

std::map<Lst, int> read_water_pene(std::string filename) {
	std::map<Lst, int> ret;
	if (filename.empty()) return ret;
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

// pick the candidate basetrack nearest (by plate) to vertex_pl
mfile0::M_Base* nearest_base(mfile0::M_Chain& c, int vertex_pl) {
	mfile0::M_Base* best = nullptr;
	int best_d = 1 << 30;
	for (auto& b : c.basetracks) {
		int pl = b.pos / 10;
		int d = std::abs(pl - vertex_pl);
		if (d < best_d) { best_d = d; best = &b; }
	}
	return best;
}

// Partner's boundary basetrack for the partner-candidate junction (NOT the
// muon-vs-partner comparison, which keeps using get_vertex_side_point/p_pl
// unchanged): if the partner chain has any basetrack with pl<=vertex_pl,
// use the one with the LARGEST such pl (closest to the vertex from
// upstream) and report upstream=true; otherwise every basetrack is
// pl>vertex_pl, so use the one with the SMALLEST pl (closest from
// downstream) and report upstream=false. The chain is guaranteed non-empty
// by the caller. Values are returned by reference (not a pointer into
// c.base) since that container may yield const elements.
void partner_boundary_point(Momentum_recon::Mom_chain& c, int vertex_pl, bool& upstream,
	int& out_pl, double& out_x, double& out_y, double& out_z, double& out_ax, double& out_ay) {
	int best_up_pl = -(1 << 30); bool have_up = false;
	double up_x = 0, up_y = 0, up_z = 0, up_ax = 0, up_ay = 0;
	int best_down_pl = 1 << 30; bool have_down = false;
	double down_x = 0, down_y = 0, down_z = 0, down_ax = 0, down_ay = 0;
	for (auto& b : c.base) {
		if (b.pl <= vertex_pl) {
			if (b.pl > best_up_pl) { best_up_pl = b.pl; up_x = b.x; up_y = b.y; up_z = b.z; up_ax = b.ax; up_ay = b.ay; have_up = true; }
		}
		else {
			if (b.pl < best_down_pl) { best_down_pl = b.pl; down_x = b.x; down_y = b.y; down_z = b.z; down_ax = b.ax; down_ay = b.ay; have_down = true; }
		}
	}
	if (have_up) { upstream = true; out_pl = best_up_pl; out_x = up_x; out_y = up_y; out_z = up_z; out_ax = up_ax; out_ay = up_ay; return; }
	upstream = false; out_pl = best_down_pl; out_x = down_x; out_y = down_y; out_z = down_z; out_ax = down_ax; out_ay = down_ay;
}

// Candidate's boundary basetrack, on the OPPOSITE side of vertex_pl from the
// partner's (partner_upstream, from partner_boundary_point above): if the
// partner's point is upstream (pl<=vertex_pl), look downstream for the
// SMALLEST candidate pl with pl>=vertex_pl+1 (the next plate past the
// vertex, or further if that exact plate has no candidate basetrack); if
// the partner's point is downstream, look upstream for the LARGEST
// candidate pl with pl<=vertex_pl. Returns false if the candidate chain has
// no basetrack on that required side at all -- the caller treats that pair
// as having no valid junction point.
bool candidate_boundary_point(mfile0::M_Chain& c, int vertex_pl, bool partner_upstream,
	int& out_pl, double& out_x, double& out_y, double& out_z, double& out_ax, double& out_ay) {
	bool found = false;
	if (partner_upstream) {
		int best_pl = 1 << 30;
		for (auto& b : c.basetracks) {
			int pl = b.pos / 10;
			if (pl >= vertex_pl + 1 && pl < best_pl) {
				best_pl = pl; out_x = b.x; out_y = b.y; out_z = b.z; out_ax = (double)b.ax; out_ay = (double)b.ay; found = true;
			}
		}
		out_pl = best_pl;
	}
	else {
		int best_pl = -(1 << 30);
		for (auto& b : c.basetracks) {
			int pl = b.pos / 10;
			if (pl <= vertex_pl && pl > best_pl) {
				best_pl = pl; out_x = b.x; out_y = b.y; out_z = b.z; out_ax = (double)b.ax; out_ay = (double)b.ay; found = true;
			}
		}
		out_pl = best_pl;
	}
	return found;
}

// pull (plate, VPH, ax, ay) out of a chain's own basetracks, for the
// VPH-vs-vertex-distance slope and the linklet-angle stddev below.
// overload for a real momch chain (Mom_basetrack -- each film's ph still
// needs the %10000 decode, same as get_vertex_side_point above; this was
// previously missing AND summed instead of averaged the two films, which
// inflated partner-side VPH by orders of magnitude relative to a candidate's
// single already-decoded mfile0::M_Base.ph)
void extract_chain_info(Momentum_recon::Mom_chain& c, std::vector<int>& pl, std::vector<double>& vph, std::vector<double>& ax, std::vector<double>& ay) {
	for (auto& b : c.base) {
		pl.push_back(b.pl);
		vph.push_back(((double)(b.m[0].ph % 10000) + (double)(b.m[1].ph % 10000)) / 2.0);
		ax.push_back((double)b.ax);
		ay.push_back((double)b.ay);
	}
}
// overload for a track_sel1-derived mfile chain (M_Base -- ph needs %10000 decode,
// same convention as vtx_output.cpp / vtx_selection.cpp / FILE_structure.cpp)
void extract_chain_info(mfile0::M_Chain& c, std::vector<int>& pl, std::vector<double>& vph, std::vector<double>& ax, std::vector<double>& ay) {
	for (auto& b : c.basetracks) {
		pl.push_back(b.pos / 10);
		vph.push_back((double)(((int64_t)b.ph) % 10000));
		ax.push_back((double)b.ax);
		ay.push_back((double)b.ay);
	}
}

// lateral/radial decomposition of the direction jump (ax0,ay0)->(ax1,ay1),
// same formula as Calc_divide_angle_lateral/radial in Check_upstream_base_analysis.cpp
void lateral_radial_diff(double ax0, double ay0, double ax1, double ay1, double& dal, double& dar) {
	double angle = std::sqrt(ax0 * ax0 + ay0 * ay0);
	double dax = ax1 - ax0, day = ay1 - ay0;
	if (angle < 0.01) { dal = dax; dar = day; }
	else { dal = (dax * ay0 - day * ax0) / angle; dar = (dax * ax0 + day * ay0) / angle; }
}

// Per-(pl0,pl1) plate-pair local alignment correction map (Delaunay-divided
// align_param2 triangles), loaded on demand from the same
// "ali_PPP_QQQ_interpolation.txt" files that mfile_to_momentum_format_partner
// (do_1_before_add_momentum.bat) uses when ORIGINALLY building a chain's
// base_pair -- see the 2026-11-11/12 discussion: m_ev_group.all/momch's
// .base coordinates only carry the single-plate ABSOLUTE alignment, not this
// finer plate-pair-specific one, so for evaluating straightness across a
// specific two-plate junction (partner<->candidate boundary) this per-pair
// map is the physically appropriate one, matching how the chain's own
// internal base_pair entries are built.
std::string g_ecc_path;
// align_param2's corr_p[3] are raw pointers INTO the raw align_param vector
// that DelaunayDivide() was given, so that raw vector must be kept alive for
// as long as the Delaunay-divided one is used -- hence both are cached here,
// never the divided one alone (which would otherwise dangle once the raw
// vector that produced it went out of scope).
std::map<std::pair<int, int>, std::vector<corrmap_3d::align_param>> g_alignmap_raw_cache;
std::map<std::pair<int, int>, std::vector<corrmap_3d::align_param2>> g_alignmap_cache;
std::set<std::pair<int, int>> g_alignmap_missing;

std::vector<corrmap_3d::align_param2>* get_alignmap(int pl0, int pl1) {
	std::pair<int, int> key{ pl0, pl1 };
	auto itr = g_alignmap_cache.find(key);
	if (itr != g_alignmap_cache.end()) return &itr->second;
	if (g_alignmap_missing.count(key)) return NULL;

	char buf[64];
	snprintf(buf, sizeof(buf), "\\Area0\\0\\align\\fine\\ali_%03d_%03d_interpolation.txt", pl0, pl1);
	std::string filename = g_ecc_path + buf;
	std::ifstream check(filename);
	if (!check) { g_alignmap_missing.insert(key); return NULL; }
	check.close();

	auto raw_ins = g_alignmap_raw_cache.insert(std::make_pair(key, corrmap_3d::read_ali_param(filename, false)));
	std::vector<corrmap_3d::align_param2> corr2 = corrmap_3d::DelaunayDivide(raw_ins.first->second);
	auto ins = g_alignmap_cache.insert(std::make_pair(key, corr2));
	return &ins.first->second;
}

// Corrects the junction's "farther" basetrack (whichever of partner/candidate
// has the larger pl) into the "nearer" one's plate-pair-local frame, using
// the SAME library routines (corrmap_3d::track_affineparam_correspondence /
// trans_base_all) that mfile_to_momentum_format_partner uses to build a real
// chain's base_pair.second. The nearer point is used as-is (matching how
// base_pair.first is left untransformed). Returns false if no alignment map
// exists for this (pl0,pl1) plate pair.
bool junction_local_quantities(int pb_pl, double pb_x, double pb_y, double pb_z, double pb_ax, double pb_ay,
	int cb_pl, double cb_x, double cb_y, double cb_z, double cb_ax, double cb_ay,
	double refz, double zr_window[2],
	double& dal, double& dar_out, double& lat_pos, double& rad_pos_out, double& pc_md) {
	bool partner_is_near = pb_pl <= cb_pl;
	int pl0 = partner_is_near ? pb_pl : cb_pl;
	int pl1 = partner_is_near ? cb_pl : pb_pl;
	if (pl0 == pl1) return false;

	std::vector<corrmap_3d::align_param2>* align = get_alignmap(pl0, pl1);
	if (align == NULL || align->empty()) return false;

	double far_z_abs = partner_is_near ? cb_z : pb_z; // real detector z, kept aside: trans_base() overwrites .z with the
	                                                  // alignment file's own small local reference value (plate-gap scale,
	                                                  // not the absolute detector z), which is not usable for pc_md/pos_diff.
	mfile0::M_Base far_pt = mfile0::M_Base();
	far_pt.pos = (partner_is_near ? cb_pl : pb_pl) * 10;
	far_pt.x = partner_is_near ? cb_x : pb_x;
	far_pt.y = partner_is_near ? cb_y : pb_y;
	far_pt.z = far_z_abs;
	far_pt.ax = (float)(partner_is_near ? cb_ax : pb_ax);
	far_pt.ay = (float)(partner_is_near ? cb_ay : pb_ay);
	far_pt.rawid = 0; far_pt.group_id = 0; far_pt.ph = 0;

	std::vector<mfile0::M_Base*> far_v; far_v.push_back(&far_pt);
	std::vector<std::pair<mfile0::M_Base*, corrmap_3d::align_param2*>> track_param =
		corrmap_3d::track_affineparam_correspondence(far_v, *align);
	corrmap_3d::trans_base_all(track_param);
	far_pt.z = far_z_abs; // restore the real absolute z that trans_base_all() overwrote

	double near_x = partner_is_near ? pb_x : cb_x, near_y = partner_is_near ? pb_y : cb_y, near_z = partner_is_near ? pb_z : cb_z;
	double near_ax = partner_is_near ? pb_ax : cb_ax, near_ay = partner_is_near ? pb_ay : cb_ay;

	double near_ax2 = near_ax, near_ay2 = near_ay, far_ax2 = (double)far_pt.ax, far_ay2 = (double)far_pt.ay;
	lateral_radial_diff(near_ax2, near_ay2, far_ax2, far_ay2, dal, dar_out);

	// position residual after extrapolating near's own trajectory by the
	// real z-gap to far's z: predicted_dx/dy = near_ax*dz, near_ay*dz, where
	// dz = far.z - near.z. Subtracting this off before the lateral/radial
	// decomposition is required for the radial component -- see the
	// 2026-11-12 within-chain self-check (mean(dar)=-0.011, essentially
	// unbiased, vs mean(rad_pos)=-4590um on raw (dx,dy) for genuinely
	// continuous same-chain linklets): un-scaled (dx,dy) is dominated by
	// dz*angle (always negative here, since z decreases with increasing pl
	// and far always has the larger pl), not real trajectory deviation.
	// The lateral component is unaffected by this correction (the removed
	// term is purely radial by construction), confirmed algebraically and
	// empirically (lateral_pos_diff showed no such bias).
	double dx = far_pt.x - near_x, dy = far_pt.y - near_y;
	double dz_real = far_pt.z - near_z;
	double pred_dx = near_ax2 * dz_real, pred_dy = near_ay2 * dz_real;
	double resid_dx = dx - pred_dx, resid_dy = dy - pred_dy;
	lateral_radial_diff(near_ax2, near_ay2, near_ax2 + resid_dx, near_ay2 + resid_dy, lat_pos, rad_pos_out);

	// minimum_distance_fixed(), not the raw matrix_3D::minimum_distance(): the
	// raw 3D skew-line formula is numerically unstable for near-parallel
	// directions (the common case for a genuine same-track junction, where
	// dal/dar are small by construction) -- its denominator (|dir0|^2|dir1|^2
	// - (dir0.dir1)^2) goes to ~0, so float noise in dal/dar gets amplified
	// into a huge, physically meaningless extrapolation distance even when
	// the true perpendicular separation (~lat_pos) is tiny. refz/zr_window
	// clamp the extrapolation to the same vertex-centered window already
	// used for partner_md/cand_md/pc_md elsewhere in this file.
	matrix_3D::vector_3D pos0, pos1, dir0, dir1;
	pos0.x = near_x; pos0.y = near_y; pos0.z = near_z;
	dir0.x = near_ax2; dir0.y = near_ay2; dir0.z = 1;
	pos1.x = far_pt.x; pos1.y = far_pt.y; pos1.z = far_pt.z;
	dir1.x = far_pt.ax; dir1.y = far_pt.ay; dir1.z = 1;
	double zr[2] = { zr_window[0], zr_window[1] };
	double extra[2];
	pc_md = minimum_distance_fixed(pos0, pos1, dir0, dir1, zr, extra, refz);
	return true;
}

// rms of the lateral/radial decomposition over a chain's own internal
// consecutive-basetrack pairs (sorted by plate) -- the chain's own intrinsic
// scattering scale, used as the reference noise level for chi2_al/chi2_ar
void chain_linklet_lateral_radial_rms(std::vector<int>& pl, std::vector<double>& ax, std::vector<double>& ay, double& rms_l, double& rms_r) {
	rms_l = NAN; rms_r = NAN;
	size_t n = pl.size();
	if (n < 3) return; // need >=2 internal pairs
	std::vector<size_t> order(n);
	for (size_t i = 0; i < n; i++) order[i] = i;
	std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return pl[a] < pl[b]; });

	double sum2_l = 0, sum2_r = 0;
	int count = 0;
	for (size_t i = 0; i + 1 < n; i++) {
		size_t a = order[i], b = order[i + 1];
		double dal, dar;
		lateral_radial_diff(ax[a], ay[a], ax[b], ay[b], dal, dar);
		sum2_l += dal * dal;
		sum2_r += dar * dar;
		count++;
	}
	if (count <= 1) return;
	rms_l = std::sqrt(sum2_l / count);
	rms_r = std::sqrt(sum2_r / count);
}

// mean/sd of a plain value list (bias-corrected sample sd, same formula as
// the old per-chain VPH helper) -- used for both partner's and a candidate's
// own per-basetrack VPH list (already merged-per-basetrack via extract_chain_info,
// so partner-self and candidate-self chi2_vph below are computed on the exact
// same quantity and are directly comparable).
void vec_mean_sd(std::vector<double>& v, double& mean, double& sd) {
	mean = -1; sd = -1;
	size_t n = v.size();
	if (n < 2) return;
	double sum1 = 0, sum2 = 0;
	for (double x : v) { sum1 += x; sum2 += x * x; }
	mean = sum1 / n;
	double var = sum2 / n - mean * mean;
	if (var <= 0) { mean = -1; sd = -1; return; }
	sd = std::sqrt(var) * std::sqrt((double)(n - 1)) / std::sqrt((double)n);
}

// self-only chi2(rad+lat) distribution over a SINGLE chain's own internal
// consecutive-basetrack pairs: each pair's lateral/radial jump normalized by
// that SAME chain's own overall rms (chain_linklet_lateral_radial_rms above).
// Called once for partner's own basetracks and, separately, once for a
// candidate's own basetracks -- the two resulting distributions never mix
// partner and candidate data, so they can be overlaid to test whether the two
// chains "look similarly straight" internally.
void chain_self_radlat_chi2(std::vector<int>& pl, std::vector<double>& ax, std::vector<double>& ay,
	std::vector<double>& chi2_al_out, std::vector<double>& chi2_ar_out) {
	chi2_al_out.clear(); chi2_ar_out.clear();
	double rms_l, rms_r;
	chain_linklet_lateral_radial_rms(pl, ax, ay, rms_l, rms_r);
	if (!(std::isfinite(rms_l) && rms_l > 0 && std::isfinite(rms_r) && rms_r > 0)) return;
	size_t n = pl.size();
	std::vector<size_t> order(n);
	for (size_t i = 0; i < n; i++) order[i] = i;
	std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return pl[a] < pl[b]; });
	for (size_t i = 0; i + 1 < n; i++) {
		size_t a = order[i], b = order[i + 1];
		double dal, dar;
		lateral_radial_diff(ax[a], ay[a], ax[b], ay[b], dal, dar);
		chi2_al_out.push_back((dal / rms_l) * (dal / rms_l));
		chi2_ar_out.push_back((dar / rms_r) * (dar / rms_r));
	}
}

// self-only chi2(vph) distribution over a single chain's own basetracks:
// each basetrack's VPH normalized by that SAME chain's own mean/sd. Same
// partner-vs-candidate independence as chain_self_radlat_chi2 above.
void chain_self_vph_chi2(std::vector<double>& vph, std::vector<double>& chi2_vph_out) {
	chi2_vph_out.clear();
	double mean, sd;
	vec_mean_sd(vph, mean, sd);
	if (!(sd > 0)) return;
	for (double v : vph) {
		double sigma = (v - mean) / sd;
		chi2_vph_out.push_back(sigma * sigma);
	}
}

double linreg_slope(std::vector<double>& x, std::vector<double>& y) {
	size_t n = x.size();
	if (n < 2) return -999.0;
	double sx = 0, sy = 0, sxx = 0, sxy = 0;
	for (size_t i = 0; i < n; i++) { sx += x[i]; sy += y[i]; sxx += x[i] * x[i]; sxy += x[i] * y[i]; }
	double denom = (double)n * sxx - sx * sx;
	if (std::abs(denom) < 1e-9) return 0.0;
	return ((double)n * sxy - sx * sy) / denom;
}

double vph_slope_vs_vertex_dist(std::vector<int>& pl, std::vector<double>& vph, int vertex_pl) {
	std::vector<double> dist(pl.size());
	for (size_t i = 0; i < pl.size(); i++) dist[i] = std::abs(pl[i] - vertex_pl);
	return linreg_slope(dist, vph);
}

double linklet_angle_stddev(std::vector<int>& pl, std::vector<double>& ax, std::vector<double>& ay) {
	size_t n = pl.size();
	if (n < 3) return -1.0; // need >=2 linklet angles for a meaningful stddev
	std::vector<size_t> order(n);
	for (size_t i = 0; i < n; i++) order[i] = i;
	std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return pl[a] < pl[b]; });

	std::vector<double> angles;
	for (size_t i = 0; i + 1 < n; i++) {
		size_t a = order[i], b = order[i + 1];
		matrix_3D::vector_3D d0, d1;
		d0.x = ax[a]; d0.y = ay[a]; d0.z = 1;
		d1.x = ax[b]; d1.y = ay[b]; d1.z = 1;
		angles.push_back(matrix_3D::opening_angle(d0, d1));
	}
	double mean = 0;
	for (double v : angles) mean += v;
	mean /= angles.size();
	double var = 0;
	for (double v : angles) var += (v - mean) * (v - mean);
	var /= angles.size();
	return std::sqrt(var);
}

int main(int argc, char** argv) {
	if (argc != 5 && argc != 6) {
		fprintf(stderr, "usage : prg m_ev_group.all momch out.txt ecc_path [water_pene.txt]\n");
		fprintf(stderr, "  ecc_path : e.g. I:\\NINJA\\E71a\\ECC5 (holds Area0\\0\\align\\fine\\ali_PPP_QQQ_interpolation.txt)\n");
		exit(1);
	}
	std::string file_group = argv[1];
	std::string file_momch = argv[2];
	std::string file_out = argv[3];
	g_ecc_path = argv[4];
	std::string file_water_pene = (argc == 6) ? argv[5] : "";

	fprintf(stderr, "reading group mfile...\n");
	mfile0::Mfile group;
	mfile1::read_mfile_extension(file_group, group);

	fprintf(stderr, "reading momch...\n");
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_momch);

	bool has_water_pene = !file_water_pene.empty();
	std::map<Lst, int> water_pene = read_water_pene(file_water_pene);

	// per-event lookups
	std::map<int, Momentum_recon::Event_information*> ev_by_gid;
	std::map<Lst, Momentum_recon::Mom_chain*> chain_by_id;
	std::map<int, std::map<PlRaw, Momentum_recon::Mom_chain*>> rawid_index_by_gid;
	std::map<int, Momentum_recon::Mom_chain*> muon_by_gid;

	for (auto& ev : momch) {
		ev_by_gid[ev.groupid] = &ev;
		auto& idx = rawid_index_by_gid[ev.groupid];
		for (auto& c : ev.chains) {
			chain_by_id[Lst{ ev.groupid, c.chainid }] = &c;
			if (is_muon(c) && muon_by_gid.count(ev.groupid) == 0) muon_by_gid[ev.groupid] = &c;
			for (auto& b : c.base) {
				idx[PlRaw{ b.pl, (int64_t)b.rawid }] = &c;
			}
		}
	}

	// group the flat chain list from m_ev_group.all by composite id (basetracks[0].group_id),
	// preserving file order: for a "...1" id, element0=muon-copy, element1=partner-copy, element2+=candidates
	std::ofstream ofs(file_out);
	std::string stem = file_out;
	{ size_t dot = stem.rfind('.'); if (dot != std::string::npos) stem = stem.substr(0, dot); }
	std::ofstream ofs_radlat(stem + "_self_radlat.txt");
	std::ofstream ofs_vph(stem + "_self_vph.txt");
	std::ofstream ofs_scale(stem + "_chain_scale.txt");
	// ONE-OFF diagnostic (2026-11-12): for each partner chain's own internal
	// consecutive basetrack pairs (same chain, no vertex crossing), the raw
	// dal/dar/lat_pos/rad_pos via the exact same lateral_radial_diff() calls
	// used for the cross-chain junction -- tests whether the systematic
	// negative skew seen in radial_pos_diff at the partner-candidate junction
	// is a sign/convention artifact of the formula itself, or specific to
	// the junction construction. groupid chainid pl0 pl1 dal dar lat_pos rad_pos
	std::ofstream ofs_selfradial(stem + "_selfradial.txt");
	// groupid chainid cand_idx  p_pl p_x p_y p_z p_ax p_ay  cand_pl cand_x
	// cand_y cand_z cand_ax cand_ay  pc_md -- raw inputs to the pc_md
	// minimum_distance_fixed() call, one row per partner x candidate pair,
	// for manual sanity-checking of the md calculation.
	std::ofstream ofs_mddebug(stem + "_md_debug.txt");
	// groupid chainid cand_idx cand_chainid wp_tag valid  p_pl cand_pl
	// pl_diff  dal_boundary dar_boundary lateral_pos_diff radial_pos_diff
	// pc_md_boundary  p_vph_mean_chain cand_vph_mean_chain vph_mean_diff
	// partner_mcsmom cand_nseg cand_npl (20 columns, 0-indexed 0..19).
	// partner_mcsmom is the partner's own MCS momentum (ecc_mcs_mom[muon or
	// proton idx]), used to test whether lower-momentum partners show larger
	// MCS-driven scattering at the junction (Highland: theta_MCS ~ 1/p).
	// cand_nseg is the candidate chain's own nseg field (number of
	// basetracks); cand_npl is (max pl - min pl + 1) spanned by those
	// basetracks -- fillfactor = cand_nseg / cand_npl measures how complete
	// (gap-free) the candidate's own reconstruction is.
	// -- the topology-aware partner-candidate junction: partner's point is
	// whichever end of the partner chain is closest to vertex_pl (largest
	// pl<=vertex_pl if it has one, else smallest pl overall); candidate's
	// point is searched on the OPPOSITE side of vertex_pl (smallest
	// pl>=vertex_pl+1, or largest pl<=vertex_pl). valid=0 (and -1 sentinels)
	// if the candidate chain has no basetrack on that required side.
	std::ofstream ofs_boundary(stem + "_boundary.txt");
	int n_boundary_valid = 0, n_boundary_invalid = 0;
	int n_groups_with_cand = 0, n_pairs = 0, n_matched = 0;

	int64_t cur_id = -1;
	std::vector<mfile0::M_Chain*> cur;
	auto flush_group = [&]() {
		if (cur.size() < 3) { cur.clear(); return; } // need muon + partner + >=1 candidate
		if (cur_id % 10 != 1) { cur.clear(); return; } // only the "...1" (candidate-overlay) variant carries candidates

		int64_t id = cur_id;
		int groupid = (int)(id / 100000);
		int chainid = (int)((id % 100000) / 10);

		Lst pk{ groupid, chainid };
		auto pc_itr = chain_by_id.find(pk);
		auto ev_itr = ev_by_gid.find(groupid);
		auto mu_itr = muon_by_gid.find(groupid);
		if (pc_itr == chain_by_id.end() || ev_itr == ev_by_gid.end() || mu_itr == muon_by_gid.end()) { cur.clear(); return; }

		Momentum_recon::Mom_chain& partner_mc = *pc_itr->second;
		Momentum_recon::Event_information& ev = *ev_itr->second;
		Momentum_recon::Mom_chain& muon_mc = *mu_itr->second;
		if (partner_mc.base.empty() || muon_mc.base.empty()) { cur.clear(); return; }

		double mu_ax, mu_ay, mu_x, mu_y, mu_z, mu_vph; int mu_pl;
		get_vertex_side_point(muon_mc, ev.vertex_pl, mu_ax, mu_ay, mu_x, mu_y, mu_z, mu_vph, mu_pl);
		double zrange[2]; zrange_for_vertex(ev.vertex_pl, zrange);

		double p_ax, p_ay, p_x, p_y, p_z, p_vph; int p_pl;
		get_vertex_side_point(partner_mc, ev.vertex_pl, p_ax, p_ay, p_x, p_y, p_z, p_vph, p_pl);

		matrix_3D::vector_3D mu_pos, mu_dir, p_pos, p_dir;
		mu_pos.x = mu_x; mu_pos.y = mu_y; mu_pos.z = mu_z;
		mu_dir.x = mu_ax; mu_dir.y = mu_ay; mu_dir.z = 1;
		p_pos.x = p_x; p_pos.y = p_y; p_pos.z = p_z;
		p_dir.x = p_ax; p_dir.y = p_ay; p_dir.z = 1;
		double extra[2]; double zr[2] = { zrange[0], zrange[1] };
		double partner_md = minimum_distance_fixed(mu_pos, p_pos, mu_dir, p_dir, zr, extra, mu_z);
		double partner_oa = matrix_3D::opening_angle(mu_dir, p_dir);

		int p_mom_idx = (partner_mc.particle_flg == 2212) ? 1 : 0;
		double partner_mcsmom = partner_mc.ecc_mcs_mom[p_mom_idx];
		double partner_rangemom = partner_mc.ecc_range_mom[p_mom_idx];

		std::vector<int> p_pl_list;
		std::vector<double> p_vph_list, p_ax_list, p_ay_list;
		extract_chain_info(partner_mc, p_pl_list, p_vph_list, p_ax_list, p_ay_list);
		double partner_vph_slope = vph_slope_vs_vertex_dist(p_pl_list, p_vph_list, ev.vertex_pl);
		double partner_linklet_std = linklet_angle_stddev(p_pl_list, p_ax_list, p_ay_list);

		{
			std::vector<Momentum_recon::Mom_basetrack> sorted_base = partner_mc.base;
			std::sort(sorted_base.begin(), sorted_base.end(), [](const Momentum_recon::Mom_basetrack& a, const Momentum_recon::Mom_basetrack& b) { return a.pl < b.pl; });
			for (size_t k = 0; k + 1 < sorted_base.size(); k++) {
				Momentum_recon::Mom_basetrack& b0 = sorted_base[k];
				Momentum_recon::Mom_basetrack& b1 = sorted_base[k + 1];
				double dal_s, dar_s;
				lateral_radial_diff((double)b0.ax, (double)b0.ay, (double)b1.ax, (double)b1.ay, dal_s, dar_s);
				double dx_s = (double)b1.x - (double)b0.x, dy_s = (double)b1.y - (double)b0.y;
				double dz_s = (double)b1.z - (double)b0.z;
				double pred_dx_s = (double)b0.ax * dz_s, pred_dy_s = (double)b0.ay * dz_s;
				double resid_dx_s = dx_s - pred_dx_s, resid_dy_s = dy_s - pred_dy_s;
				double lat_s, rad_s;
				lateral_radial_diff((double)b0.ax, (double)b0.ay, (double)b0.ax + resid_dx_s, (double)b0.ay + resid_dy_s, lat_s, rad_s);
				ofs_selfradial << groupid << " " << chainid << " " << b0.pl << " " << b1.pl << " "
					<< std::fixed << std::setprecision(5) << dal_s << " " << dar_s << " "
					<< std::setprecision(2) << lat_s << " " << rad_s << std::endl;
			}
		}

		double rms_l, rms_r;
		chain_linklet_lateral_radial_rms(p_pl_list, p_ax_list, p_ay_list, rms_l, rms_r);
		double p_vph_mean, p_vph_sd;
		vec_mean_sd(p_vph_list, p_vph_mean, p_vph_sd);

		int wp_tag = -1;
		if (has_water_pene) {
			auto itr = water_pene.find(pk);
			wp_tag = (itr != water_pene.end()) ? itr->second : 0;
		}

		// partner's boundary point for the topology-aware junction (see
		// partner_boundary_point/candidate_boundary_point above), computed
		// once per partner -- independent of p_pos/p_dir/p_pl above, which
		// keep the old vertex_pl-nearest convention for the muon comparison.
		bool pb_upstream; int pb_pl; double pb_x, pb_y, pb_z, pb_ax, pb_ay;
		partner_boundary_point(partner_mc, ev.vertex_pl, pb_upstream, pb_pl, pb_x, pb_y, pb_z, pb_ax, pb_ay);

		n_groups_with_cand++;
		auto& rawid_index = rawid_index_by_gid[groupid];

		// partner-self chi2 distributions (role=0), written once per partner
		std::vector<double> p_chi2_al_self, p_chi2_ar_self, p_chi2_vph_self;
		chain_self_radlat_chi2(p_pl_list, p_ax_list, p_ay_list, p_chi2_al_self, p_chi2_ar_self);
		chain_self_vph_chi2(p_vph_list, p_chi2_vph_self);
		for (size_t k = 0; k < p_chi2_al_self.size(); k++) {
			ofs_radlat << groupid << " " << chainid << " " << 0 << " " << -1 << " " << -1 << " "
				<< wp_tag << " " << (int)k << " "
				<< std::fixed << std::setprecision(2) << p_chi2_al_self[k] << " " << p_chi2_ar_self[k] << std::endl;
		}
		for (size_t k = 0; k < p_chi2_vph_self.size(); k++) {
			ofs_vph << groupid << " " << chainid << " " << 0 << " " << -1 << " " << -1 << " "
				<< wp_tag << " " << (int)k << " "
				<< std::fixed << std::setprecision(2) << p_chi2_vph_self[k] << std::endl;
		}

		// raw per-chain scattering scale (role=0, partner): rms_l/rms_r and
		// vph mean/sd computed above from ALL of the partner's own
		// basetracks. Unlike the chi2 columns above, these are NOT
		// self-normalized, so they are the quantities actually worth
		// comparing partner-vs-candidate (a chi2 built from a chain's own
		// rms, then averaged back over that same chain, is tautologically
		// ~1 per chi2 regardless of track quality -- the raw scale is not).
		ofs_scale << groupid << " " << chainid << " " << 0 << " " << -1 << " " << -1 << " "
			<< wp_tag << " " << (int)p_pl_list.size() << " "
			<< std::fixed << std::setprecision(5) << rms_l << " " << rms_r << " "
			<< std::setprecision(2) << p_vph_mean << " " << p_vph_sd << std::endl;

		for (size_t ci = 2; ci < cur.size(); ci++) {
			mfile0::M_Chain& cand = *cur[ci];
			mfile0::M_Base* nb = nearest_base(cand, ev.vertex_pl);
			if (nb == nullptr) continue;
			n_pairs++;

			double cand_vph = (double)(((int64_t)nb->ph) % 10000);
			int cand_pl = nb->pos / 10;

			matrix_3D::vector_3D c_pos, c_dir;
			c_pos.x = nb->x; c_pos.y = nb->y; c_pos.z = nb->z;
			c_dir.x = (double)nb->ax; c_dir.y = (double)nb->ay; c_dir.z = 1;
			double pos_diff = matrix_3D::distance(p_pos, c_pos);
			double angle_diff = matrix_3D::opening_angle(p_dir, c_dir);

			double extra2[2]; double zr2[2] = { zrange[0], zrange[1] };
			double cand_md = minimum_distance_fixed(mu_pos, c_pos, mu_dir, c_dir, zr2, extra2, mu_z);

			double extra3[2]; double zr3[2] = { zrange[0], zrange[1] };
			double pc_md = minimum_distance_fixed(p_pos, c_pos, p_dir, c_dir, zr3, extra3, mu_z);

			// debug dump of the raw basetrack-level inputs to pc_md above,
			// for manual inspection: the plate/position/direction of the
			// single basetrack used on each side (partner's vertex-side
			// basetrack via get_vertex_side_point; candidate's vertex-side
			// basetrack via nearest_base), plus the resulting pc_md.
			ofs_mddebug << groupid << " " << chainid << " " << (int)(ci - 2) << " "
				<< p_pl << " " << std::fixed << std::setprecision(2) << p_x << " " << p_y << " " << p_z << " "
				<< std::setprecision(5) << p_ax << " " << p_ay << " "
				<< cand_pl << " " << std::setprecision(2) << nb->x << " " << nb->y << " " << nb->z << " "
				<< std::setprecision(5) << (double)nb->ax << " " << (double)nb->ay << " "
				<< std::setprecision(2) << pc_md << std::endl;

			std::vector<int> c_pl_list;
			std::vector<double> c_vph_list, c_ax_list, c_ay_list;
			extract_chain_info(cand, c_pl_list, c_vph_list, c_ax_list, c_ay_list);
			double cand_vph_slope = vph_slope_vs_vertex_dist(c_pl_list, c_vph_list, ev.vertex_pl);
			double cand_linklet_std = linklet_angle_stddev(c_pl_list, c_ax_list, c_ay_list);

			// topology-aware junction (see ofs_boundary header comment above)
			{
				int cb_pl; double cb_x, cb_y, cb_z, cb_ax, cb_ay;
				bool cb_found = candidate_boundary_point(cand, ev.vertex_pl, pb_upstream, cb_pl, cb_x, cb_y, cb_z, cb_ax, cb_ay);
				int b_cand_chainid = -1;
				{
					auto bmi = rawid_index.find(PlRaw{ cand_pl, (int64_t)nb->rawid });
					if (bmi != rawid_index.end()) b_cand_chainid = bmi->second->chainid;
				}
				int cand_npl = 0;
				if (!c_pl_list.empty()) {
					int plmin = c_pl_list[0], plmax = c_pl_list[0];
					for (size_t k = 0; k < c_pl_list.size(); k++) { if (c_pl_list[k] < plmin) plmin = c_pl_list[k]; if (c_pl_list[k] > plmax) plmax = c_pl_list[k]; }
					cand_npl = plmax - plmin + 1;
				}
				int cand_nseg = cand.nseg;

				if (!cb_found) {
					n_boundary_invalid++;
					ofs_boundary << groupid << " " << chainid << " " << (int)(ci - 2) << " " << b_cand_chainid << " "
						<< wp_tag << " " << 0 << " "
						<< pb_pl << " " << -1 << " " << -1 << " "
						<< -1 << " " << -1 << " " << -1 << " " << -1 << " " << -1 << " "
						<< std::fixed << std::setprecision(2) << p_vph_mean << " " << -1 << " " << -1 << " "
						<< std::setprecision(1) << partner_mcsmom << " " << cand_nseg << " " << cand_npl << std::endl;
				}
				else {
					int pl_diff = cb_pl - pb_pl;

					// local/linklet-frame junction quantities (see
					// junction_local_quantities header comment above) --
					// replaces the earlier raw-mfile-frame dal/lat_pos/md.
					double dal_b, dar_b, lat_pos_b, rad_pos_b, pcb_md;
					bool jl_ok = junction_local_quantities(pb_pl, pb_x, pb_y, pb_z, pb_ax, pb_ay,
						cb_pl, cb_x, cb_y, cb_z, cb_ax, cb_ay, mu_z, zrange, dal_b, dar_b, lat_pos_b, rad_pos_b, pcb_md);
					if (!jl_ok) {
						n_boundary_invalid++;
						ofs_boundary << groupid << " " << chainid << " " << (int)(ci - 2) << " " << b_cand_chainid << " "
							<< wp_tag << " " << 0 << " "
							<< pb_pl << " " << cb_pl << " " << pl_diff << " "
							<< -1 << " " << -1 << " " << -1 << " " << -1 << " " << -1 << " "
							<< std::fixed << std::setprecision(2) << p_vph_mean << " " << -1 << " " << -1 << " "
							<< std::setprecision(1) << partner_mcsmom << " " << cand_nseg << " " << cand_npl << std::endl;
					}
					else {
						n_boundary_valid++;
						double cb_vph_mean, cb_vph_sd;
						vec_mean_sd(c_vph_list, cb_vph_mean, cb_vph_sd);
						double vph_mean_diff = std::abs(p_vph_mean - cb_vph_mean);

						ofs_boundary << groupid << " " << chainid << " " << (int)(ci - 2) << " " << b_cand_chainid << " "
							<< wp_tag << " " << 1 << " "
							<< pb_pl << " " << cb_pl << " " << pl_diff << " "
							<< std::fixed << std::setprecision(5) << dal_b << " " << dar_b << " "
							<< std::setprecision(2) << lat_pos_b << " " << rad_pos_b << " " << pcb_md << " "
							<< std::setprecision(2) << p_vph_mean << " " << cb_vph_mean << " " << vph_mean_diff << " "
							<< std::setprecision(1) << partner_mcsmom << " " << cand_nseg << " " << cand_npl << std::endl;
					}
				}
			}

			double dal, dar;
			lateral_radial_diff(p_ax, p_ay, (double)nb->ax, (double)nb->ay, dal, dar);
			double chi2_al = std::isfinite(rms_l) && rms_l > 0 ? (dal / rms_l) * (dal / rms_l) : -1;
			double chi2_ar = std::isfinite(rms_r) && rms_r > 0 ? (dar / rms_r) * (dar / rms_r) : -1;
			double sigma_vph = (p_vph_sd > 0 && cand_vph > 0) ? (cand_vph - p_vph_mean) / p_vph_sd : -1;
			double chi2_vph = (sigma_vph != -1) ? sigma_vph * sigma_vph : -1;

			// try to recover the candidate's true momch identity via (pl,rawid), trying the
			// nearest-vertex basetrack first, then falling back to any basetrack in the candidate
			Momentum_recon::Mom_chain* matched_chain = nullptr;
			{
				auto mi = rawid_index.find(PlRaw{ nb->pos / 10, (int64_t)nb->rawid });
				if (mi != rawid_index.end()) matched_chain = mi->second;
			}
			if (matched_chain == nullptr) {
				for (auto& b : cand.basetracks) {
					auto mi = rawid_index.find(PlRaw{ b.pos / 10, (int64_t)b.rawid });
					if (mi != rawid_index.end()) { matched_chain = mi->second; break; }
				}
			}

			int cand_matched = (matched_chain != nullptr) ? 1 : 0;
			int cand_chainid = cand_matched ? matched_chain->chainid : -1;
			double cand_mcsmom = -1, cand_rangemom = -1, mom_diff = -1;
			if (cand_matched) {
				n_matched++;
				int c_mom_idx = (matched_chain->particle_flg == 2212) ? 1 : 0;
				cand_mcsmom = matched_chain->ecc_mcs_mom[c_mom_idx];
				cand_rangemom = matched_chain->ecc_range_mom[c_mom_idx];
				mom_diff = std::abs(partner_mcsmom - cand_mcsmom);
			}

			// candidate-self chi2 distributions (role=1), skipping candidates
			// that turn out to be the partner's own chain re-listed (same
			// "self-match" condition the per-partner plot uses)
			if (!(cand_matched && cand_chainid == chainid)) {
				std::vector<double> c_chi2_al_self, c_chi2_ar_self, c_chi2_vph_self;
				chain_self_radlat_chi2(c_pl_list, c_ax_list, c_ay_list, c_chi2_al_self, c_chi2_ar_self);
				chain_self_vph_chi2(c_vph_list, c_chi2_vph_self);
				int cidx = (int)(ci - 2);
				for (size_t k = 0; k < c_chi2_al_self.size(); k++) {
					ofs_radlat << groupid << " " << chainid << " " << 1 << " " << cidx << " " << cand_chainid << " "
						<< wp_tag << " " << (int)k << " "
						<< std::fixed << std::setprecision(2) << c_chi2_al_self[k] << " " << c_chi2_ar_self[k] << std::endl;
				}
				for (size_t k = 0; k < c_chi2_vph_self.size(); k++) {
					ofs_vph << groupid << " " << chainid << " " << 1 << " " << cidx << " " << cand_chainid << " "
						<< wp_tag << " " << (int)k << " "
						<< std::fixed << std::setprecision(2) << c_chi2_vph_self[k] << std::endl;
				}

				// raw per-chain scattering scale (role=1, this one candidate),
				// same quantities as the partner row above -- the ones
				// actually meaningful to compare partner-vs-candidate.
				double c_rms_l, c_rms_r;
				chain_linklet_lateral_radial_rms(c_pl_list, c_ax_list, c_ay_list, c_rms_l, c_rms_r);
				double c_vph_mean, c_vph_sd;
				vec_mean_sd(c_vph_list, c_vph_mean, c_vph_sd);
				ofs_scale << groupid << " " << chainid << " " << 1 << " " << cidx << " " << cand_chainid << " "
					<< wp_tag << " " << (int)c_pl_list.size() << " "
					<< std::fixed << std::setprecision(5) << c_rms_l << " " << c_rms_r << " "
					<< std::setprecision(2) << c_vph_mean << " " << c_vph_sd << std::endl;
			}

			ofs << std::right << std::fixed
				<< std::setw(6) << groupid << " "
				<< std::setw(5) << chainid << " "
				<< std::setw(7) << std::setprecision(0) << p_vph << " "
				<< std::setw(9) << std::setprecision(1) << partner_mcsmom << " "
				<< std::setw(9) << std::setprecision(1) << partner_rangemom << " "
				<< std::setw(9) << std::setprecision(2) << partner_md << " "
				<< std::setw(7) << std::setprecision(4) << partner_oa << " "
				<< std::setw(3) << (int)(ci - 2) << " "
				<< std::setw(7) << std::setprecision(0) << cand_vph << " "
				<< std::setw(4) << cand_pl << " "
				<< std::setw(10) << std::setprecision(1) << pos_diff << " "
				<< std::setw(7) << std::setprecision(4) << angle_diff << " "
				<< std::setw(9) << std::setprecision(2) << cand_md << " "
				<< std::setw(2) << cand_matched << " "
				<< std::setw(6) << cand_chainid << " "
				<< std::setw(9) << std::setprecision(1) << cand_mcsmom << " "
				<< std::setw(9) << std::setprecision(1) << cand_rangemom << " "
				<< std::setw(9) << std::setprecision(1) << mom_diff << " "
				<< std::setw(2) << wp_tag << " "
				<< std::setw(9) << std::setprecision(2) << pc_md << " "
				<< std::setw(9) << std::setprecision(3) << partner_vph_slope << " "
				<< std::setw(9) << std::setprecision(3) << cand_vph_slope << " "
				<< std::setw(9) << std::setprecision(4) << partner_linklet_std << " "
				<< std::setw(9) << std::setprecision(4) << cand_linklet_std << " "
				<< std::setw(9) << std::setprecision(5) << dal << " "
				<< std::setw(9) << std::setprecision(5) << dar << " "
				<< std::setw(9) << std::setprecision(5) << rms_l << " "
				<< std::setw(9) << std::setprecision(5) << rms_r << " "
				<< std::setw(10) << std::setprecision(2) << chi2_al << " "
				<< std::setw(10) << std::setprecision(2) << chi2_ar << " "
				<< std::setw(8) << std::setprecision(3) << sigma_vph << " "
				<< std::setw(10) << std::setprecision(2) << chi2_vph << std::endl;
		}
		cur.clear();
	};

	for (auto& c : group.chains) {
		if (c.basetracks.empty()) continue;
		int64_t id = c.basetracks[0].group_id;
		if (id != cur_id) {
			flush_group();
			cur_id = id;
		}
		cur.push_back(&c);
	}
	flush_group();

	std::cout << "* #of partner groups with >=1 attach candidate = " << n_groups_with_cand << std::endl;
	std::cout << "* #of partner x candidate pairs                = " << n_pairs << std::endl;
	std::cout << "* #of candidates matched back to a momch chain = " << n_matched << std::endl;
	std::cout << "* #of boundary junctions valid / invalid       = " << n_boundary_valid << " / " << n_boundary_invalid << std::endl;
}
