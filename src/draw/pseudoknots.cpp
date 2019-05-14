//
// Created by David Hoksza on 14.05.19.
//

#include "pseudoknots.hpp"
#include "convex_hull.hpp"

using namespace std;

vector<pseudoknot_segment> find_pseudoknot_segments(rna_tree::pre_post_order_iterator begin, rna_tree::pre_post_order_iterator end){

    vector<pair<rna_tree::pre_post_order_iterator, rna_tree::pre_post_order_iterator>> pn_pairs;

    int ix_begin = 0;
    while(begin != end)
    {
        for (int i = 0; i < begin->size(); ++i)
//        for(auto&& l: begin->labels)
        {
            auto l = (*begin)[i];
            if(!l.pseudoknot.empty())
            {
                auto rest = begin;
                auto ix_rest = ix_begin;
                ++rest; ++ix_rest;

                while(rest != end)
                {
                    for (int j = 0; j < rest->size(); ++j)
//                    for(auto&& ll: rest->labels)
                    {
                        auto ll = (*rest)[j];
                        if(ll.pseudoknot == l.pseudoknot)
                        {
                            pn_pairs.push_back(make_pair(begin, rest));
                        }
                    }
                    ++rest;
                }
            }
        }

        ++begin; ++ix_begin;
    }

    vector<pseudoknot_segment> segments;


    if (pn_pairs.size() > 0){
        pseudoknot_segment s = {make_pair(pn_pairs[0].first, pn_pairs[0].first), make_pair(pn_pairs[0].second, pn_pairs[0].second)};
        for (int i = 1; i <= pn_pairs.size(); ++i){
            auto next1 = rna_tree::pre_post_order_iterator(s.interval1.second)++;
            auto next2 = rna_tree::pre_post_order_iterator(s.interval2.second)++;
            if (next1 == pn_pairs[i].first && next2 == pn_pairs[i].second){
                //if the first and second residue in the considered pseudoknot pair both directly extend the last pseudoknot segment, let's extend the segment
                s.interval1.second = pn_pairs[i].first;
                s.interval2.second = pn_pairs[i].second;
            } else {
                segments.push_back(s);
                s = {make_pair(pn_pairs[i].first, pn_pairs[i].first), make_pair(pn_pairs[i].second, pn_pairs[i].second)};
            }
        }
    }

    return segments;

}



vector<line> get_lines_from_points(vector<point> &points){
    vector<line> lines;
    for (int i = 1; i < points.size(); i++) {
        lines.emplace_back(points[i-1], points[i]);
//        lines.push_back(make_pair(points[i-1], points[i]));
    }
    if (points.size() > 1) {
        lines.emplace_back(points[points.size()-1], points[0]);
//        lines.push_back(make_pair(points[points.size()-1], points[0]));
    }
    return lines;
}

pair<int, point> get_closest_hull_intersection(vector<line> hull_lines, point p) {
    auto direction_lines = {
            make_pair(p, p + point(1,0)),
            make_pair(p, p + point(0,1))
    };

    double min_dist = numeric_limits<double>::max();
    point intersection_point = point::bad_point();
    int min_ix = -1;
    for (line dl:direction_lines){
        for (int ix = 0; ix< hull_lines.size(); ix++){
            line hl = hull_lines[ix];
            point i = lines_intersection(dl.first, dl.second, hl.first, hl.second);
            if (!i.bad()){
                auto dist = distance(p, i);
                if (dist < min_dist) {
                    min_dist = dist;
                    min_ix = ix;
                    intersection_point = i;
                }
            }
        }
    }
    assert(min_ix >= 0);
    return make_pair(min_ix, intersection_point);


}

vector<line> get_pseudoknot_curves(pseudoknot_segment pn, vector<point> hull){

    vector<line> curves;

    auto hull_lines = get_lines_from_points(hull);

    point begin = pn.interval1.first->at(pn.interval1.first.label_index()).p;
    point end = pn.interval2.first->at(pn.interval2.first.label_index()).p;

    auto hull_intersection_begin = get_closest_hull_intersection(hull_lines, begin);
    auto hull_intersection_end = get_closest_hull_intersection(hull_lines, end);

    int ix_begin = hull_intersection_begin.first;
    point intersection_begin = hull_intersection_begin.second;
    curves.emplace_back(begin, intersection_begin);

    int ix_end = hull_intersection_end.first;
    point intersection_end = hull_intersection_end.second;

    if (ix_begin == ix_end){
        curves.emplace_back(intersection_begin, intersection_end);

    } else {
        // We need to iterate through all lines from first to second intersection and compute the accumulated distance
        // Then we check whether this "clockwise" distance is lower then the anticlockwise distance and based on
        // that we add the lines to the curve

        vector<line> aux_curve;
        aux_curve.emplace_back(intersection_begin, hull_lines[ix_begin].second);
        int ix = ix_begin + 1;
        while (ix != ix_end ){
            aux_curve.emplace_back(hull_lines[ix].first, hull_lines[ix].second);
            ix = (ix + 1) % hull_lines.size();

        }
        aux_curve.emplace_back(hull_lines[ix_end].first, intersection_end);

        double hull_perimeter = 0;
        for (line hl: hull_lines) {
            hull_perimeter += distance(hl.first, hl.second);
        }

        double dist = 0;
        for (line c: aux_curve) {
            dist += distance(c.first, c.second);
        }

        if (dist <= hull_perimeter / 2) {
            curves.insert(curves.end(), aux_curve.begin(), aux_curve.end());
        } else {
            curves.emplace_back(intersection_begin, hull_lines[ix_begin].first);
            ix = ix_begin - 1;
            while (ix != ix_end ){
                curves.emplace_back(hull_lines[ix].second, hull_lines[ix].first);
                ix--;
                if (ix == -1) ix = hull_lines.size() - 1;
            }
            curves.emplace_back(hull_lines[ix_end].second, intersection_end);
        }
    }

    curves.emplace_back(intersection_end, end);

    return curves;
}

pseudoknots::pseudoknots(rna_tree &rna) {
    this->segments = find_pseudoknot_segments(rna.begin_pre_post(), rna.end_pre_post());

    auto points = rna.get_points();
    auto h = convex_hull(points);

    for (auto s:this->segments){
        add_padding(h, rna.get_pairs_distance());
        vector<line> lines = get_pseudoknot_curves(s, h);
        for (line l:lines){
            this->lines.push_back(l);
        }
    }

}



