/*
 * File: document_writer.cpp
 *
 * Copyright (C) 2016 Richard Eliáš <richard.elias@matfyz.cz>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
 * USA.
 */


#include "document_writer.hpp"
#include "svg_writer.hpp"
#include "ps_writer.hpp"
#include "traveler_writer.hpp"

using namespace std;

#define COLOR_DELETE        RGB::GRAY
#define COLOR_INSERT        RGB::RED
#define COLOR_REINSERT      RGB::BLUE
#define COLOR_EDIT          RGB::GREEN
#define COLOR_ROTATE        RGB::BROWN
#define COLOR_DEFAULT       RGB::BLACK

// initialize RGB constants:
const RGB RGB::RED = RGB(1., 0., 0., "red");
const RGB RGB::GREEN = RGB(0., 1., 0., "green");
const RGB RGB::BLUE = RGB(0., 0., 1., "blue");
const RGB RGB::BLACK = RGB(0., 0., 0., "black");
const RGB RGB::GRAY = RGB(0.8, 0.8, 0.8, "gray");
const RGB RGB::BROWN = RGB(0.83, 0.41, 0.12, "brown");

RGB::RGB(
         double _red,
         double _green,
         double _blue,
         const std::string& _name)
: red(_red), green(_green), blue(_blue), name(_name)
{ }

bool RGB::operator==(
                     const RGB& other) const
{
    return name == other.name;
}

struct pseudoknot_segment{
    pair<rna_tree::pre_post_order_iterator, rna_tree::pre_post_order_iterator> interval1;
    pair<rna_tree::pre_post_order_iterator, rna_tree::pre_post_order_iterator> interval2;


//    void set_interval1(int begin, int end){
//        interval1 = make_pair(begin, end);
//    }
//
//    void set_interval2(int begin, int end){
//        interval2 = make_pair(begin, end);
//    }
//
//    pair<int, int> get_interval1(){
//        return interval1;
//    }
//
//    pair<int, int> get_interval2(){
//        return interval2;
//    }
};



/* static */ image_writers document_writer::get_writers(
                                                        bool use_colors)
{
    image_writers vec;
    vec.emplace_back(new svg_writer());
    vec.emplace_back(new ps_writer());
    vec.emplace_back(new traveler_writer());
    
    for (const auto& writer : vec)
        writer->use_colors(use_colors);
    
    return vec;
}

/* static */ unique_ptr<document_writer> document_writer::get_traveler_writer()
{
    return unique_ptr<document_writer>(new traveler_writer());
}

std::string document_writer::get_edge_formatted(
                                                point from,
                                                point to,
                                                bool is_base_pair) const
{
    if (from.bad() || to.bad())
    {
        WARN("Cannot draw line between bad points");
        return "";
    }
    if (is_base_pair)
    {
        point tmp = rna_tree::base_pair_edge_point(from, to);
        to = rna_tree::base_pair_edge_point(to, from);
        from = tmp;
    }
    
    return get_line_formatted(from, to, RGB::BLACK);
}

std::string document_writer::get_label_formatted(
                                                 rna_tree::pre_post_order_iterator it, label_info li) const
{
    if (!it->initiated_points())
        return "";
    
    ostringstream out;
    
    out
    << get_label_formatted(it->at(it.label_index()), get_default_color(it->status), li);
    
    if (it->paired() &&
        it.preorder() &&
        it->initiated_points() &&
        !rna_tree::is_root(it))
    {
        out
        << get_edge_formatted(it->at(0).p, it->at(1).p, true);
    }

    auto x = out.str();
    return out.str();
}

const RGB& document_writer::get_default_color(
                                              rna_pair_label::status_type status) const
{
    if (!colored)
        return RGB::BLACK;
    
    switch (status)
    {
#define switchcase(status, rgb) \
case rna_pair_label::status: \
return COLOR_ ## rgb;
            
            switchcase(deleted, DELETE);
            switchcase(edited, EDIT);
            switchcase(inserted, INSERT);
            switchcase(reinserted, REINSERT);
            switchcase(rotated, ROTATE);
            switchcase(touched, DEFAULT);
            switchcase(untouched, DEFAULT);
        default:
            abort();
#undef switchcase
    }
}

void document_writer::print_to_stream(
                                      const std::string& text)
{
    out << text;
    validate_stream();
}

void document_writer::seek_from_current_pos(
                                            off_type offset)
{
    out.seekp(offset, fstream::cur);
    validate_stream();
}

void document_writer::validate_stream() const
{
    if (out.fail())
        throw io_exception("Writing document failed");
}

std::string document_writer::get_rna_subtree_formatted(
                                                       rna_tree::iterator root) const
{
    ostringstream out;

    int seq_ix = 0;
    auto print =
    [&out, &seq_ix, this](rna_tree::pre_post_order_iterator it)
    {
        out << get_label_formatted(it, {seq_ix, it->at(it.label_index()).tmp_label});
        seq_ix++;
    };
    
    rna_tree::for_each_in_subtree(root, print);
    
    return out.str();
}

double document_writer::get_scaling_ratio() const{
    return scaling_ratio;
}

void document_writer::set_scaling_ratio(rna_tree& rna){
    auto bp_dist = rna.get_pair_base_distance();
    scaling_ratio = 20 / bp_dist;
};

std::string document_writer::get_rna_background_formatted(
                                                          rna_tree::pre_post_order_iterator begin,
                                                          rna_tree::pre_post_order_iterator end) const
{
    rna_tree::pre_post_order_iterator prev;
    ostringstream out;
    
    while (++rna_tree::pre_post_order_iterator(begin) != end)
    {
        prev = begin++;
        
        point p1 = prev->at(prev.label_index()).p;
        point p2 = begin->at(begin.label_index()).p;
        
        if (p1.bad() || p2.bad())
            continue;

        point diff_orig = p2 - p1;

        point tmp = rna_tree::base_pair_edge_point(p1, p2, get_scaling_ratio());
        p2 = rna_tree::base_pair_edge_point(p2, p1, get_scaling_ratio());
        p1 = tmp;

        point diff_edge = p2 - p1;

        point diff = diff_orig * diff_edge;

        //If the edge points cross, then the line should not be drawn at all
        if (diff.x > 0 && diff.y > 0) out << get_line_formatted(p1, p2, RGB::GRAY);
    }
    
    return out.str();
}

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

std::string document_writer::find_pseudoknots(rna_tree::pre_post_order_iterator begin, rna_tree::pre_post_order_iterator end) const
{
    ostringstream out;

    auto pn_segments = find_pseudoknot_segments(begin, end);
    for (auto s:pn_segments){

        auto l = s.interval1.first->at(s.interval1.first.label_index());
        auto ll = s.interval2.first->at(s.interval2.first.label_index());

        out << get_line_formatted(l.p, ll.p, RGB::RED);
    }
    

    return out.str();
}

std::string document_writer::get_rna_formatted(
                                               rna_tree rna) const
{
    return get_rna_subtree_formatted(rna.begin())
    + get_rna_background_formatted(rna.begin_pre_post(), rna.end_pre_post())
    + find_pseudoknots(rna.begin_pre_post(), rna.end_pre_post());
}

void document_writer::init(
                           const std::string& filename,
                           const std::string& suffix)
{
    APP_DEBUG_FNAME;
    assert(!filename.empty());
    
    string file = filename + suffix;
    INFO("Opening document %s for writing RNA", file);
    
    out.close();
    
    // create file & truncate
    out.open(file, ios::out);
    out.close();
    
    // open in normal mode
    out.open(file, ios::out | ios::in);
    out
    << std::unitbuf
    << std::scientific;
    
    if (!out.good())
        throw io_exception("Cannot open output file %s for writing.", filename);
    assert(out.good());
}

void document_writer::seek(
                           streampos pos)
{
    out.seekp(pos);
    
    assert(out.good());
}

void document_writer::seek_end()
{
    out.seekp(0, out.end);
    
    assert(out.good());
}

streampos document_writer::get_pos()
{
    streampos pos = out.tellp();
    assert(pos != -1);
    
    return pos;
}

size_t document_writer::fill(
                             char ch)
{
    streampos pos, end;
    
    pos = get_pos();
    seek_end();
    end = get_pos();
    
    size_t n = end - pos;
    
    if (n != 0)
    {
        seek(pos);
        out
        << string(n - 1, ch)
        << endl;
    }
    seek(pos);
    return n;
}

void document_writer::use_colors(
                                 bool _colored)
{
    colored = _colored;
}
