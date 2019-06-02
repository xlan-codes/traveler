#ifndef TRAVELER_PSEUDOKNOTS_HPP

#include <vector>
#include "rna_tree.hpp"
#include "geometry.hpp"

typedef std::pair<point, point> line;
typedef std::vector<line> curve;


struct pseudoknot_segment{
    std::pair<rna_tree::pre_post_order_iterator, rna_tree::pre_post_order_iterator> interval1;
    std::pair<rna_tree::pre_post_order_iterator, rna_tree::pre_post_order_iterator> interval2;

    curve connecting_curve;

    std::string get_label() const;

};

struct pseudoknots {

    std::vector<pseudoknot_segment> segments;

    pseudoknots(rna_tree &rna);

};


#define TRAVELER_PSEUDOKNOTS_HPP

#endif //TRAVELER_PSEUDOKNOTS_HPP
