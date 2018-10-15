/*
 * File: crw_extractor.cpp
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


#include <fstream>

#include "crw_extractor.hpp"
#include "types.hpp"


using namespace std;


void crw_extractor::extract(
                            const std::string& filename)
{
    APP_DEBUG_FNAME;
    
    labels.clear();
    points.clear();
    
    ifstream in(filename);
    
    regex regexp_base_line = create_regex(
                                          msprintf("^\\(%s\\)\\s+%s\\s+%s\\s+lwstring\\s*$",
                                                   BASE_REGEX, DOUBLE_REGEX, DOUBLE_REGEX));
    
    smatch match;
    string line;
    point p;
    string base;

    bool non_letter_line = false;
    
    while (true)
    {
        getline(in, line);
        if (in.fail())
            break;
        
        if (regex_search(line, match, regexp_base_line))
        {
            if (non_letter_line) {
                labels.clear();
                points.clear();
                non_letter_line = false;
            }

            // is base line
            stringstream str;
            
            str << match[1] << " " << match[2] << " " << match[4];
            
            str
            >> base
            >> p.x
            >> p.y;
            
            assert(!str.fail() && str.eof() && base.size() == 1);
            
            labels.push_back(base[0]);
            points.push_back(p);
        } else if (line.find("setfont") == std::string::npos){
            non_letter_line = true;
        }
    }
}

void crw_extractor::adjust_residues_lists(const int cnt){

    int sz = labels.size();

    assert(sz == points.size())

    if (cnt < sz) {
        labels.erase(0, sz - cnt);
        points.erase(points.begin(), points.begin() + sz - cnt);
    }
}
