/*
 * File: geometry.hpp
 *
 * Copyright (C) 2019 David Hoksza <david.hoksza@gmail.com>
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

#ifndef TRAVELER_GEOMETRY_HPP
#define TRAVELER_GEOMETRY_HPP

#include "point.hpp"

int orientation(point p, point q, point r);
bool lines_intersect(point p1, point q1, point p2, point q2);
point lines_intersection(point p1, point q1, point p2, point q2);


#endif //TRAVELER_GEOMETRY_HPP
