/*
 * File: tests.cpp
 *
 * Copyright (C) 2015 Richard Eliáš <richard.elias@matfyz.cz>
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 */


#include "tests.hpp"
#include "rna_tree.hpp"
#include "utils.hpp"
#include "app.hpp"
#include "compact.hpp"

#include "write_ps_document.hpp"

#define INDIR           (string("precomputed/"))
#define OUTDIR_OP       (string("build/files/run-op/"))

#define FILES   (std::vector<string>({"1.hairpin", "2.interior", "3.multibranch", "4.fullbranch"}))

#define print(string)   psout.print((string), true)
#define parent(iter)    rna_tree::parent(iter)


#define FILEIN(index)  (INDIR + (FILES.at(index)))
#define FILEOUT(index) (OUTDIR_OP + (FILES.at(index)) + ".out")


static ps_writer psout;
rna_pair_label lbl_leaf("I");
rna_pair_label lbl_pair = lbl_leaf + lbl_leaf;

using namespace std;

#ifdef NODEF

void test::save_seq_fold_subtree(iterator it, string name)
{
    APP_DEBUG_FNAME;

#define OUTDIR_PART     (string("build/files/run-part/"))
    name = OUTDIR_PART + name;

    write_file(name + ".seq", rna_tree::get_labels(it));
    write_file(name + ".fold", rna_tree::get_brackets(it));
}

void test::generate()
{
    APP_DEBUG_FNAME;

    size_t n = 100;
    size_t i = 0;
    rna_tree rna = get_rna(INDIR + "frog");

    for (i = 0; i < 19; ++i)
    {
        iterator it = rna.begin();
        it = plusplus(it, n);

        while (rna_tree::is_leaf(it))
            ++it;

        string name = OUTDIR_PART + "rna_part." + to_string(i);
        save_seq_fold_subtree(it, name);
        save_to_psout(name + ".ps", it);

        n += 100;
    }
}

#endif

void test::save_seq_fold(rna_tree rna, std::string name)
{
    APP_DEBUG_FNAME;

    write_file(name + ".seq", rna.get_labels());
    write_file(name + ".fold", rna.get_brackets());
}





rna_tree get_rna(const string& name)
{
    DEBUG("get_rna(%s)", to_cstr(name));

    string l, b;
    l = read_file(name + ".seq");
    b = read_file(name + ".fold");

    ps_document doc(name + ".ps");

    return rna_tree(b, l, doc.points, name);
}

std::vector<std::string> test::create_app_arguments(const std::string& file_in, const std::string& file_out)
{
    vector<string> vec = {
        "program_name",
        "-tt",
            file_in + ".ps",
            file_in + ".fold",
            "--name", file_in + "_templ",
        "-mt",
            file_out + ".seq",
            file_out + ".fold",
            "--name", file_out + "_match",
        //"-r",
            //"--strategies", file_out + ".matched.rted",
        //"-g",
            //"--ted-out", file_out + ".matched.ted",
            //"--mapping", file_out + ".matched.map",
        //"--ps",
            //"--mapping", file_out + ".matched.map",
            //file_out + ".matched.ps",
        "-a",
            file_out + ".ps",
    };

    return vec;
}

string test::ending_3_5_strings(iterator it)
{
    if (rna_tree::is_root(it))
        return "";
    auto get_direction = [](iterator it) {
        assert(!rna_tree::is_leaf(it));

        point p1, p2, p, ch;

        p1 = it->at(0).p;
        p2 = it->at(1).p;

        ch = it.begin()->centre();

        p = -orthogonal(p2 - p1, ch - p1);

        return p;
    };

    point dir = get_direction(it) * BASES_DISTANCE * 1.5;
    point p1, p2;
    p1 = it->at(0).p;
    p2 = it->at(1).p;
    ostringstream out;

    out
        << ps_writer::sprint(p1 + dir, "5'")
        << ps_writer::sprint(p2 + dir, "3'")
        << ps_writer::sprint_edge(p1, p1 + dir, true)
        << ps_writer::sprint_edge(p2, p2 + dir, true);

    return out.str();
}

void test::save_to_psout(const std::string& filename, iterator it)
{
    APP_DEBUG_FNAME;

    psout.init_default(filename, it);

    print(psout.sprint_subtree(it));
    print(ending_3_5_strings(it));
}






void test::run()
{
    APP_DEBUG_FNAME;

    run_hairpin();
    run_interior();
}


void test::insert_hairpin(
                rna_tree& rna,
                sibling_iterator ch,
                size_t n)
{
    int i = 4;
    while (--i != 0)
        ch = rna.insert(ch, lbl_leaf, 0);
    ch = rna.insert(ch, lbl_pair, 3);

    while (--n != 0)
        ch = rna.insert(ch, lbl_pair, 1);
}



void test::run_app(std::string filetempl, std::string fileother)
{
    APP_DEBUG_FNAME;

    auto args = create_app_arguments(filetempl, fileother);

    app app;
    app.run(args);

    rna_tree rna = get_rna(fileother);
    save_to_psout(fileother + ".ps", rna.begin());
}

void test::run_hairpin()
{
    APP_DEBUG_FNAME;

#define FILEINDEX 0

    auto run_insert_bulge = [this](size_t n)
    {
        string filein = FILEIN(FILEINDEX);
        string fileout = FILEOUT(FILEINDEX) + ".bulge." + to_string(n);

        rna_tree rna = get_rna(filein);
        sibling_iterator it = plusplus(rna.begin(), 2);

        while (--n != 0)
            it = rna.insert(it, lbl_leaf, 0);

        rna.print_tree();
        WAIT;
        save_seq_fold(rna, fileout);
        run_app(filein, fileout);
    };

    run_insert_bulge(4);


    //iterator it;

    //it = plusplus(rna.begin(), 2);

    //rna.print_tree();

    //insert_hairpin(rna, it.begin(), 4);
    //insert_hairpin(rna, it.end(), 4);

    //rna.print_tree();

    //save_seq_fold(rna, FILEOUT(file_index));
    //run_app(FILEIN(file_index), FILEOUT(file_index));
}

void test::run_interior()
{
    APP_DEBUG_FNAME;

    size_t file_index = 1;
    rna_tree rna = get_rna(FILEIN(file_index));
    iterator it;
    sibling_iterator ch;

    it = plusplus(rna.begin(), 6);

    rna.print_tree();

    insert_hairpin(rna, plusplus(it.begin(), 2), 4);
    insert_hairpin(rna, plusplus(it.begin(), 5), 4);

    rna.print_tree();

    save_seq_fold(rna, FILEOUT(file_index));
    run_app(FILEIN(file_index), FILEOUT(file_index));
}

