module maskerade.bed;

import intervaltree;
import intervaltree.iitree;
import intervaltree.cgranges: cr_init;
import std.stdio;
import std.array: split;
import std.conv: to;

IITree!BasicInterval getIITreefromBed(string fn){
    auto tree = IITree!BasicInterval(cr_init()); 
    auto file =  File(fn);
    string line;
    while ((line = file.readln()) !is null){
        auto fields = line[0..$-1].split("\t");
        tree.insert(fields[0],BasicInterval(fields[1].to!int,fields[2].to!int));
    }
    tree.index;
    return tree;
}

