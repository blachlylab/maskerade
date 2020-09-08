module maskerade.bam;
import dhtslib.sam;
import dhtslib.cigar;
import dhtslib.htslib.sam: bam_get_seq;
import intervaltree;
import intervaltree.iitree;
import std.algorithm: map, sum, filter;
import std.array: array;
import maskerade.nmask;
import std.stdio;

void processReads(SAMReader * bamr, SAMWriter * bamw, IITree!BasicInterval * tree, bool invert, bool eject){
    if(invert && eject) processReads!(true,true)(bamr,bamw,tree);
    else if(invert) processReads!(true,false)(bamr,bamw,tree);
    else if(eject) processReads!(false,true)(bamr,bamw,tree);
    else processReads!(false,false)(bamr,bamw,tree);
}

void processReads(bool invert, bool eject)(SAMReader * bamr, SAMWriter * bamw, IITree!BasicInterval * tree){
    auto contigs = bamr.target_names;
    foreach(SAMRecord rec; bamr.all_records){
        if(rec.isSecondary || !rec.isMapped) continue;
        auto contig = contigs[rec.tid];
        auto cigar = Cigar(rec.cigar.ops.dup);
        auto interval = BasicInterval(rec.pos, rec.pos + cigar.ref_bases_covered);
        auto overlaps = tree.findOverlapsWith(contig,interval).map!(x=> *(cast(BasicInterval *) x.interval)).array;
        if(!overlaps.length){
            static if(!eject) bamw.write(&rec);
            continue;
        }
        auto refRegions = overlaps.map!(x => BasicInterval(x.start - rec.pos, x.end - rec.pos)).array;
        if(refRegions[0].start < 0) refRegions[0].start = 0;
        if(refRegions[$-1].end >= cigar.ref_bases_covered) refRegions[$-1].end = cigar.ref_bases_covered;

        // writeln(refRegions);
        // debug verifyRegions(refRegions, cigar.ref_bases_covered);
        
        auto readRegions = convertRposToQpos(rec.getAlignedCoordinates,refRegions);
        // writeln(readRegions);
        // debug verifyRegions(readRegions, rec.length);

        static if(invert){
            auto invertedRegions = invertIntervals(readRegions,rec.length).filter!(x => x.start != x.end).array;
            debug verifyRegions(invertedRegions, rec.length);
            nmaskRead(&rec, invertedRegions);
        }else{
            readRegions = readRegions.filter!(x => x.start != x.end).array;
            debug verifyRegions(readRegions, rec.length);
            nmaskRead(&rec, readRegions);
        }
        // writeln(invertedRegions);
        
        bamw.write(&rec);
    }
}

void verifyRegions(BasicInterval[] regions, int length){
    foreach (reg; regions)
    {
        assert(reg.start >= 0, reg.toString);
        assert(reg.end <= length, reg.toString);
        assert(reg.end - reg.start > 0, reg.toString);
    }
}

BasicInterval[] convertRposToQpos(Range)(Range alignedCoords, BasicInterval[] regions){
    BasicInterval[] readRegions = new BasicInterval[regions.length];
    foreach (i,BasicInterval region; regions)
    {
        while(alignedCoords.front.rpos < region.start){
            alignedCoords.popFront;
        }
        readRegions[i].start = alignedCoords.front.qpos;

        while(alignedCoords.front.rpos < region.end){
            if(alignedCoords.empty) break;
            alignedCoords.popFront;
        }
        readRegions[i].end = alignedCoords.front.qpos;
    }
    return readRegions;
}

BasicInterval[] invertIntervals(BasicInterval[] regions, int length){
    BasicInterval[] newRegions;
    newRegions.reserve(regions.length);
    if(regions[0].start != 0) newRegions ~= BasicInterval(0, regions[0].start);
    for(auto i = 1; i < regions.length; i++){
        newRegions ~= BasicInterval(regions[i-1].end,regions[i].start);
    }
    if(regions[$-1].end != length) newRegions ~= BasicInterval(regions[$-1].end, length);
    return newRegions;
}

void nmaskRead(SAMRecord * rec, BasicInterval[] regions){
    ubyte[] seq = bam_get_seq(rec.b)[0 .. (rec.b.core.l_qseq >> 1) + (rec.b.core.l_qseq & 1)];
    foreach(region; regions){
        nmask(seq, region.start, region.end);
    }
}