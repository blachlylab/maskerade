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

void processReads(SAMReader * bamr, SAMWriter * bamw, IITree!BasicInterval * tree){
    auto contigs = bamr.target_names;
    foreach(SAMRecord rec; bamr.all_records){
        if(rec.isSecondary || !rec.isMapped) continue;
        auto contig = contigs[rec.tid];
        auto cigar = Cigar(rec.cigar.ops.dup);
        auto interval = BasicInterval(rec.pos, rec.pos + cigar.ref_bases_covered);
        auto overlaps = tree.findOverlapsWith(contig,interval).map!(x=> *(cast(BasicInterval *) x.interval)).array;
        if(!overlaps.length) continue;
        auto refRegions = overlaps.map!(x => BasicInterval(x.start - rec.pos, x.end - rec.pos)).array;
        if(refRegions[0].start < 0) refRegions[0].start = 0;
        if(refRegions[$-1].end >= cigar.ref_bases_covered) refRegions[$-1].end = cigar.ref_bases_covered;

        // writeln(refRegions);
        // debug verifyRegions(refRegions, cigar.ref_bases_covered);
        
        auto readRegions = convertRposToQpos(cigar,refRegions);
        // writeln(readRegions);
        // debug verifyRegions(readRegions, rec.length);

        auto invertedRegions = invertIntervals(readRegions,rec.length).filter!(x => x.start != x.end).array;
        // writeln(invertedRegions);
        debug verifyRegions(invertedRegions, rec.length);

        nmaskRead(&rec, invertedRegions);
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

BasicInterval[] convertRposToQpos(Cigar cigar, BasicInterval[] regions){
    auto ops = cigar.ops;
    int qpos, rpos;
    BasicInterval[] readRegions = new BasicInterval[regions.length];
    if(ops[0].op == Ops.HARD_CLIP) ops = ops[1..$];
    if(ops[0].op == Ops.SOFT_CLIP) qpos += ops[0].length, ops = ops[1..$];
    if(ops[$-1].op == Ops.HARD_CLIP) ops = ops[0..$-1];
    foreach (i,BasicInterval region; regions)
    {
        CigarOp last_op;
        bool iterated = false;
        while(rpos < region.start){
            iterated = true;
            auto refConsuming = ops[0].is_reference_consuming;
            auto queryConsuming = ops[0].is_query_consuming;
            if(refConsuming){
                rpos += ops[0].length;
            }
            if(queryConsuming){
                qpos += ops[0].length;
            }
            last_op = ops[0];
            ops = ops[1..$];
        }
        if(iterated && last_op.is_query_consuming) readRegions[i].start = qpos - (rpos - region.start);
        else readRegions[i].start = qpos;
        // writeln(readRegions[i].start," ",rpos," ",qpos," ",region.end," ",last_op.is_reference_consuming," ",iterated);
        last_op = last_op.init;
        iterated = false;

        while(rpos < region.end){
            iterated = true;
            auto refConsuming = ops[0].is_reference_consuming;
            auto queryConsuming = ops[0].is_query_consuming;
            if(refConsuming){
                rpos += ops[0].length;
            }
            if(queryConsuming){
                qpos += ops[0].length;
            }
            last_op = ops[0];
            ops = ops[1..$];
        }
        // writeln(readRegions[i].start," ",rpos," ",qpos," ",region.end," ",last_op.is_reference_consuming," ",iterated);
        if(iterated && last_op.is_query_consuming) readRegions[i].end = qpos - (rpos - region.end);
        else readRegions[i].end = qpos;
        // writeln(readRegions[i]);
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