module maskerade.nmask;
import htslib.hts;

pragma(inline, true)
void nmaskSimple(ubyte[] seq, bool oddStart = false, bool oddEnd = false){
    ubyte MASK = 0xFF;
    ubyte SMASK = oddStart ? 0xF : MASK;
    ubyte EMASK = oddEnd ? 0xF0 : MASK;
    seq[0] |= SMASK;
    for(auto i = 1; i < seq.length - 1; i++){
        seq[i] |= MASK;
    }
    seq[$-1] |= EMASK;
}

pragma(inline, true)
void nmask(ubyte[] seq, long start, long end){
    assert(start >= 0);
    assert((end >> 1) + (end & 1) <= seq.length);
    nmaskSimple(seq[start >> 1 .. (end >> 1) + (end & 1)],start & 1, end & 1);
}

private ubyte[] nibbleEncoding(const(char)[] seq){
    ubyte[] res = new ubyte[(seq.length >> 1) + (seq.length & 1)];
    for(auto i=0; i < (seq.length >> 1);i++){
        res[i] = cast(ubyte)(seq_nt16_table[seq[i << 1]] << 4);
        res[i] |= cast(ubyte)(seq_nt16_table[seq[(i << 1) + 1]]);
    }
    if(seq.length & 1){
        res[$-1] = cast(ubyte)(seq_nt16_table[seq[$-1]] << 4);
    }
    return res;
}

private char[] nibbleDecoding(ubyte[] en, ulong seqlen){
    assert(en.length == ((seqlen >> 1) + (seqlen & 1)));
    char[] ret;
    ret.length = seqlen;
    for(auto i=0; i < (seqlen >> 1);i++){
        ret[i << 1] = seq_nt16_str[en[i] >> 4];
        ret[(i << 1) + 1] = seq_nt16_str[en[i] & 0x0F];
    }
    if(seqlen & 1){
        ret[$-1] = seq_nt16_str[en[$-1] >> 4];
    }
    return ret;
}

unittest
{
    import std.stdio;
    string seq = "GATCGATCGATCGCTACG";
    auto en = nibbleEncoding(seq);
    nmask(en, 1, 3);
    assert(nibbleDecoding(en,seq.length) == "GNNCGATCGATCGCTACG");
    nmask(en, 1, 4);
    assert(nibbleDecoding(en,seq.length) == "GNNNGATCGATCGCTACG");
    nmask(en, 0, 4);
    assert(nibbleDecoding(en,seq.length) == "NNNNGATCGATCGCTACG");
    nmask(en, 5, 7);
    assert(nibbleDecoding(en,seq.length) == "NNNNGNNCGATCGCTACG");
}
