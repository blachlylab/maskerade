import std.stdio;
import std.getopt;
import maskerade.bed;
import maskerade.bam;
import dhtslib.sam;
import intervaltree;
import std.algorithm: map;
import std.range: takeExactly, generate;
import std.array: array, join;
import std.random;

string logo = "                                        /
                                        /
           NNNNNNNNNNNNNNNNNN           /
        NNNNNNNNNNNNNNNNNNNNNNNN        /
       NNNNNNNNNNNNNNNNNNNNNNNNNN       /
      NNNNN    NNNNNNNNNN    NNNNN      /
      NNN       NNNNNNNN       NNN      /
     NNNNNNnnnnNNNNNNNNNNnnnnNNNNNN     /
      NNNNNNNNNNNNNNNNNNNNNNNNNNNN      /
       NNNNNNNNN        NNNNNNNNN       /
                                        /
                                        /";

auto bases = ["\033[0;31mt\033[0m","\033[0;33mg\033[0m","\033[0;34mc\033[0m","\033[0;32ma\033[0m"];


void main(string[] args)
{
    auto background = generate!(() => uniform(0, 4)).takeExactly(logo.length).array.map!(x=> bases[x]).array.join;
    string editedLogo;
    foreach (i,c; logo)
    {
        if(c =='n') editedLogo~="\033[1;30mn\033[0m";
        else if(c =='N') editedLogo~="\033[1;30mN\033[0m";
        else if(c=='\n') editedLogo~=c;
        else editedLogo~=background[(i * 12) .. ((i * 12) + 12)];
    }
    stderr.writeln(editedLogo);
	int threads = 0;
    bool u;
    bool s;
    bool b;
	bool invert;
    bool eject;
    GetoptResult res;
    try
    {
        res = getopt(args, 
					"threads|t", "Threads for decompression/compression (will be split)", &threads,
					"bam|b", "output bam", &b, 
					"ubam|u", "output uncompressed bam", &u, 
					"sam|s", "output sam", &s,
					"invert|v", "invert n-masking",&invert,
                    "eject|e", "eject reads that don't overlap with any provided regions",&eject);
    }
    catch (GetOptException e)
    {
        stderr.writeln(e.msg);
        stderr.writeln("Run with --help/-h for detailed usage information");
    }
    if (res.helpWanted || args.length < 3)
    {
        defaultGetoptPrinter("\nmaskerade usage: ./maskerade [options] [bam/sam] [bed file of regions] (bam/sam out)\n",
                res.options);
        stderr.writeln();
        return;
    }
    auto bamr = SAMReader(args[1], (threads >> 1) + (threads & 1));
	auto tree = getIITreefromBed(args[2]);
    SAMWriter bamw;
    if (args.length > 3)
    {
        bamw = SAMWriter(args[3], bamr.header, SAMWriterTypes.DEDUCE, threads >> 1);
    }
    else
    {
        switch ((b << 2) | (u << 1) | (s))
        {
        case 0b100:
            bamw = SAMWriter(stdout, bamr.header, SAMWriterTypes.BAM, threads >> 1);
            break;
        case 0b10:
            bamw = SAMWriter(stdout, bamr.header, SAMWriterTypes.UBAM, threads >> 1);
            break;
        case 0b1:
        case 0:
            bamw = SAMWriter(stdout, bamr.header, SAMWriterTypes.SAM, threads >> 1);
            break;
        default:
            stderr.writeln("Odd combination of output flags");
            return;
        }
    }
	processReads(&bamr, &bamw, &tree,invert,eject);
	bamw.close;
}
