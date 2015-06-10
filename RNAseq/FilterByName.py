import argparse
 
 
def main(parser):
 
    alreadySeen = set([])
 
    args = parser.parse_args()
    ifile = args.input
    ofile = args.output
    verbose = args.verbose
 
    print("Filtering transcripts from file {}, writing to {}".format(ifile.name, ofile.name))
    
    for l in ifile:
        if l.startswith('>'):
            toks = l.split()
            tname = toks[0][1:]
            if tname in alreadySeen:
                if verbose:
                    print("saw duplicate {}".format(tname))
            else:
                ofile.write(l)
                alreadySeen.add(tname)
        else:
            ofile.write(l)
 
    ofile.close()
    ifile.close()
 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter transcripts with duplicate names.')
    parser.add_argument('--input', type=argparse.FileType('r'))
    parser.add_argument('--output', type=argparse.FileType('w'))
    parser.add_argument('--verbose', action='store_true')
    main(parser)
