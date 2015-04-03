import re
import sys
from argparse import ArgumentParser

def get_new_id(m):
    if m == "cel-miR-1834":
        return "cel-miR-58b-5p"


def get_mir_protein(l, db_type):
    if db_type == "mirdb":
        cols = l.strip().split("\t")
        mir = cols[0]
        protein = cols[1].lower()
    elif db_type == "mirtarbase":
        cols = l.strip().split(",")
        mir = cols[1]
        protein = cols[3].lower()
    else:
        raise Exception("need to be mirdb or mirtarbase!")
    mir = mir.replace("mir", "miR")

    return mir, protein

if __name__ == "__main__":
    parser = ArgumentParser(description="Map miRNA-mRNA")
    parser.add_argument("--mirna", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--mapping")
    args = parser.parse_args()

    if args.mapping:
        gene_map = {}
        with open(args.mapping) as in_handle:
            for line in in_handle:
                cols = line.strip().split()
                if len(cols) > 1:
                    gene_map[cols[0].lower()] = cols[1]

    mirna_map = {}
    mimap = open(args.mirna, 'r')
    for l in mimap:
        if l.startswith("#"):
            continue
        cols = l.strip().split("\t")
        if cols[0].startswith("hsa"):
            mirna_map[cols[0]] = 0
    mimap.close()

    newfile = {}
    with open(args.targets, 'r') as targets:
        for l in targets:
                mir, protein = get_mir_protein(l, args.db)
                if not mir.startswith("hsa"):
                    continue
                if protein in gene_map:
                    protein = gene_map[protein]
                if mir in mirna_map:
                    print "%s\t%s" % (mir, protein)
                elif mir+"-5p" in mirna_map:
                    print "%s-5p\t%s" % (mir, protein)
                elif mir+"a-5p" in mirna_map:
                    print "%sa-5p\t%s" % (mir, protein)
                elif mir+"a" in mirna_map:
                    print "%sa\t%s" % (mir, protein)
                elif mir.replace(".1", "-5p") in mirna_map:
                    print "%sa\t%s" % (mir.replace(".1", "5p"), protein)
                elif mir.replace(".2", "-5p") in mirna_map:
                    print "%sa\t%s" % (mir.replace(".2", "5p"), protein)
                elif re.sub("-[1-9]-", "-", mir) in mirna_map:
                    print "%s\t%s" % (mir, protein)
                #  elif get_new_id(mir):
                #    print "%s %s" % (get_new_id(mir), protein)
                else:
                    print "%s\t%s" % (mir, protein)
                #    print("-" + mir + "-")
                #    raise Exception("notfound!%s" % mir)
