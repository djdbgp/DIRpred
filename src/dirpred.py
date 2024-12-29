import sys
import utils as dc
import argparse
import os
import logging
import datetime
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
import matplotlib
import json

# good luck~

def read_input_file(args):
    """
    this function is to parse and read the input file.
    everything will be stored in parser arguments

    parse param added:
    args.ligands    -   dictionary of ligand class
    args.msa        -   multiple sequence alignment inside MSA class
    args.msta       -   multiple structure alignment inside MSA class
    args.receptor   -   receptor class
    :param args:
    :return:
    """

    # read standard input file
    ligands = {}
    ref = args.reference

    with open(args.input_file) as infile:
        for line in infile:
            if line[0] != "#":
                line = line.strip().split()
                type = line[0]

                # new ligand
                if type == "LIGAND":
                    p, name, pdb = line[1:]
                    this_lig = dc.Ligand(name, pdb=pdb)
                    this_lig.prot_msa = dc.create_msa(open(p))
                    assert name not in ligands, "I am not overwriting ligand {}".format(name)
                    ligands[name] = this_lig
                # new MSA
                if type == "MSA":
                    p = line[1]
                    msa = dc.MSA("MSA", reference=ref)
                    msa.msa = dc.create_msa(open(p))
                    args.msa = msa
                # new MSTA
                if type == "MSTA":
                    p = line[1]
                    msta = dc.MSA("MSTA", reference=ref)
                    msta.msa = dc.create_msa(open(p))
                    args.msta = msta
                # new receptor
                if type == "RECEPTOR":
                    name, p = line[1:]
                    receptor = dc.Receptor(name, pdb=pdb, reference="Hsap")
                    receptor.prot_msa = dc.create_msa(open(p))
                    args.receptor = receptor

        # add ligands
        args.ligands = ligands


def dir_prediction(args):
    """
    Divergence Inducing Residue Prediction
    Basically like cross-conservation but averaging over all paralogs and correcting for coevolution

    :param args:
    :return:
    """
    rp = args.reference          # reference ligand name
    prot = args.ligands[rp]             # reference ligand
    seq = prot.get_ref_seq()            # reference ligand sequence
    msta = args.msta                    # multiple structure alignment
    msa = args.msa                      # multiple sequence alignment
    contest = args.conservation_test    # conservation test type
    algtest = args.alignment_test       # conservation for ligand alignment type

    # FOR NOW HARDCODE
    args.msa_cols = ["EVO cons MSA", "MSA", "COEVO receptor MSA", "COEVO ligands MSA"]
    args.msta_cols = ["EVO cons MSTA","MSTA", "COEVO receptor MSTA","COEVO ligands MSTA"]
    # data preparation
    cs = {}

    # ligands alignment scorings
    cs["MSA"] = [round(x, 2) for x in msa.get_cons_scores(algtest)]
    cs["MSTA"] = [round(x, 2) for x in msta.get_cons_scores(algtest)]
    # now I need one evo cons from MSA alignment average and one from MSTA alignment average
    relative_conservation = ligand_table(args)
    cs["EVO cons MSA"] = average_score(relative_conservation,"MSA")
    cs["EVO cons MSTA"] = average_score(relative_conservation,"MSTA")

    # now let's do the coevolution score
    rel_intra_coevol, rel_inter_coevol = coevol_test(args)
    cs["COEVO receptor MSA"] = average_score(rel_inter_coevol,"MSA")
    cs["COEVO receptor MSTA"] = average_score(rel_inter_coevol,"MSTA")
    cs["COEVO ligands MSA"] = average_score(rel_intra_coevol,"MSA")
    cs["COEVO ligands MSTA"] = average_score(rel_intra_coevol,"MSTA")

    # now make it DataFrame
    cs = pd.DataFrame(cs)
    # insert relative positions
    cs.insert(0, 'POS', [str(x) + seq[x - 1] for x in range(1, 1 + len(cs))])

    # get the values of abcd:
    abcd = [1/int(x) for x in args.abcd.split(",")]
    # average MSA/MSTA then distance with
    avg_dirp, part_scores = calculate_dirp_score(cs,args.msa_cols,args.msta_cols, abcd)
    cs["avg DIR score"] = avg_dirp
    cs["partial scores"] = part_scores

    makeplots(cs, args)
    # plotting


def makeplots(cs, args):

    relabel_cols = {
      "EVO cons MSA": "Evolutionary conservation",
        "EVO cons MSTA": "Evolutionary conservation",
        "MSA": "Paralogs conservation",
        "MSTA": "Paralogs conservation",
        "COEVO receptor MSA":"Ligand-Receptor co-evolution",
        "COEVO receptor MSTA":"Ligand-Receptor co-evolution",
        "COEVO ligands MSA":"Ligands internal co-evolution",
        "COEVO ligands MSTA":"Ligands internal co-evolution"
    }
    # plots folder
    plot_output = args.output_path + "/plots/"
    if not os.path.exists(plot_output):
        os.makedirs(plot_output)
    logging.warning("running plotting module")

    # palette
    sb.set_palette(sb.color_palette("Dark2"))

    # distribution plots
    df_dist = cs[args.msa_cols + args.msta_cols].melt(var_name='rows')
    df_dist['type'] = ["MSA" if "MSA" in x else "MSTA" for x in df_dist['rows']]
    df_dist['rows'] = [" ".join(x.split(" ")[:-1]) if len(x) > 4 else "Ligands alignment" for x in df_dist['rows']]

    # rename
    df_dist = df_dist.rename(columns=relabel_cols)

    # plot distributions
    g = sb.FacetGrid(df_dist, row='rows', col='type', hue='rows')
    g = (g.map(sb.histplot, 'value')).set_titles("{col_name} | {row_name}")
    plt.savefig("{}/scores_distribution.png".format(plot_output), dpi=400)
    plt.close()

    # stacked box plot
    matplotlib.rcParams.update({'font.size': 10})
    # format data
    ilab = [args.msa_cols, args.msta_cols]
    df_stack = cs[['POS', 'partial scores']].copy()
    for laidx in range(2):
        for idx in range(4):
            this_col = []
            for x in df_stack['partial scores']:
                this_col.append(x[laidx][idx])
            df_stack["{}".format(ilab[laidx][idx])] = this_col
    del (df_stack['partial scores'])
    #rename
    df_msa = df_stack[args.msa_cols].rename(columns={k:v for k,v in relabel_cols.items() if k in args.msa_cols})
    df_msta = df_stack[args.msta_cols].rename(columns={k:v for k,v in relabel_cols.items() if k in args.msta_cols})
    # figure out how to do this
    fig, axes = plt.subplots(2, 1, figsize=(12, 24))
    df_msa.iloc[::-1].plot.barh(stacked=True,ax=axes[0])
    axes[0].set_yticklabels(cs['POS'][::-1])
    axes[0].set_title('MSA',fontsize=24)
    df_msta.iloc[::-1].plot.barh(stacked=True,ax=axes[1],legend=False)
    axes[1].set_yticklabels(cs['POS'][::-1])
    axes[1].set_title('MSTA',fontsize=24)

    # handle legends
    patches, labels = axes[0].get_legend_handles_labels()
    axes[0].legend(patches,labels,bbox_to_anchor=(0.45,-0.08),ncol=2,loc="center", prop={'size': 18})
    plt.savefig("{}/combined_score.png".format(plot_output), dpi=400)
    plt.close()

    # same but sorted
    fig, axes = plt.subplots(2, 1, figsize=(12, 24))
    msa_sort = df_stack[args.msa_cols + ["POS"]].copy()
    msa_sort["sum"] = msa_sort[args.msa_cols].sum(axis=1)
    msa_sort = msa_sort.sort_values(by=['sum'])
    del(msa_sort['sum'])
    msa_lab = msa_sort.pop('POS')
    msa_sort = msa_sort.rename(columns=relabel_cols)
    msa_sort.plot.barh(stacked=True,ax=axes[0])
    axes[0].set_yticklabels(msa_lab)
    axes[0].set_title('MSA',fontsize=32)
    #msta
    msta_sort = df_stack[args.msta_cols + ["POS"]].copy()
    msta_sort["sum"] = msta_sort[args.msta_cols].sum(axis=1)
    msta_sort = msta_sort.sort_values(by=['sum'])
    del(msta_sort['sum'])
    msta_lab = msta_sort.pop('POS')
    msta_sort.plot.barh(stacked=True,ax=axes[1],legend=False)
    axes[1].set_yticklabels(msta_lab)
    axes[1].set_title('MSTA',fontsize=32)
    # handle legends
    patches, labels = axes[0].get_legend_handles_labels()
    axes[0].legend(patches,labels,bbox_to_anchor=(0.45,-0.08),ncol=2,loc="center", prop={'size': 18})

    plt.savefig("{}/combined_score_sorted.png".format(plot_output), dpi=400)
    plt.close()

    # save top 12
    matplotlib.rcParams.update({'font.size': 18})
    fig, axes = plt.subplots(2, 1, figsize=(12, 24))
    msa_sort.tail(12).plot.barh(stacked=True,ax=axes[0])
    axes[0].set_yticklabels(msa_lab[-12:])
    axes[0].set_title('MSA',fontsize=32)
    msta_sort.tail(12).plot.barh(stacked=True,ax=axes[1],legend=False)
    axes[1].set_yticklabels(msta_lab[-12:])
    axes[1].set_title('MSTA',fontsize=32)
    # handle legends
    patches, labels = axes[0].get_legend_handles_labels()
    axes[0].legend(patches,labels,bbox_to_anchor=(0.45,-0.09),ncol=2,loc="center", prop={'size': 18})

    plt.savefig("{}/combined_score_sorted_top12.png".format(plot_output), dpi=400)
    plt.close()


    # save text output also but add average
    #index plus one
    cs.index = cs.index + 1
    with open("{}/DIRpred.csv".format(args.output_path), "w") as outfile:
        cs.to_csv(outfile)
    # save also sorted version
    sorted_cs = cs.sort_values("avg DIR score", ascending=False)
    sorted_cs.index = [i for i in range(1,len(sorted_cs.index)+1)]
    with open("{}/DIRpred_sorted.csv".format(args.output_path), "w") as outfile:
        sorted_cs.to_csv(outfile)


def calculate_dirp_score(df, msa_ids, msta_ids, abcd=(0.25,0.25,0.25,0.25)):
    """sum method, v1.0
    the 4 scores are
    I: avg evolutionary conservation (positive)
    II: ligands alignment conservation (negative)
    III: avg ligand-receptor max coevolution (positive)
    IV: avg ligand-ligand max coevolution (negative)

    the 4 scores will be weighted by a parameter each, by default 0.25, and summed to compose the final score

    returns:
    average MSA,MSTA dirp score array = a(I) + b(1-II) + c(III) + d(1-IV)
    partial scores array = same as before but comma instead of sum, and (MSAps,MSTAps)
    """

    a,b,c,d = abcd
    MSA_score = []
    MSTA_score = []
    partial_score = []

    for i,r in df.iterrows():
        # msa
        I,II,III,IV = list(r[msa_ids])
        msa_partial_score = a * I, b * (1-II), c * III, d * (1 - IV)
        MSA_score.append(sum(msa_partial_score))
        # msta
        I,II,III,IV = list(r[msta_ids])
        msta_partial_score = a * I, b * (1-II), c * III, d * (1 - IV)
        MSTA_score.append(sum(msta_partial_score))
        partial_score.append((msa_partial_score, msta_partial_score))

    return np.mean([MSA_score, MSTA_score], axis=0).round(2), partial_score


def average_score(d,alg_type):
    """returns the average value from the dictionary of relative scores"""
    l = []
    for k in d:
        if alg_type in k:
            a = [np.nan if x == "-" else float(x) for x in d[k]]
            l.append(a)
    mean_out = np.nanmean(l,axis=0).round(2)
    return mean_out


def ligand_table(args):
    """
    realigns the paralogs to a reference sequence and save the results along with the reference-aligned conservation, or sort of...
    :param args:
    :return:
    """
    rp = args.reference          # reference ligand name
    prot = args.ligands[rp]             # reference ligand
    seq = prot.get_ref_seq()            # reference ligand sequence
    msta = args.msta                    # multiple structure alignment
    msa = args.msa                      # multiple sequence alignment
    contest = args.conservation_test    # conservation test type
    algtest = args.alignment_test       # ligand alignment test type

    # data preparation
    cs = {}
    rel_cons = {}  # relative conservation is stored here

    # add referenced positions
    for m in [msa, msta]:
        for n in m.get_names():
            logging.warning("processing {}".format(n))
            assert n in args.ligands, "{} not found in ligands list: {}".format(n,args.ligands.keys())
            this_ligand = args.ligands[n]
            this_ligand_score = this_ligand.get_cons_scores(contest)
            this_name = "{}_{}".format(m.name, n)
            cs[this_name] = m.get_referenced_positions(n)
            rel_cons[this_name] = m.get_referenced_scores(n, this_ligand, this_ligand_score)

    # add also cons score
    cs["MSA"] = [round(x, 2) + 0.04 for x in msa.get_cons_scores(algtest)]
    cs["MSTA"] = [round(x, 2) - 0.04 for x in msta.get_cons_scores(algtest)]
    cs[prot.name] = prot.get_cons_scores(contest)
    for n in rel_cons:
        cs[n+"_cons_{}".format(contest)] = rel_cons[n]

    cs = pd.DataFrame(cs)
    # add egf pos
    cs.insert(0, 'POS', [str(x) + seq[x - 1] for x in range(1, 1 + len(cs))])

    with open("{}/ligand_table.csv".format(args.output_path), "w") as outfile:
        cs.to_csv(outfile)

    return rel_cons


def coevol_test(args):
    """
    has to return two dictionaries, intra and inter coevolution
    as key msa_ligand, as values, an array based on reference with the max coevolution score
    :param args:
    :return:
    """
    # preparation
    rp = args.reference          # reference ligand name
    prot = args.ligands[rp]             # reference ligand
    seq = prot.get_ref_seq()            # reference sequence
    receptor = args.receptor            # receptor MSA
    msta = args.msta  # multiple structure alignment
    msa = args.msa  # multiple sequence alignment
    coevo_test = args.coevolution_test
    coevol_output = args.output_path + "/coevolution/"
    if not os.path.exists(coevol_output):
        os.makedirs(coevol_output)
    logging.warning("running coevolution module")

    # receptor intra coevol
    _ = dc.intra_coevol(receptor, "{}/{}".format(coevol_output, receptor.name), coevo_test)
    # every ligand coevol
    # now try to combine all ligands together using msa and msta
    # all ligands
    intra_rc,inter_rc = {},{}
    report_intra, report_inter = {},{}
    for m in [msa, msta]:
        lig_rc = dc.all_ligands_coevol(args, m, "{}/all_ligands".format(coevol_output), coevo_test)
        this_name = "{}_all_ligands".format(m.name)
        intra_rc[this_name] = lig_rc
        report_intra[this_name] = lig_rc

        for n in m.get_names():
            logging.warning("processing {}".format(n))
            assert n in args.ligands, "{} not found in ligands list: {}".format(n, args.ligands.keys())
            this_ligand = args.ligands[n]
            this_ligand_inter_score = dc.inter_coevol(this_ligand, receptor, "{}/{}u{}intercoevo".format(coevol_output, this_ligand.name, receptor.name), coevo_test)
            this_name = "{}_{}".format(m.name, n)
            report_inter[this_name] = m.get_referenced_positions(n)
            inter_rc[this_name] = m.get_referenced_scores(n, this_ligand, this_ligand_inter_score)

    report_intra["MSA"] = dc.intra_coevol(msa,"{}/{}".format(coevol_output, msa.name), coevo_test)
    report_intra["MSTA"] = dc.intra_coevol(msta, "{}/{}".format(coevol_output, msta.name), coevo_test)
#                                                 args.coevolution_test)
    report_inter[prot.name] = dc.inter_coevol(prot, receptor, "{}/inter_{}".format(coevol_output, prot.name), coevo_test)
    for n in intra_rc:
        report_intra[n+"_coevol_{}".format(coevo_test)] = intra_rc[n]
    for n in inter_rc:
        report_inter[n+"_coevol_{}".format(coevo_test)] = inter_rc[n]
    # make df
    report_intra_df = pd.DataFrame(report_intra)
    report_inter_df = pd.DataFrame(report_inter)
    # add egf pos
    report_intra_df.insert(0, 'POS', [str(x) + seq[x - 1] for x in range(1, 1 + len(report_intra_df))])
    report_inter_df.insert(0, 'POS', [str(x) + seq[x - 1] for x in range(1, 1 + len(report_inter_df))])

    # save reports
    with open("{}/intra_coevol_report.csv".format(coevol_output), "w") as outfile:
        report_intra_df.to_csv(outfile)
    with open("{}/inter_coevol_report.csv".format(coevol_output), "w") as outfile:
        report_inter_df.to_csv(outfile)

    return intra_rc, inter_rc


def main(argv):

    # PARSER
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@', description="parser for handling input files of cons-cons pipeline")

    parser.add_argument("--input-file", required=True, help="input file contains all ligands with their MSA location"
                                                        "the MStA path and its format is very important, check example")
    parser.add_argument("--output-path", required=True, help="path to output folder, will be created if not ex")
    parser.add_argument("--test-type", default="DIRP", help="""The type of test:
     DIRP - average multiple concon, or Divergence Inducing Residue Prediction
     """)
    parser.add_argument("--reference", default="", help="the reference protein for DIRpred test.")
    parser.add_argument("--conservation-test", default="id", help="""type of conservation measure:"
     id - identity
     blos - blosum62
     jsdw - jensen shannon divergence with window (as in concavity)""")
    parser.add_argument("--alignment-test", default="", help="type of conservation measure, used for alignment. By default, like conservation_test")
    parser.add_argument("--coevolution-test", default="MI",
                        help="coevolution measure. By default, MI: mutual information. can be MIp, APC corrected MI")
    parser.add_argument("--abcd", default="4,4,4,4", help="the individual contributions for each score, as in 1 over x, by defaultis: 4,4,4,4. meaning 1 over 4, four times.")
    args = parser.parse_args()

    # CREATE LOG DETAILS
    if not os.path.exists("logs"):
        os.mkdir("logs")
    logging.basicConfig(level=logging.DEBUG, filename="logs/{}.txt".format(str(datetime.datetime.now()).replace(" ", "_")))
    logging.warning("LOGFILE INITIATED\nPIPELINE STARTED\nOk, I am back to life and ready to conquer the wo... oh.. damn..\n")

    # OUTPUT FOLDER CREATION
    o = args.output_path
    if not os.path.exists(o):
        os.makedirs(o)
        with open('{}/args.txt'.format(o), 'w') as f:
            json.dump(args.__dict__, f, indent=2)
        logging.warning("~~~\nCREATING OUTPUT FOLDER at {}\n~~~".format(o))

    # CONSERVATION check
    if not args.alignment_test:
        args.alignment_test = args.conservation_test

    # INPUT FILE READER
    read_input_file(args)

    # REPORT ARGS
    logging.warning("~~~~~pipeline has the following args:\n{}\n~~~~~".format(args))

    # DIRP
    if args.test_type == "DIRP":
        dir_prediction(args)

    if args.test_type == "COEVOL":
        coevol_test(args)

    if args.test_type == "LIGTAB":
        _ = ligand_table(args)

if __name__ == '__main__':
    main(sys.argv)

